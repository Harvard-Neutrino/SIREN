"""Tests for siren.download -- file download, integrity verification, and caching."""
from __future__ import annotations

import hashlib
import http.server
import os
import threading
import zipfile
from pathlib import Path

import pytest

# Import from the source tree so tests work without installing
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import download as dl


# ======================================================================
# Fixtures
# ======================================================================

@pytest.fixture()
def tmp(tmp_path):
    """Short alias for tmp_path."""
    return tmp_path


@pytest.fixture()
def sample_file(tmp):
    """Create a small sample file and return (path, sha256)."""
    path = tmp / "sample.dat"
    content = b"hello siren\n" * 100
    path.write_bytes(content)
    sha = hashlib.sha256(content).hexdigest()
    return path, sha


@pytest.fixture()
def sample_zip(tmp):
    """Create a zip with a directory structure mimicking Zenodo archives."""
    zip_path = tmp / "archive.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("data/alpha/file1.txt", "alpha-1")
        zf.writestr("data/alpha/file2.txt", "alpha-2")
        zf.writestr("data/beta/file1.txt", "beta-1")
        zf.writestr("data/beta/sub/deep.txt", "beta-deep")
    return zip_path


class _LocalServer:
    """Minimal HTTP server that serves files from a directory."""

    def __init__(self, serve_dir: Path):
        self.serve_dir = serve_dir
        handler = http.server.SimpleHTTPRequestHandler
        self.httpd = http.server.HTTPServer(
            ("127.0.0.1", 0),
            lambda *args, **kwargs: handler(*args, directory=str(serve_dir), **kwargs))
        self.port = self.httpd.server_address[1]
        self.url = f"http://127.0.0.1:{self.port}"
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.thread.start()

    def stop(self):
        self.httpd.shutdown()


@pytest.fixture()
def local_server(tmp):
    """Start a local HTTP server serving from tmp, return (server, serve_dir)."""
    srv = _LocalServer(tmp)
    yield srv, tmp
    srv.stop()


# ======================================================================
# download_file
# ======================================================================

class TestDownloadFile:

    def test_basic_download(self, local_server, tmp):
        srv, serve_dir = local_server
        (serve_dir / "hello.txt").write_text("world")
        dest = tmp / "out" / "hello.txt"
        dl.download_file(f"{srv.url}/hello.txt", str(dest), show_progress=False)
        assert dest.read_text() == "world"

    def test_creates_parent_dirs(self, local_server, tmp):
        srv, serve_dir = local_server
        (serve_dir / "a.txt").write_text("data")
        dest = tmp / "deep" / "nested" / "dir" / "a.txt"
        dl.download_file(f"{srv.url}/a.txt", str(dest), show_progress=False)
        assert dest.exists()

    def test_sha256_pass(self, local_server, tmp):
        srv, serve_dir = local_server
        content = b"verify me"
        (serve_dir / "v.bin").write_bytes(content)
        sha = hashlib.sha256(content).hexdigest()
        dest = tmp / "v.bin"
        dl.download_file(f"{srv.url}/v.bin", str(dest), sha256=sha,
                         show_progress=False)
        assert dest.read_bytes() == content

    def test_sha256_mismatch(self, local_server, tmp):
        srv, serve_dir = local_server
        (serve_dir / "bad.bin").write_bytes(b"actual content")
        dest = tmp / "output" / "bad.bin"
        with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
            dl.download_file(f"{srv.url}/bad.bin", str(dest),
                             sha256="0" * 64, show_progress=False)
        # Failed file should be cleaned up
        assert not dest.exists()

    def test_atomic_write_on_failure(self, tmp):
        dest = tmp / "never.txt"
        with pytest.raises(Exception):
            dl.download_file("http://127.0.0.1:1/nonexistent", str(dest),
                             show_progress=False)
        assert not dest.exists()
        assert not (tmp / "never.txt.tmp").exists()

    def test_404_raises(self, local_server, tmp):
        srv, serve_dir = local_server
        dest = tmp / "missing.txt"
        with pytest.raises(Exception):
            dl.download_file(f"{srv.url}/no_such_file.txt", str(dest),
                             show_progress=False)


# ======================================================================
# ensure_files
# ======================================================================

class TestEnsureFiles:

    def test_skips_existing(self, tmp):
        existing = tmp / "already.txt"
        existing.write_text("present")
        dl.ensure_files([{"path": str(existing), "url": "http://should.not.be.called"}])
        assert existing.read_text() == "present"

    def test_downloads_missing(self, local_server, tmp):
        srv, serve_dir = local_server
        (serve_dir / "needed.txt").write_text("fetched")
        dest = tmp / "out" / "needed.txt"
        dl.ensure_files([{
            "path": str(dest),
            "url": f"{srv.url}/needed.txt",
        }])
        assert dest.read_text() == "fetched"

    def test_missing_no_url_raises(self, tmp):
        with pytest.raises(FileNotFoundError, match="no download URL"):
            dl.ensure_files([{"path": str(tmp / "ghost.txt")}])

    def test_mixed_existing_and_missing(self, local_server, tmp):
        srv, serve_dir = local_server
        existing = tmp / "old.txt"
        existing.write_text("old")
        (serve_dir / "new.txt").write_text("new")
        dl.ensure_files([
            {"path": str(existing)},
            {"path": str(tmp / "new.txt"), "url": f"{srv.url}/new.txt"},
        ])
        assert existing.read_text() == "old"
        assert (tmp / "new.txt").read_text() == "new"


# ======================================================================
# ensure_zenodo_archive -- with a local mock
# ======================================================================

class TestEnsureZenodoArchive:
    """Test zip caching and prefix extraction using a local server."""

    def _serve_zip(self, local_server, sample_zip):
        """Copy the sample zip into the serve directory."""
        srv, serve_dir = local_server
        import shutil
        served = serve_dir / "test_archive.zip"
        shutil.copy2(sample_zip, served)
        return srv, served

    def test_full_extract(self, local_server, sample_zip, tmp, monkeypatch):
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "extract_dest"
        dest.mkdir()

        # Mock zenodo_file_url to return our local URL
        monkeypatch.setattr(dl, "zenodo_file_url",
                            lambda *a, **kw: f"{srv.url}/test_archive.zip")

        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest))
        assert (dest / "data" / "alpha" / "file1.txt").read_text() == "alpha-1"
        assert (dest / "data" / "beta" / "file1.txt").read_text() == "beta-1"

    def test_prefix_extract(self, local_server, sample_zip, tmp, monkeypatch):
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "prefix_dest"
        dest.mkdir()

        monkeypatch.setattr(dl, "zenodo_file_url",
                            lambda *a, **kw: f"{srv.url}/test_archive.zip")

        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                 prefix="data/alpha")
        assert (dest / "data" / "alpha" / "file1.txt").exists()
        assert not (dest / "data" / "beta").exists()

    def test_skip_when_sentinel_present(self, local_server, sample_zip, tmp, monkeypatch):
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "skip_dest"
        alpha_dir = dest / "data" / "alpha"
        alpha_dir.mkdir(parents=True)
        # Write the sentinel that ensure_zenodo_archive looks for
        (alpha_dir / dl._EXTRACTED_SENTINEL).write_text("already done")

        called = []
        real_zenodo_file_url = dl.zenodo_file_url
        def mock_url(*a, **kw):
            called.append(1)
            return real_zenodo_file_url(*a, **kw)
        monkeypatch.setattr(dl, "zenodo_file_url", mock_url)

        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                 prefix="data/alpha")
        assert not called, "should have skipped download entirely"

    def test_does_not_skip_when_dir_has_only_py_files(self, local_server, sample_zip, tmp, monkeypatch):
        """Regression: directories with .py loaders but no data should NOT be skipped."""
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "py_only_dest"
        alpha_dir = dest / "data" / "alpha"
        alpha_dir.mkdir(parents=True)
        (alpha_dir / "loader.py").write_text("# fake loader")

        monkeypatch.setattr(dl, "zenodo_file_url",
                            lambda *a, **kw: f"{srv.url}/test_archive.zip")

        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                 prefix="data/alpha")
        # Data files should have been extracted despite loader.py being present
        assert (dest / "data" / "alpha" / "file1.txt").read_text() == "alpha-1"

    def test_cache_reuse(self, local_server, sample_zip, tmp, monkeypatch):
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "cache_dest"
        dest.mkdir()

        monkeypatch.setattr(dl, "zenodo_file_url",
                            lambda *a, **kw: f"{srv.url}/test_archive.zip")

        # First call: downloads and caches
        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                 prefix="data/alpha")

        # Verify cache exists with SHA256 in filename
        cache_dir = dest / ".download_cache"
        assert cache_dir.is_dir()
        cached = list(cache_dir.glob("99999_test_archive.zip.*.zip"))
        assert len(cached) == 1
        cache_name = cached[0].name
        # Name format: {tag}.{base64_sha256}.zip
        parts = cache_name.split(".")
        assert len(parts) == 4  # tag, sha_b64, zip (with tag containing a dot)
        sha_part = parts[-2]
        assert len(sha_part) == 43  # URL-safe base64 of 32 bytes, no padding

        # Second call with different prefix: should NOT re-download
        # (remove the served file to prove it)
        served.unlink()
        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                 prefix="data/beta")
        assert (dest / "data" / "beta" / "file1.txt").read_text() == "beta-1"

    def test_bad_prefix_raises(self, local_server, sample_zip, tmp, monkeypatch):
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "bad_prefix_dest"
        dest.mkdir()

        monkeypatch.setattr(dl, "zenodo_file_url",
                            lambda *a, **kw: f"{srv.url}/test_archive.zip")

        with pytest.raises(FileNotFoundError, match="not found"):
            dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                     prefix="data/nonexistent")

    def test_sentinel_written_after_extract(self, local_server, sample_zip, tmp, monkeypatch):
        srv, served = self._serve_zip(local_server, sample_zip)
        dest = tmp / "sentinel_dest"
        dest.mkdir()

        monkeypatch.setattr(dl, "zenodo_file_url",
                            lambda *a, **kw: f"{srv.url}/test_archive.zip")

        dl.ensure_zenodo_archive("99999", "test_archive.zip", str(dest),
                                 prefix="data/alpha")
        sentinel = dest / "data" / "alpha" / dl._EXTRACTED_SENTINEL
        assert sentinel.is_file()
        assert "99999" in sentinel.read_text()

    def test_zip_slip_rejected(self, tmp):
        """Zip entries with '../' paths must be rejected."""
        malicious_zip = tmp / "evil.zip"
        with zipfile.ZipFile(malicious_zip, "w") as zf:
            zf.writestr("../../../etc/passwd", "pwned")

        dest = tmp / "victim"
        dest.mkdir()
        with zipfile.ZipFile(malicious_zip) as zf:
            with pytest.raises(ValueError, match="outside"):
                dl._safe_extract(zf, str(dest))


# ======================================================================
# _sha256_file
# ======================================================================

class TestSha256File:

    def test_known_hash(self, tmp):
        path = tmp / "known.bin"
        content = b"deterministic content"
        path.write_bytes(content)
        import base64
        expected = base64.urlsafe_b64encode(
            hashlib.sha256(content).digest()
        ).rstrip(b"=").decode("ascii")
        assert dl._sha256_file(str(path)) == expected

    def test_empty_file(self, tmp):
        path = tmp / "empty.bin"
        path.write_bytes(b"")
        import base64
        expected = base64.urlsafe_b64encode(
            hashlib.sha256(b"").digest()
        ).rstrip(b"=").decode("ascii")
        assert dl._sha256_file(str(path)) == expected
