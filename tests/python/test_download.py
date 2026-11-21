"""Tests for siren.download -- file download, integrity verification, and caching."""
from __future__ import annotations

import hashlib
import http.server
import io
import os
import tarfile
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


@pytest.fixture()
def sample_tar_xz(tmp):
    """Create an xz-compressed tar archive with a nested data tree."""
    source = tmp / "tar_source"
    files = {
        "data/alpha/file1.txt": "alpha-1",
        "data/alpha/file2.txt": "alpha-2",
        "data/beta/file1.txt": "beta-1",
        "data/beta/sub/deep.txt": "beta-deep",
    }
    for rel_path, content in files.items():
        path = source / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    archive = tmp / "archive.tar.xz"
    with tarfile.open(archive, "w:xz") as tf:
        for path in sorted(source.rglob("*")):
            if path.is_file():
                tf.add(path, arcname=path.relative_to(source))
    return archive


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
# ensure_tar_archive
# ======================================================================

class TestEnsureTarArchive:

    def test_downloads_verifies_and_extracts(self, local_server, sample_tar_xz,
                                             tmp):
        srv, _ = local_server
        sha = hashlib.sha256(sample_tar_xz.read_bytes()).hexdigest()
        dest = tmp / "generic_extract"

        dl.ensure_tar_archive(
            f"{srv.url}/{sample_tar_xz.name}", sample_tar_xz.name,
            str(dest), sha)

        assert (dest / "data" / "alpha" / "file1.txt").read_text() == "alpha-1"
        cached = list((dest / ".download_cache").glob(f"{sha}.*"))
        assert len(cached) == 1
        sentinels = list(dest.glob(f".*{dl._ARCHIVE_EXTRACTED_SENTINEL}"))
        assert len(sentinels) == 1
        assert f"sha256:{sha}" in sentinels[0].read_text()

    def test_sentinel_skips_second_download(self, local_server, sample_tar_xz,
                                            tmp, monkeypatch):
        srv, _ = local_server
        sha = hashlib.sha256(sample_tar_xz.read_bytes()).hexdigest()
        dest = tmp / "generic_cached"
        url = f"{srv.url}/{sample_tar_xz.name}"

        dl.ensure_tar_archive(url, sample_tar_xz.name, str(dest), sha)
        monkeypatch.setattr(
            dl, "download_file",
            lambda *a, **kw: pytest.fail("archive was downloaded twice"))
        dl.ensure_tar_archive(url, sample_tar_xz.name, str(dest), sha)

    def test_bad_sha_does_not_extract(self, local_server, sample_tar_xz, tmp):
        srv, _ = local_server
        dest = tmp / "generic_bad_sha"

        with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
            dl.ensure_tar_archive(
                f"{srv.url}/{sample_tar_xz.name}", sample_tar_xz.name,
                str(dest), "0" * 64)

        assert not (dest / "data").exists()
        assert not list(dest.glob(f".*{dl._ARCHIVE_EXTRACTED_SENTINEL}"))

    def test_path_traversal_rejected_before_extract(self, local_server, tmp):
        srv, serve_dir = local_server
        archive = serve_dir / "generic_evil.tar.xz"
        content = b"pwned"
        with tarfile.open(archive, "w:xz") as tf:
            member = tarfile.TarInfo("../../../escaped.txt")
            member.size = len(content)
            tf.addfile(member, io.BytesIO(content))
        sha = hashlib.sha256(archive.read_bytes()).hexdigest()
        dest = tmp / "generic_tar_slip"

        with pytest.raises(ValueError, match="outside"):
            dl.ensure_tar_archive(
                f"{srv.url}/{archive.name}", archive.name, str(dest), sha)

        assert not (tmp.parent / "escaped.txt").exists()

    def test_links_are_rejected(self, local_server, tmp):
        srv, serve_dir = local_server
        archive = serve_dir / "generic_link.tar.xz"
        with tarfile.open(archive, "w:xz") as tf:
            member = tarfile.TarInfo("link")
            member.type = tarfile.SYMTYPE
            member.linkname = "target"
            tf.addfile(member)
        sha = hashlib.sha256(archive.read_bytes()).hexdigest()

        with pytest.raises(ValueError, match="links are not allowed"):
            dl.ensure_tar_archive(
                f"{srv.url}/{archive.name}", archive.name,
                str(tmp / "generic_link"), sha)


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


class TestResolveDataPath:
    """resolve_data_path prefers an existing install-dir copy, else the download dir."""

    def test_prefers_install_dir_when_present(self, tmp):
        install = tmp / "install"
        download = tmp / "download"
        install.mkdir()
        download.mkdir()
        (install / "f.dat").write_bytes(b"x")
        # Present in install_dir -> returned directly (no redundant download)
        assert dl.resolve_data_path(str(install), str(download), "f.dat") == str(install / "f.dat")

    def test_falls_back_to_download_dir_when_missing(self, tmp):
        install = tmp / "install"
        download = tmp / "download"
        install.mkdir()
        download.mkdir()
        # Absent from install_dir -> returns the download_dir path (regardless of existence there)
        assert dl.resolve_data_path(str(install), str(download), "missing.dat") == str(download / "missing.dat")


# ======================================================================
# atomic_output_path
# ======================================================================

class TestAtomicOutputPath:
    """Concurrent writers must not share a temp filename.

    Regression guard for a real failure: several processes calling
    siren.load_detector() against a cold cache each wrote the composite GDML
    through a fixed '<dest>.tmp' name. The first to rename pulled the file out
    from under the rest, and most of them died with FileNotFoundError.
    """

    def test_writes_and_renames(self, tmp):
        dest = tmp / "out.dat"
        with dl.atomic_output_path(str(dest)) as t:
            Path(t).write_text("payload")
            assert not dest.exists(), "must not publish before the rename"
        assert dest.read_text() == "payload"

    def test_creates_missing_parent_directory(self, tmp):
        dest = tmp / "nested" / "deeper" / "out.dat"
        with dl.atomic_output_path(str(dest)) as t:
            Path(t).write_text("x")
        assert dest.read_text() == "x"

    def test_error_leaves_dest_untouched_and_cleans_up(self, tmp):
        dest = tmp / "out.dat"
        dest.write_text("original")
        with pytest.raises(RuntimeError):
            with dl.atomic_output_path(str(dest)) as t:
                Path(t).write_text("half written")
                raise RuntimeError("boom")
        assert dest.read_text() == "original"
        assert list(tmp.iterdir()) == [dest], "temp file must be cleaned up"

    def test_temp_name_is_unique_per_call(self, tmp):
        dest = tmp / "out.dat"
        with dl.atomic_output_path(str(dest)) as a:
            with dl.atomic_output_path(str(dest)) as b:
                assert a != b, "overlapping writers must not share a temp name"
                Path(a).write_text("a")
                Path(b).write_text("b")
        assert dest.exists()

    def test_concurrent_writers_all_succeed(self, tmp):
        """Every writer completes and the destination holds one intact payload."""
        dest = tmp / "composite.gdml"
        nwriters = 12
        # Large enough that writers are still writing while others rename.
        def payload(tag):
            return "".join(f"{tag:03d}" + "x" * 80 + "\n" for _ in range(20000))

        errors: list[BaseException] = []

        def writer(tag):
            try:
                with dl.atomic_output_path(str(dest)) as t:
                    Path(t).write_text(payload(tag))
            except BaseException as e:  # noqa: BLE001 - recorded and re-raised below
                errors.append(e)

        threads = [threading.Thread(target=writer, args=(i,))
                   for i in range(nwriters)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert not errors, f"writers failed: {errors[:3]}"
        text = dest.read_text()
        tags = {line[:3] for line in text.splitlines() if line}
        assert len(tags) == 1, "destination interleaves several writers' output"
        assert text == payload(int(tags.pop())), "destination is truncated"
        assert [p.name for p in tmp.iterdir()] == [dest.name], "temp files left behind"
