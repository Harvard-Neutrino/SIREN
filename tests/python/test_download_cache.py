"""The --download-cache wrapper serves and stores pinned files by digest."""
from __future__ import annotations

import hashlib

from download_cache import cached_ensure_files

PAYLOAD = b"<gdml/>\n"
DIGEST = hashlib.sha256(PAYLOAD).hexdigest()


def _downloader(calls):
    """Stand-in for ensure_files: 'downloads' each missing file and records it."""
    def ensure_files(specs):
        for spec in specs:
            dest = spec["path"]
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(PAYLOAD)
                calls.append(dest)
    return ensure_files


def _spec(dest, digest=DIGEST):
    return [{"path": dest, "url": "https://example.invalid/a.gdml", "sha256": digest}]


def test_cold_cache_downloads_and_stores(tmp_path):
    calls, cache, dest = [], tmp_path / "cache", tmp_path / "run" / "gdml" / "a.gdml"
    cached_ensure_files(_downloader(calls), cache)(_spec(dest))
    assert calls == [dest]
    assert (cache / DIGEST).read_bytes() == PAYLOAD


def test_warm_cache_serves_without_downloading(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / DIGEST).write_bytes(PAYLOAD)
    calls, dest = [], tmp_path / "run" / "gdml" / "a.gdml"
    cached_ensure_files(_downloader(calls), cache)(_spec(dest))
    assert calls == []
    assert dest.read_bytes() == PAYLOAD


def test_corrupt_entry_is_ignored_and_replaced(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / DIGEST).write_bytes(b"truncated")
    calls, dest = [], tmp_path / "run" / "a.gdml"
    cached_ensure_files(_downloader(calls), cache)(_spec(dest))
    assert calls == [dest]
    assert dest.read_bytes() == PAYLOAD
    assert (cache / DIGEST).read_bytes() == PAYLOAD


def test_a_wrong_local_file_does_not_enter_the_cache(tmp_path):
    cache, dest = tmp_path / "cache", tmp_path / "run" / "a.gdml"
    dest.parent.mkdir(parents=True)
    dest.write_bytes(b"edited by hand")
    cached_ensure_files(_downloader([]), cache)(_spec(dest))
    assert not (cache / DIGEST).exists()


def test_unpinned_files_bypass_the_cache(tmp_path):
    calls, cache, dest = [], tmp_path / "cache", tmp_path / "run" / "a.gdml"
    cached_ensure_files(_downloader(calls), cache)(_spec(dest, digest=""))
    assert calls == [dest]
    assert list(cache.iterdir()) == []
