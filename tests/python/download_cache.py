"""Keep pinned test downloads in a local directory between runs.

``pytest --download-cache DIR`` wraps ``siren.download.ensure_files`` so that a
file with a ``sha256`` pin is copied from ``DIR/<sha256>`` instead of being
downloaded, and is stored there after a download. Entries are named by their
digest and checked against it in both directions, so a stale or corrupt entry
is never used and a cache shared between branches cannot serve the wrong file.
Unpinned files bypass the cache.
"""
from __future__ import annotations

import hashlib
import shutil
from pathlib import Path


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _verified(path: Path, digest: str) -> bool:
    return path.is_file() and sha256_of(path) == digest


def cached_ensure_files(ensure_files, cache: Path):
    """Return *ensure_files* reading pinned files from, and storing them in, *cache*."""
    cache = Path(cache)

    def wrapper(specs):
        pinned = [(Path(s["path"]), s["sha256"]) for s in specs if s.get("sha256")]
        for dest, digest in pinned:
            if not dest.exists() and _verified(cache / digest, digest):
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(cache / digest, dest)
        ensure_files(specs)
        cache.mkdir(parents=True, exist_ok=True)
        for dest, digest in pinned:
            if _verified(dest, digest) and not _verified(cache / digest, digest):
                shutil.copyfile(dest, cache / digest)

    return wrapper
