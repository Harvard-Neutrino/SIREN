"""
siren.download -- unified file download and integrity verification.

Core primitives for downloading data files needed by SIREN resource
loaders (detectors, fluxes, cross-section splines). Each resource
loader declares its own file specs; this module provides the shared
download/verify/extract machinery.

Public API
----------
writable_data_dir(module_dir)
    Return module_dir if writable, else a fallback under ~/.siren/.

download_file(url, dest, sha256="")
    Download a single file with atomic write and optional SHA-256 check.

ensure_files(specs)
    Download any missing files from a list of {path, url, sha256} specs.

ensure_zenodo_files(record_id, files, dest_dir, token="")
    Download individual files from a Zenodo record into dest_dir.

ensure_zenodo_archive(record_id, filename, dest_dir, token="", prefix="")
    Download and extract a zip archive from a Zenodo record.

zenodo_file_url(record_id, filename, token="")
    Resolve the download URL for one file in a Zenodo record.

CLI
---
siren-download --list              List installed resources
siren-download --fetch-all         Pre-fetch data for all installed resources
siren-download --fetch NAME        Pre-fetch data for one resource
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
import zipfile
from typing import Any


def writable_data_dir(module_dir: str) -> str:
    """Return a writable directory for storing downloaded data.

    If *module_dir* (typically ``os.path.dirname(__file__)`` in a loader
    script) is writable, returns it directly. Otherwise falls back to a
    mirror path under ``~/.siren/resources/``, preserving the relative
    directory structure from the ``resources/`` root downward.

    This handles the case where SIREN is installed into a read-only
    location (e.g. system site-packages).
    """
    try:
        os.makedirs(module_dir, exist_ok=True)
        probe = os.path.join(module_dir, ".write_probe")
        with open(probe, "w") as f:
            f.write("")
        os.remove(probe)
        return module_dir
    except OSError:
        pass

    # Fall back to ~/.siren, mirroring the path from resources/ downward
    norm = os.path.normpath(module_dir)
    parts = norm.split(os.sep)
    try:
        idx = parts.index("resources")
        rel = os.path.join(*parts[idx:])
    except ValueError:
        rel = os.path.basename(norm)

    fallback = os.path.join(os.path.expanduser("~"), ".siren", rel)
    os.makedirs(fallback, exist_ok=True)
    return fallback


def resolve_data_path(install_dir: str, download_dir: str,
                      filename: str) -> str:
    """Resolve a data file path, preferring *install_dir* when the file
    already exists there.

    Resource loaders download files into *download_dir* (the writable
    location returned by :func:`writable_data_dir`).  When the package
    is installed into a read-only tree the writable directory differs
    from the original install directory.  If the file was shipped with
    the package (or cached from a previous writable install), it may
    still be present at *install_dir* -- use it directly and avoid a
    redundant download.
    """
    local = os.path.join(install_dir, filename)
    if os.path.isfile(local):
        return local
    return os.path.join(download_dir, filename)


# ======================================================================
# Core download primitives
# ======================================================================

DEFAULT_TIMEOUT = 30  # seconds


def download_file(url: str, dest: str, sha256: str = "",
                  label: str = "", show_progress: bool = True,
                  timeout: float = DEFAULT_TIMEOUT) -> None:
    """Download a file from *url* to *dest*.

    - Atomic write: writes to a .tmp file, renames on success.
    - Follows HTTP redirects (GitHub release assets, Zenodo CDN).
    - Optional SHA-256 verification.
    - Progress bar on stdout.
    """
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    tmp = dest + ".tmp"

    if not label:
        label = os.path.basename(dest)

    try:
        if show_progress:
            print(f"  Downloading {label} ...")

        req = urllib.request.Request(url, headers={"User-Agent": "siren"})
        h = hashlib.sha256()
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0

            with open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    out.write(chunk)
                    h.update(chunk)
                    downloaded += len(chunk)

                    if show_progress and total > 0:
                        pct = min(100, downloaded * 100 // total)
                        filled = pct // 2
                        bar = "#" * filled + "-" * (50 - filled)
                        mb_done = downloaded / 1e6
                        mb_total = total / 1e6
                        print(
                            f"\r  [{bar}] {pct:3d}%  "
                            f"{mb_done:.1f}/{mb_total:.1f} MB",
                            end="", flush=True)

            if show_progress and total > 0:
                print()

        if sha256:
            got = h.hexdigest()
            if got != sha256:
                raise RuntimeError(
                    f"SHA-256 mismatch for {label}:\n"
                    f"  expected: {sha256}\n"
                    f"  got:      {got}")

        os.replace(tmp, dest)

    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def ensure_files(specs: list[dict[str, Any]]) -> None:
    """Download any missing files described by *specs*.

    Each spec is a dict with keys:
        path   -- absolute local file path (required)
        url    -- download URL (optional; error if missing and file absent)
        sha256 -- expected SHA-256 hex digest (optional)

    Files that already exist on disk are skipped.
    """
    missing_no_url: list[str] = []

    for spec in specs:
        path = spec["path"]
        if os.path.isfile(path):
            continue
        url = spec.get("url")
        if not url:
            missing_no_url.append(path)
            continue
        download_file(url, path, sha256=spec.get("sha256", ""))

    if missing_no_url:
        raise FileNotFoundError(
            "Missing data files with no download URL configured:\n"
            + "\n".join(f"  - {p}" for p in missing_no_url)
            + "\n\nPlace or symlink the files manually.")


# ======================================================================
# Zenodo helpers
# ======================================================================

def _zenodo_api(record_id: str, token: str = "") -> list[dict]:
    """Query the Zenodo files API for a record. Returns file entry dicts."""
    url = f"https://zenodo.org/api/records/{record_id}/files"
    if token:
        url += f"?token={token}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
        if isinstance(data, list):
            return data
        return data.get("entries", data.get("files", []))
    except urllib.error.HTTPError as e:
        if e.code == 404 and not token:
            raise FileNotFoundError(
                f"Zenodo record {record_id} not found (may be restricted; "
                "pass token= for access)") from e
        raise
    except OSError as e:
        raise ConnectionError(f"Network error querying Zenodo: {e}") from e


def zenodo_file_url(record_id: str, filename: str,
                    token: str = "") -> str:
    """Resolve the download URL for one file in a Zenodo record."""
    entries = _zenodo_api(record_id, token)
    for entry in entries:
        if entry.get("key") == filename:
            links = entry.get("links", {})
            url = links.get("content") or links.get("self")
            if url:
                return url
            raise ValueError(
                f"Zenodo entry for '{filename}' has no download link")
    available = sorted(e.get("key", "?") for e in entries)
    raise FileNotFoundError(
        f"'{filename}' not found in Zenodo record {record_id}.\n"
        f"Available: {', '.join(available)}")


def ensure_zenodo_files(record_id: str,
                        files: list[dict[str, str]],
                        dest_dir: str,
                        token: str = "") -> None:
    """Download individual files from a Zenodo record.

    Each entry in *files* is a dict with keys:
        name   -- filename in the Zenodo record
        path   -- relative path under dest_dir to save to (optional; defaults to name)
        sha256 -- expected hash (optional)

    Only downloads files not already present at dest_dir/path.
    """
    entries_cache: list[dict] | None = None

    for fspec in files:
        name = fspec["name"]
        rel_path = fspec.get("path", name)
        dest = os.path.join(dest_dir, rel_path)
        if os.path.isfile(dest):
            continue

        # Lazy-fetch the API listing
        if entries_cache is None:
            entries_cache = _zenodo_api(record_id, token)

        url = None
        for entry in entries_cache:
            if entry.get("key") == name:
                links = entry.get("links", {})
                url = links.get("content") or links.get("self")
                break
        if not url:
            raise FileNotFoundError(
                f"'{name}' not found in Zenodo record {record_id}")

        download_file(url, dest, sha256=fspec.get("sha256", ""),
                      label=name)


def _sha256_file(path: str) -> str:
    """Compute URL-safe base64-encoded SHA-256 digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 16), b""):
            h.update(block)
    return base64.urlsafe_b64encode(h.digest()).rstrip(b"=").decode("ascii")


def _zip_cache_dir(dest_dir: str) -> str:
    return os.path.join(dest_dir, ".download_cache")


def _cached_zip(dest_dir: str, record_id: str, filename: str,
                token: str = "") -> str:
    """Return the path to a cached zip, downloading if necessary.

    The cache filename is the full SHA-256 of the zip content (URL-safe
    base64, 43 chars) so that different versions never collide.
    """
    cache_dir = _zip_cache_dir(dest_dir)
    os.makedirs(cache_dir, exist_ok=True)

    # Check for an existing cached file for this record/filename
    tag = f"{record_id}_{filename}"
    existing = sorted(
        (f for f in os.listdir(cache_dir)
         if f.startswith(tag + ".") and f.endswith(".zip")),
        key=lambda f: os.path.getmtime(os.path.join(cache_dir, f)),
        reverse=True,
    )
    if existing:
        return os.path.join(cache_dir, existing[0])

    # Download, hash, rename into cache
    url = zenodo_file_url(record_id, filename, token)
    tmp_path = os.path.join(cache_dir, f"{tag}.{os.getpid()}.tmp")
    try:
        download_file(url, tmp_path, label=filename)
        sha = _sha256_file(tmp_path)
        cache_path = os.path.join(cache_dir, f"{tag}.{sha}.zip")
        os.replace(tmp_path, cache_path)
        return cache_path
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


_EXTRACTED_SENTINEL = ".zenodo_extracted"


def _safe_extract(zf: zipfile.ZipFile, dest_dir: str,
                  members: list[str] | None = None) -> None:
    """Extract zip members, rejecting paths that escape dest_dir (Zip Slip).

    Extracts one member at a time, re-validating each member's resolved
    path before extracting it (realpath reflects symlinks created by
    earlier members) so that earlier entries cannot redirect later ones
    outside dest_dir.
    """
    abs_dest = os.path.realpath(dest_dir)
    names = members if members is not None else zf.namelist()
    for name in names:
        target = os.path.realpath(os.path.join(dest_dir, name))
        if not target.startswith(abs_dest + os.sep) and target != abs_dest:
            raise ValueError(
                f"Zip entry '{name}' would extract outside {dest_dir}")
        zf.extract(name, dest_dir)


def ensure_zenodo_archive(record_id: str, filename: str, dest_dir: str,
                          token: str = "", prefix: str = "") -> None:
    """Download a zip archive from a Zenodo record and extract it.

    Downloads *filename* from Zenodo record *record_id*, extracts into
    *dest_dir*. The zip is cached locally (keyed by SHA-256) so that
    multiple loaders sharing the same archive only download once.

    If *prefix* is given (e.g. ``"processes/DipoleHNLDISSplines"``),
    only zip entries under that prefix are extracted. A sentinel file
    is written after successful extraction so subsequent calls skip
    the download even when the directory already contains other files
    (like loader scripts).

    If *prefix* is empty, the entire zip is extracted and extraction
    is skipped when the sentinel exists.
    """
    if prefix:
        sentinel = os.path.join(dest_dir, prefix, _EXTRACTED_SENTINEL)
    else:
        sentinel = os.path.join(dest_dir, _EXTRACTED_SENTINEL)
    if os.path.isfile(sentinel):
        return

    zip_path = _cached_zip(dest_dir, record_id, filename, token)

    os.makedirs(dest_dir, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        if prefix:
            pfx = prefix.rstrip("/") + "/"
            members = [m for m in zf.namelist()
                       if m == prefix or m.startswith(pfx)]
            if not members:
                available = sorted(set(
                    n.split("/")[0] for n in zf.namelist() if "/" in n))
                raise FileNotFoundError(
                    f"Prefix '{prefix}' not found in {filename}.\n"
                    f"Top-level directories: {', '.join(available)}")
            print(f"  Extracting {prefix} from {filename} "
                  f"({len(members)} entries) ...")
            _safe_extract(zf, dest_dir, members)
        else:
            print(f"  Extracting {filename} ...")
            _safe_extract(zf, dest_dir)

    os.makedirs(os.path.dirname(sentinel), exist_ok=True)
    with open(sentinel, "w") as f:
        f.write(f"zenodo:{record_id}/{filename}\n")


# ======================================================================
# CLI
# ======================================================================

def _cmd_list() -> int:
    """List installed resources that have loader scripts."""
    try:
        from . import _util
    except ImportError:
        print("Cannot import siren._util; is siren installed?",
              file=sys.stderr)
        return 1

    for category, lister in [("detectors", _util.list_detectors),
                             ("fluxes", _util.list_fluxes),
                             ("processes", _util.list_processes)]:
        names = lister()
        if names:
            print(f"{category}:")
            for n in names:
                print(f"  {n}")
            print()
    return 0


def _cmd_fetch(names: list[str]) -> int:
    """Pre-fetch data for the named resources.

    Imports each resource module and calls its ``fetch_data()`` function
    if defined. This triggers ``ensure_zenodo_archive`` /
    ``ensure_files`` without requiring domain-specific load arguments
    (like flux tags or mass parameters).
    """
    try:
        from . import _util
    except ImportError:
        print("Cannot import siren._util; is siren installed?",
              file=sys.stderr)
        return 1

    # Build a map of name -> resource_type for all installed resources.
    # A name may appear in multiple categories; store all of them.
    resource_map: dict[str, list[str]] = {}
    for resource_type, lister in [("detector", _util.list_detectors),
                                  ("flux", _util.list_fluxes),
                                  ("processes", _util.list_processes)]:
        for n in lister():
            resource_map.setdefault(n, []).append(resource_type)

    if not names:
        names = list(resource_map.keys())

    ret = 0
    for name in names:
        print(f"--- {name} ---")
        resource_types = resource_map.get(name)
        if resource_types is None:
            print(f"  (unknown resource '{name}')")
            ret = 1
            continue
        for resource_type in resource_types:
            try:
                mod = _util.import_resource(resource_type, name)
                if mod is None:
                    print(f"  (no {resource_type} loader script)")
                    continue
                if hasattr(mod, "fetch_data"):
                    mod.fetch_data()
                    print(f"  OK ({resource_type})")
                else:
                    print(f"  (no fetch_data function for {resource_type}, skipping)")
            except Exception as e:
                print(f"  Error ({resource_type}): {e}", file=sys.stderr)
                ret = 1
    return ret


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="siren-download",
        description="Manage SIREN resource data files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  siren-download --list                 List installed resources
  siren-download --fetch-all            Pre-fetch data for all resources
  siren-download --fetch SBN            Pre-fetch data for one resource
""")
    parser.add_argument("--list", action="store_true",
                        help="List installed resource loaders.")
    parser.add_argument("--fetch-all", action="store_true",
                        help="Pre-fetch data for all installed resources.")
    parser.add_argument("--fetch", metavar="NAME", nargs="+",
                        help="Pre-fetch data for specific resource(s).")

    args = parser.parse_args()

    if not any([args.list, args.fetch_all, args.fetch]):
        parser.print_help()
        return 0

    try:
        if args.list:
            return _cmd_list()

        names = args.fetch or []
        if args.fetch_all:
            names = []  # empty = all
        return _cmd_fetch(names)
    except (OSError, ConnectionError, urllib.error.URLError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
