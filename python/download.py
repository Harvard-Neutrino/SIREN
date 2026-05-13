"""siren-download: fetch SIREN resource data from Zenodo."""

import os
import sys
import json
import zipfile
import argparse
import tempfile
import urllib.request
import urllib.error

ZENODO_RECORD_ID = "20129082"
ZENODO_FILES_API = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}/files"


def _install_dir(override=None):
    if override:
        return override
    try:
        from siren._util import resource_package_dir, appdata_dir
        pkg_dir = resource_package_dir()
        # Prefer the package resources dir; fall back to ~/.siren if not writable
        try:
            os.makedirs(pkg_dir, exist_ok=True)
            test = os.path.join(pkg_dir, ".write_test")
            open(test, "w").close()
            os.remove(test)
            return pkg_dir
        except OSError:
            print(
                f"  Package resources dir {pkg_dir} is not writable; "
                "installing to ~/.siren instead.",
                file=sys.stderr,
            )
            return appdata_dir("siren")
    except ImportError:
        base = os.environ.get("LEPTONINJECTOR_USERDIR", os.path.expanduser("~"))
        path = os.path.join(base, ".siren")
        os.makedirs(path, exist_ok=True)
        return path


def _query_files(token=None):
    """Return list of file dicts from Zenodo API, or None on error."""
    url = ZENODO_FILES_API
    if token:
        url += f"?token={token}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
        # Zenodo API v2 returns {"entries": [...]} or a bare list
        if isinstance(data, list):
            return data
        return data.get("entries", data.get("files", []))
    except urllib.error.HTTPError as e:
        print(f"Zenodo API error: {e.code} {e.reason}", file=sys.stderr)
        if e.code == 404 and not token:
            print(
                "  Record may be restricted. Retry with --token <your-token>.",
                file=sys.stderr,
            )
        return None
    except OSError as e:
        print(f"Network error: {e}", file=sys.stderr)
        return None


def _download(url, dest, token=None, label=""):
    """Download url to dest, printing a progress bar."""
    if token:
        sep = "&" if "?" in url else "?"
        url += f"{sep}token={token}"

    done = [0]

    def _hook(count, block, total):
        done[0] = count * block
        if total > 0:
            pct = min(100, done[0] * 100 // total)
            filled = pct // 2
            bar = "#" * filled + "-" * (50 - filled)
            mb_done = done[0] / 1e6
            mb_total = total / 1e6
            print(
                f"\r  [{bar}] {pct:3d}%  {mb_done:.1f}/{mb_total:.1f} MB",
                end="",
                flush=True,
            )

    if label:
        print(f"Downloading {label} ...")
    urllib.request.urlretrieve(url, dest, _hook)
    print()  # end progress line


def _extract(zip_path, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(dest_dir)


def _file_download_url(entry, token=None):
    """Extract the content download URL from a Zenodo file entry dict."""
    links = entry.get("links", {})
    # InvenioRDM (new Zenodo) uses links.content; legacy uses links.self
    url = links.get("content") or links.get("self")
    if not url:
        raise ValueError(f"Cannot find download URL in entry: {entry}")
    return url


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_list(args):
    entries = _query_files(args.token)
    if entries is None:
        return 1
    if not entries:
        print("No files found in this Zenodo record.")
        return 0
    print(f"Files in Zenodo record {ZENODO_RECORD_ID}:")
    for e in sorted(entries, key=lambda x: x.get("key", "")):
        size_mb = e.get("size", 0) / 1e6
        print(f"  {e['key']:<50s}  {size_mb:8.1f} MB")
    return 0


def _download_matching(args, resource_type, name_filter=None):
    """
    Download Zenodo files whose key starts with resource_type (e.g. 'fluxes'
    or 'processes'), optionally filtered by name_filter substring.
    Extracts each zip into <dest>/<resource_type>/.
    """
    install_dir = _install_dir(args.dest)
    resource_dir = os.path.join(install_dir, resource_type)

    entries = _query_files(args.token)
    if entries is None:
        return 1

    # Match entries: key must contain the resource_type prefix
    matches = [
        e for e in entries
        if e.get("key", "").startswith(resource_type)
        and e.get("key", "").endswith(".zip")
        and (name_filter is None or name_filter.lower() in e["key"].lower())
    ]

    if not matches:
        desc = resource_type
        if name_filter:
            desc += f"/{name_filter}"
        print(f"No Zenodo files matched '{desc}'. Run `siren-download --list` to see what's available.")
        return 1

    for entry in matches:
        key = entry["key"]
        size_mb = entry.get("size", 0) / 1e6
        url = _file_download_url(entry, args.token)

        suffix = os.path.splitext(key)[1] or ".zip"
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp_path = tmp.name
        try:
            _download(url, tmp_path, args.token, label=f"{key} ({size_mb:.1f} MB)")
            print(f"  Extracting into {resource_dir}/ ...")
            _extract(tmp_path, resource_dir)
            print(f"  Done.")
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="siren-download",
        description="Download SIREN resource data files from Zenodo.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  siren-download --list
  siren-download --processes
  siren-download --flux
  siren-download --flux atmospheric
  siren-download --all
  siren-download --dest /path/to/resources --all
  siren-download --token <zenodo-token> --all
""",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available files on Zenodo without downloading.",
    )
    parser.add_argument(
        "--processes", action="store_true",
        help="Download all cross-section/process tables (~10 MB).",
    )
    parser.add_argument(
        "--flux", metavar="NAME", nargs="?", const="all",
        help=(
            "Download flux data. Omit NAME to download all flux files. "
            "When the Zenodo record is split by model, NAME can be a specific "
            "model (e.g. 'atmospheric', 'accelerator')."
        ),
    )
    parser.add_argument(
        "--all", dest="download_all", action="store_true",
        help="Download all resources (processes + all fluxes).",
    )
    parser.add_argument(
        "--dest", metavar="DIR",
        help="Override the install directory (default: siren package resources/ dir, "
             "or ~/.siren if the package dir is not writable).",
    )
    parser.add_argument(
        "--token", metavar="TOKEN",
        help="Zenodo access token for restricted or unpublished records.",
    )

    args = parser.parse_args()

    if not any([args.list, args.processes, args.flux, args.download_all]):
        parser.print_help()
        return 0

    if args.list:
        return cmd_list(args)

    ret = 0
    if args.download_all:
        ret |= _download_matching(args, "processes")
        ret |= _download_matching(args, "fluxes")
    else:
        if args.processes:
            ret |= _download_matching(args, "processes")
        if args.flux:
            name_filter = None if args.flux.lower() == "all" else args.flux
            ret |= _download_matching(args, "fluxes", name_filter=name_filter)

    return ret


if __name__ == "__main__":
    sys.exit(main())
