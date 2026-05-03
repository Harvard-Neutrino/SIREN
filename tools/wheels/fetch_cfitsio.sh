#!/bin/bash
# Download the cfitsio source tarball with mirror fallback and integrity
# validation. Used both by the GitHub Actions populate_cache job and as a
# defensive backstop inside cibw_before_all.sh when the cache misses.
#
# Usage: fetch_cfitsio.sh <version> <output-path>
#
# The heasarc mirror is the official source. The github mirror serves the
# same source tree (different wrapper directory and a few extra
# .github/.gitignore files, normalized at extract time via tar
# --strip-components=1 by the caller).
#
# Implemented as a portable bash retry loop because curl's --retry-all-errors
# was added in 7.71 and the manylinux2014 base image still ships curl 7.29.
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <version> <output-path>" >&2
    exit 2
fi

CFITSIO_VERSION="$1"
OUT="$2"

CFITSIO_TARBALL="cfitsio-${CFITSIO_VERSION}.tar.gz"
CFITSIO_MIRRORS=(
    "https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/${CFITSIO_TARBALL}"
    "https://github.com/HEASARC/cfitsio/archive/refs/tags/${CFITSIO_TARBALL}"
)

fetch_with_fallback() {
    local out="$1"
    shift
    local max_attempts=3
    local delay=10
    local url attempt rc
    for url in "$@"; do
        echo "Trying $url" >&2
        for attempt in $(seq 1 $max_attempts); do
            rc=0
            curl --fail --location --connect-timeout 30 "$url" --output "$out" || rc=$?
            if [ "$rc" -eq 0 ] && tar -tzf "$out" >/dev/null 2>&1; then
                return 0
            fi
            if [ "$rc" -eq 0 ]; then
                echo "  downloaded file is not a valid gzip tarball" >&2
            else
                echo "  curl exited $rc on attempt $attempt/$max_attempts" >&2
            fi
            rm -f "$out"
            if [ "$attempt" -lt "$max_attempts" ]; then
                sleep "$delay"
            fi
        done
    done
    echo "All cfitsio mirrors exhausted" >&2
    return 1
}

mkdir -p "$(dirname "$OUT")"
fetch_with_fallback "$OUT" "${CFITSIO_MIRRORS[@]}"
