#!/bin/bash
# Download the GSL source tarball with mirror fallback.
#  Used both by the GitHub Actions populate_cache job and as a
# defensive backstop inside cibw_before_all.sh when the cache misses.
#
# Usage: fetch_gsl.sh <version> <output-path>
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <version> <output-path>" >&2
    exit 2
fi

GSL_VERSION="$1"
OUT="$2"

GSL_TARBALL="gsl-${GSL_VERSION}.tar.gz"
GSL_MIRRORS=(
    "https://ftpmirror.gnu.org/gsl/${GSL_TARBALL}"
    "https://ftp.gnu.org/gnu/gsl/${GSL_TARBALL}"
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
    echo "All GSL mirrors exhausted" >&2
    return 1
}

mkdir -p "$(dirname "$OUT")"
fetch_with_fallback "$OUT" "${GSL_MIRRORS[@]}"
