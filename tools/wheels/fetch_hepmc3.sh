#!/bin/bash
# Download the HepMC3 source tarball with mirror fallback.
#  Used both by the GitHub Actions populate_cache job and as a
# defensive backstop inside cibw_before_all.sh when the cache misses.
#
# Usage: fetch_hepmc3.sh <version> <output-path>
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <version> <output-path>" >&2
    exit 2
fi

HEPMC3_VERSION="$1"
OUT="$2"

HEPMC3_MIRRORS=(
    "https://hepmc.web.cern.ch/hepmc/releases/HepMC3-${HEPMC3_VERSION}.tar.gz"
    "https://gitlab.cern.ch/hepmc/HepMC3/-/archive/${HEPMC3_VERSION}/HepMC3-${HEPMC3_VERSION}.tar.gz"
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
    echo "All HepMC3 mirrors exhausted" >&2
    return 1
}

mkdir -p "$(dirname "$OUT")"
fetch_with_fallback "$OUT" "${HEPMC3_MIRRORS[@]}"
