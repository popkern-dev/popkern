#!/usr/bin/env bash
# Line count summary for popkern workspace.
# Usage: tools/loc.sh [crate-name]
#   tools/loc.sh           — full workspace summary
#   tools/loc.sh bool      — popkern-bool only
#   tools/loc.sh core      — popkern-core only

set -euo pipefail

if [ $# -eq 0 ]; then
    echo "=== Full Workspace ==="
    cloc crates/ --quiet
    echo ""
    echo "=== Per Crate ==="
    for d in crates/popkern-*/src; do
        crate=$(basename "$(dirname "$d")")
        code=$(cloc "$d" --quiet --csv 2>/dev/null | tail -1 | cut -d, -f5)
        printf "  %-20s %s lines\n" "$crate" "${code:-0}"
    done
else
    crate="popkern-$1"
    dir="crates/$crate/src"
    if [ ! -d "$dir" ]; then
        echo "Unknown crate: $crate (no $dir)"
        exit 1
    fi
    echo "=== $crate ==="
    cloc "$dir"
fi
