#!/usr/bin/env bash
# Blast radius estimator for popkern.
# Usage:
#   tools/blast-radius.sh                    # analyze uncommitted changes (git diff)
#   tools/blast-radius.sh --staged           # analyze staged changes
#   tools/blast-radius.sh src/foo.rs ...     # analyze specific files
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

# ---------------------------------------------------------------------------
# 1. Determine changed files
# ---------------------------------------------------------------------------

STAGED=false
FILES=()

if [[ $# -eq 0 ]]; then
    mapfile -t FILES < <(git diff --name-only 2>/dev/null || true)
elif [[ "${1:-}" == "--staged" ]]; then
    STAGED=true
    mapfile -t FILES < <(git diff --cached --name-only 2>/dev/null || true)
else
    FILES=("$@")
fi

# Filter to only .rs files under crates/
RS_FILES=()
for f in "${FILES[@]}"; do
    # Normalise: strip leading ./ if present
    f="${f#./}"
    if [[ "$f" == crates/*.rs ]]; then
        RS_FILES+=("$f")
    fi
done

if [[ ${#RS_FILES[@]} -eq 0 ]]; then
    echo "=== Blast Radius Analysis ==="
    echo ""
    echo "No changed .rs files detected under crates/."
    if [[ $STAGED == false && $# -eq 0 ]]; then
        echo "Try: tools/blast-radius.sh --staged   (for staged changes)"
        echo "     tools/blast-radius.sh path/to/file.rs"
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# 2. Identify directly affected crates
# ---------------------------------------------------------------------------

declare -A DIRECT_CRATES   # crate_name -> 1

for f in "${RS_FILES[@]}"; do
    # crates/<crate-name>/...
    crate=$(echo "$f" | sed 's|^crates/\([^/]*\)/.*|\1|')
    DIRECT_CRATES["$crate"]=1
done

# ---------------------------------------------------------------------------
# 3. Build dependency map: crate -> space-separated list of its dependencies
#    Then find downstream dependents (crates that depend on affected crates).
# ---------------------------------------------------------------------------

# dep_of[A] = "B C"  means A depends on B and C
declare -A DEP_OF

for cargo_toml in "$ROOT"/crates/*/Cargo.toml; do
    crate=$(grep -m1 '^name' "$cargo_toml" | sed 's/name *= *"\(.*\)"/\1/')
    # Collect all popkern-* dependencies listed in the file
    deps=$(grep -oE 'popkern-[a-z]+' "$cargo_toml" | grep -v "^${crate}$" | sort -u | tr '\n' ' ' || true)
    DEP_OF["$crate"]="$deps"
done

# Build transitive closure of "affected" starting from DIRECT_CRATES.
# A crate is affected if it IS a directly changed crate, or it depends (directly
# or transitively) on an affected crate.

declare -A AFFECTED         # all affected crates (direct + downstream)
for c in "${!DIRECT_CRATES[@]}"; do
    AFFECTED["$c"]=1
done

# Iterate to fixpoint
changed=true
while $changed; do
    changed=false
    for crate in "${!DEP_OF[@]}"; do
        if [[ -n "${AFFECTED[$crate]+set}" ]]; then
            continue   # already marked
        fi
        deps="${DEP_OF[$crate]}"
        for dep in $deps; do
            if [[ -n "${AFFECTED[$dep]+set}" ]]; then
                AFFECTED["$crate"]=1
                changed=true
                break
            fi
        done
    done
done

# Downstream = AFFECTED minus DIRECT
declare -A DOWNSTREAM
for c in "${!AFFECTED[@]}"; do
    if [[ -z "${DIRECT_CRATES[$c]+set}" ]]; then
        DOWNSTREAM["$c"]=1
    fi
done

# ---------------------------------------------------------------------------
# 4. Count files and tests per crate
# ---------------------------------------------------------------------------

crate_file_count() {
    local crate="$1"
    local count
    count=$(find "$ROOT/crates/$crate/src" -name "*.rs" 2>/dev/null | wc -l)
    echo "${count// /}"
}

crate_test_count() {
    local crate="$1"
    local count
    count=$(grep -r '#\[test\]' "$ROOT/crates/$crate/src" 2>/dev/null | wc -l) || count=0
    echo "${count// /}"
}

# ---------------------------------------------------------------------------
# 5. Detect public API changes in the diff
# ---------------------------------------------------------------------------

if [[ $STAGED == true ]]; then
    DIFF_TEXT=$(git diff --cached 2>/dev/null || true)
else
    DIFF_TEXT=$(git diff 2>/dev/null || true)
fi

# If files were passed explicitly, diff them against HEAD
if [[ $# -gt 0 && "${1:-}" != "--staged" ]]; then
    DIFF_TEXT=$(git diff HEAD -- "${RS_FILES[@]}" 2>/dev/null || true)
    # For untracked/new files read them directly
    for f in "${RS_FILES[@]}"; do
        if [[ -f "$ROOT/$f" ]]; then
            DIFF_TEXT+=$(cat "$ROOT/$f" 2>/dev/null || true)
        fi
    done
fi

PUB_API_CHANGED=false
PUB_API_REASON=""
if echo "$DIFF_TEXT" | grep -qE '^\+.*pub[[:space:]]+(trait|enum[[:space:]]+\w+Error)' 2>/dev/null; then
    PUB_API_CHANGED=true
    PUB_API_REASON="pub trait or pub enum Error in diff"
fi

# ---------------------------------------------------------------------------
# 6. Rigor level
# ---------------------------------------------------------------------------

num_direct=${#DIRECT_CRATES[@]}
num_downstream=${#DOWNSTREAM[@]}
num_total_crates=$(( num_direct + num_downstream ))
num_changed_files=${#RS_FILES[@]}

RIGOR="Light"
RIGOR_REASON=""

if [[ -n "${AFFECTED[popkern-core]+set}" ]]; then
    RIGOR="Full"
    RIGOR_REASON="popkern-core affected"
elif $PUB_API_CHANGED; then
    RIGOR="Full"
    RIGOR_REASON="public API changes detected"
elif [[ $num_total_crates -ge 2 ]]; then
    RIGOR="Heavy"
    RIGOR_REASON="2+ crates affected"
elif [[ $num_changed_files -ge 4 ]]; then
    RIGOR="Medium"
    RIGOR_REASON="4+ files changed in 1 crate"
else
    RIGOR="Light"
    RIGOR_REASON="1 crate, 1-3 files"
fi

# Upgrade to Full if both Heavy conditions AND pub API changed
if [[ "$RIGOR" == "Heavy" ]] && $PUB_API_CHANGED; then
    RIGOR="Full"
    RIGOR_REASON="2+ crates affected + public API changes"
fi

# ---------------------------------------------------------------------------
# 7. Output
# ---------------------------------------------------------------------------

echo "=== Blast Radius Analysis ==="
echo ""
echo "Changed files: ${#RS_FILES[@]}"
for f in "${RS_FILES[@]}"; do
    echo "  $f"
done

echo ""
echo "Directly affected crates: $num_direct"
for crate in $(echo "${!DIRECT_CRATES[@]}" | tr ' ' '\n' | sort); do
    fc=$(crate_file_count "$crate")
    tc=$(crate_test_count "$crate")
    echo "  $crate ($fc files, $tc tests)"
done

if [[ $num_downstream -gt 0 ]]; then
    echo ""
    echo "Downstream dependents: $num_downstream"
    for crate in $(echo "${!DOWNSTREAM[@]}" | tr ' ' '\n' | sort); do
        deps="${DEP_OF[$crate]:-}"
        # Filter dep list to only directly affected crates (for readability)
        affected_deps=()
        for dep in $deps; do
            if [[ -n "${DIRECT_CRATES[$dep]+set}" ]]; then
                affected_deps+=("${dep#popkern-}")
            fi
        done
        if [[ ${#affected_deps[@]} -gt 0 ]]; then
            dep_str=$(printf '%s, ' "${affected_deps[@]}")
            dep_str="${dep_str%, }"
            echo "  $crate (depends on: $dep_str)"
        else
            echo "  $crate"
        fi
    done
else
    echo ""
    echo "Downstream dependents: 0"
fi

# Total test count across all affected crates
total_tests=0
for crate in "${!AFFECTED[@]}"; do
    tc=$(crate_test_count "$crate")
    total_tests=$(( total_tests + tc ))
done

echo ""
echo "Total blast radius: $num_total_crates crates, ${total_tests}+ tests"
echo ""

if $PUB_API_CHANGED; then
    echo "Public API changes detected: YES ($PUB_API_REASON)"
else
    echo "Public API changes detected: no"
fi

echo ""
echo "Recommended rigor: $RIGOR"
echo "  Reason: $RIGOR_REASON"
