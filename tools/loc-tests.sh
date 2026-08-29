#!/usr/bin/env bash
# Line count split: test code vs non-test code.
#
# Bucketing follows popkern's conventions:
#   * Non-test code   = files under crates/*/src/ NOT matching *_tests.rs
#   * Test (src)      = crates/*/src/*_tests.rs (the #[path] attribute pattern
#                       used to keep source files under the 500-line cap;
#                       see CLAUDE.md "Known Gotchas")
#   * Integration     = everything under crates/*/tests/
#
# Inline `#[cfg(test)] mod tests { … }` blocks count as non-test under this
# scheme — fine in practice because the codebase already extracts test
# modules to *_tests.rs once they grow.
#
# Usage:
#   tools/loc-tests.sh           — workspace summary table
#   tools/loc-tests.sh fillet    — popkern-blend detail (per-bucket cloc)

set -euo pipefail

# Run cloc and return the SUM row's "code" column (5th field).  Returns 0
# when the path is empty / no files match the filter.
code_lines() {
    local out code
    out=$(cloc --quiet --csv "$@" 2>/dev/null || true)
    code=$(printf '%s\n' "$out" | tail -1 | cut -d, -f5)
    [[ "$code" =~ ^[0-9]+$ ]] && echo "$code" || echo 0
}

count_crate() {
    local crate_dir=$1
    local non_test=0 test_src=0 integ=0
    if [ -d "$crate_dir/src" ]; then
        non_test=$(code_lines "$crate_dir/src" --not-match-f='_tests\.rs$')
        test_src=$(code_lines "$crate_dir/src" --match-f='_tests\.rs$')
    fi
    if [ -d "$crate_dir/tests" ]; then
        integ=$(code_lines "$crate_dir/tests")
    fi
    echo "$non_test $test_src $integ"
}

if [ $# -eq 0 ]; then
    printf "%-22s %12s %12s %14s %12s %8s\n" \
        "crate" "non-test" "test (src)" "integration" "total" "test%"
    printf "%-22s %12s %12s %14s %12s %8s\n" \
        "-----" "--------" "----------" "-----------" "-----" "-----"
    sum_nt=0; sum_ts=0; sum_in=0
    for d in crates/popkern-*; do
        [ -d "$d/src" ] || continue
        crate=$(basename "$d")
        read -r nt ts integ < <(count_crate "$d")
        total=$((nt + ts + integ))
        if [ "$total" -gt 0 ]; then
            test_pct=$(awk -v t="$ts" -v i="$integ" -v tot="$total" \
                'BEGIN { printf "%.0f%%", 100 * (t + i) / tot }')
        else
            test_pct="-"
        fi
        printf "%-22s %12d %12d %14d %12d %8s\n" \
            "$crate" "$nt" "$ts" "$integ" "$total" "$test_pct"
        sum_nt=$((sum_nt + nt))
        sum_ts=$((sum_ts + ts))
        sum_in=$((sum_in + integ))
    done
    printf "%-22s %12s %12s %14s %12s %8s\n" \
        "-----" "--------" "----------" "-----------" "-----" "-----"
    total=$((sum_nt + sum_ts + sum_in))
    if [ "$total" -gt 0 ]; then
        test_pct=$(awk -v t="$sum_ts" -v i="$sum_in" -v tot="$total" \
            'BEGIN { printf "%.0f%%", 100 * (t + i) / tot }')
    else
        test_pct="-"
    fi
    printf "%-22s %12d %12d %14d %12d %8s\n" \
        "TOTAL" "$sum_nt" "$sum_ts" "$sum_in" "$total" "$test_pct"
else
    crate="popkern-$1"
    dir="crates/$crate"
    if [ ! -d "$dir" ]; then
        echo "Unknown crate: $crate (no $dir)" >&2
        exit 1
    fi
    echo "=== $crate — non-test source (src/, excluding *_tests.rs) ==="
    cloc "$dir/src" --not-match-f='_tests\.rs$' --quiet
    echo
    echo "=== $crate — test source (src/*_tests.rs) ==="
    cloc "$dir/src" --match-f='_tests\.rs$' --quiet
    if [ -d "$dir/tests" ]; then
        echo
        echo "=== $crate — integration tests (tests/) ==="
        cloc "$dir/tests" --quiet
    fi
fi
