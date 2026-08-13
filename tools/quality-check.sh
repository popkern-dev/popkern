#!/usr/bin/env bash
# Full quality gate for popkern.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mode="${1:---full}"

case "$mode" in
    --full|--fast|--perf)
        ;;
    -h|--help)
        cat <<'USAGE'
Usage: tools/quality-check.sh [--full|--fast|--perf]

--full  Run fmt, semantic-name scan, structural tests, clippy, and all tests.
--fast  Run fmt, semantic-name scan, structural tests, and cargo check.
--perf  Run the release-mode performance budget gate (PK-190): the
        wall-clock budget suites and the §17 cancellation-latency
        test — the same gate CI's perf job enforces.
USAGE
        exit 0
        ;;
    *)
        echo "Usage: tools/quality-check.sh [--full|--fast|--perf]" >&2
        exit 2
        ;;
esac

run_step() {
    local name="$1"
    shift
    echo "== $name =="
    "$@"
    echo "PASS: $name"
    echo
}

if [[ "$mode" == "--perf" ]]; then
    run_step "fillet perf budgets + §17 cancellation latency" \
        cargo test --release -p popkern-fillet \
        --test perf_budget --test cancellation_latency --test perf_acceptance \
        -- --nocapture
    run_step "smooth-mesh perf budgets" \
        cargo test --release -p popkern-part --test perf_budget -- --nocapture
    echo "Quality gate passed ($mode)."
    exit 0
fi

run_step "format" cargo fmt --check
run_step "semantic names" tools/check-semantic-names.sh
run_step "structural tests" cargo test -p popkern-structural-tests

if [[ "$mode" == "--fast" ]]; then
    run_step "workspace check" cargo check --workspace
else
    run_step "clippy" cargo clippy --workspace -- -D warnings
    run_step "workspace tests" cargo test --workspace
fi

echo "Quality gate passed ($mode)."
