#!/usr/bin/env bash
# Full quality gate for popkern.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mode="${1:---full}"

case "$mode" in
    --full|--fast)
        ;;
    -h|--help)
        cat <<'USAGE'
Usage: tools/quality-check.sh [--full|--fast]

--full  Run fmt, semantic-name scan, structural tests, clippy, and all tests.
--fast  Run fmt, semantic-name scan, structural tests, and cargo check.
USAGE
        exit 0
        ;;
    *)
        echo "Usage: tools/quality-check.sh [--full|--fast]" >&2
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
