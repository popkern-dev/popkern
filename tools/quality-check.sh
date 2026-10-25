#!/usr/bin/env bash
# Full quality gate for popkern.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mode="${1:---full}"

case "$mode" in
    --full|--fast|--perf|--physics)
        ;;
    -h|--help)
        cat <<'USAGE'
Usage: tools/quality-check.sh [--full|--fast|--perf|--physics]

--full     Run fmt, semantic-name scan, structural tests, clippy, and the
           kernel suite (every crate's tests except the physics crates',
           which keep only the BRep adapter's unit tests as a smoke set).
--fast     Run fmt, semantic-name scan, structural tests, and cargo check.
--perf     Run the release-mode performance budget gate (PK-190): the
           wall-clock budget suites and the §17 cancellation-latency
           test — the same gate CI's perf job enforces.
--physics  Run the physics suite (PK-948): the analysis crates' tests
           and doctests and the solver's wasm32 clippy, in their own
           target dir (target/physics) so it can run beside the kernel
           suite without swapping its test binaries.
USAGE
        exit 0
        ;;
    *)
        echo "Usage: tools/quality-check.sh [--full|--fast|--perf|--physics]" >&2
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
        cargo test --release -p popkern-blend \
        --test perf_budget --test cancellation_latency --test perf_acceptance \
        -- --nocapture
    run_step "smooth-mesh perf budgets + thread determinism gate" \
        cargo test --release -p popkern-part \
        --test perf_budget --test thread_determinism_gate -- --nocapture
    run_step "seam certification perf budget (PK-562)" \
        cargo test --release -p popkern-core --test certify_perf_budget -- --nocapture
    run_step "KP20 cap assembly perf budget (PK-600)" \
        cargo test --release -p popkern-core --test kp20_cap_perf_budget -- --nocapture
    run_step "boolean perf budgets" \
        cargo test --release -p popkern-bool --test perf_budget -- --nocapture
    run_step "tessellation perf budgets" \
        cargo test --release -p popkern-mesh --test perf_budget -- --nocapture
    run_step "offset/shell perf budgets + determinism pins" \
        cargo test --release -p popkern-examples \
        --test perf_budget --test offset_determinism_gates -- --nocapture
    echo "Quality gate passed ($mode)."
    exit 0
fi

# The physics crates (PK-948): their suite runs on its own, and the
# kernel suite builds none of their test binaries but the adapter's.
physics_packages=(popkern-analysis popkern-analysis-cli popkern-analysis-setup popkern-analysis-fixtures popkern-analysis-wasm-bench)
physics_only=()
kernel_excludes=()
for package in "${physics_packages[@]}"; do
    physics_only+=(-p "$package")
    kernel_excludes+=(--exclude "$package")
done

if [[ "$mode" == "--physics" ]]; then
    export CARGO_TARGET_DIR="$ROOT/target/physics"
    if rustup target list --installed 2>/dev/null | grep -qx wasm32-unknown-unknown; then
        run_step "solver clippy on wasm32" \
            cargo clippy --target wasm32-unknown-unknown -p popkern-analysis -- -D warnings
    fi
    run_step "physics tests (nextest)" \
        cargo nextest run --profile physics "${physics_only[@]}" -p popkern-analysis-brep
    run_step "physics doc tests" \
        cargo test --doc "${physics_only[@]}" -p popkern-analysis-brep
    echo "Quality gate passed ($mode)."
    exit 0
fi

run_step "format" cargo fmt --check
run_step "semantic names" tools/check-semantic-names.sh
run_step "structural tests" cargo test -p popkern-structural-tests

if [[ "$mode" == "--fast" ]]; then
    run_step "workspace check" cargo check --workspace
else
    # PK-270: --all-targets matches CI's clippy leg — the bare
    # --workspace variant skipped tests/benches/examples, so local
    # and CI gated different target sets and shared no fingerprints.
    run_step "clippy" cargo clippy --workspace --all-targets -- -D warnings
    # PK-267: nextest schedules all test binaries into one parallel
    # run (cargo test runs them strictly sequentially) and isolates
    # each test in its own process. nextest does not run doctests, so
    # a --doc leg always follows.
    # PK-270: explicit --lib --bins --tests matches CI's test leg.
    if command -v cargo-nextest >/dev/null 2>&1; then
        run_step "workspace tests (nextest)" \
            cargo nextest run --workspace "${kernel_excludes[@]}" --lib --bins --tests
    else
        echo "note: cargo-nextest not installed (https://nexte.st) — falling back to serial cargo test" >&2
        run_step "workspace tests" \
            cargo test --workspace "${kernel_excludes[@]}" --lib --bins --tests
    fi
    run_step "doc tests" cargo test --doc --workspace "${kernel_excludes[@]}" --exclude popkern-analysis-brep
fi

echo "Quality gate passed ($mode)."
