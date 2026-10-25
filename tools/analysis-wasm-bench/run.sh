#!/usr/bin/env bash
# DT1 (native = wasm bits, in node) and E3 (solve time in headless
# Chrome). Needs the wasm32-unknown-unknown target, wasm-bindgen 0.2.114,
# node and google-chrome.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
# Headless Chrome can trap on exit after writing its output; never leave
# its core dumps in the repository.
ulimit -c 0
out=target/wasm-bench
target=${CARGO_TARGET_DIR:-target}
cargo test --release -p popkern-analysis-wasm-bench --test dt1 -- --include-ignored
cargo build --release -p popkern-analysis-wasm-bench --lib --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir "$out/pkg" \
  "$target"/wasm32-unknown-unknown/release/popkern_analysis_wasm_bench.wasm
node tools/analysis-wasm-bench/dt1.mjs
if [[ "${1:-}" == "--e3" ]]; then
  wasm-bindgen --target no-modules --out-dir "$out/pkg-classic" \
    "$target"/wasm32-unknown-unknown/release/popkern_analysis_wasm_bench.wasm
  budgets=(10000 15000 20000 30000 100000)
  cargo test --release -p popkern-analysis-wasm-bench --test dt1 native_e3_timing \
    -- --ignored --nocapture 2>&1 | grep '^native' 
  python3 tools/analysis-wasm-bench/e3.py "$out" "${budgets[@]}"
  for b in "${budgets[@]}"; do
    profile=$(mktemp -d)
    google-chrome --headless=new --disable-gpu --no-sandbox \
      --user-data-dir="$profile" --dump-dom "file://$PWD/$out/e3-$b.html" 2>/dev/null \
      | sed -n 's:.*<pre id="out">\(.*\)</pre>.*:\1:p'
    rm -rf "$profile"
  done
fi
