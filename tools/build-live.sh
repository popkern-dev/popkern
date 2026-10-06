#!/usr/bin/env bash
# Build popkern live's wasm bundle (design §6.5, PK-410).
#
#   tools/build-live.sh            # web + node bundles, unoptimised
#   tools/build-live.sh --opt      # then wasm-opt (level $WASM_OPT_LEVEL, default -O3; measured slower — see below)
#   tools/build-live.sh --node     # node bundle only (the smoke tests)
#   tools/build-live.sh --web      # web bundle only (the page)
#
# Outputs:
#   web/live/pkg/            wasm-bindgen --target web   (the page loads this)
#   target/live/node/        wasm-bindgen --target nodejs (smoke tests)
#
# Tools: rustup target wasm32-unknown-unknown; wasm-bindgen-cli at the
# exact version of the wasm-bindgen crate in Cargo.lock (the ABI is
# version-locked and mismatches fail obscurely, so this script checks
# first); wasm-opt from binaryen for --opt. Each missing tool prints
# its install command and fails.
set -euo pipefail
cd "$(dirname "$0")/.."

OPT=0; WEB=1; NODE=1
# The feature set rustc targets (wasm32 baseline + these); NOT
# --all-features, which let binaryen emit post-MVP types node rejects
# ("invalid value type 0x0"). Measured 2026-09-02 (PK-414): every
# wasm-opt level ran the boss plate ~20% SLOWER under node (355 ms
# unoptimised vs 430-460 ms at -O2/-O3/-O4/-Oz) for 1-3% less gzip, so
# the shipped bundle is unoptimised and --opt stays an option.
WASM_OPT_FEATURES="--enable-bulk-memory --enable-sign-ext --enable-nontrapping-float-to-int --enable-mutable-globals --enable-reference-types --enable-multivalue --enable-simd"
for arg in "$@"; do
  case "$arg" in
    --opt) OPT=1 ;;
    --node) WEB=0; NODE=1 ;;
    --web) WEB=1; NODE=0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

LOCKED=$(awk '/^name = "wasm-bindgen"$/{getline; sub(/version = "/,""); sub(/"/,""); print; exit}' Cargo.lock)
if [[ -z "$LOCKED" ]]; then
  echo "cannot read the wasm-bindgen version from Cargo.lock" >&2; exit 1
fi
if ! rustup target list --installed 2>/dev/null | grep -q '^wasm32-unknown-unknown$'; then
  echo "missing target: run  rustup target add wasm32-unknown-unknown" >&2; exit 1
fi
if ! command -v wasm-bindgen >/dev/null; then
  echo "missing wasm-bindgen: run  cargo install wasm-bindgen-cli --version $LOCKED --locked" >&2; exit 1
fi
CLI=$(wasm-bindgen --version | awk '{print $2}')
if [[ "$CLI" != "$LOCKED" ]]; then
  echo "wasm-bindgen-cli $CLI does not match the locked crate $LOCKED: run  cargo install wasm-bindgen-cli --version $LOCKED --locked --force" >&2; exit 1
fi
if [[ $OPT == 1 ]] && ! command -v wasm-opt >/dev/null; then
  echo "missing wasm-opt: install binaryen (apt install binaryen) or  cargo install wasm-opt" >&2; exit 1
fi

cargo build --profile live --target wasm32-unknown-unknown -p popkern-live
# Cargo honours CARGO_TARGET_DIR; so must the packaging step, or a
# relocated build silently re-packages a stale ./target wasm (measured
# 2026-10-06: two gallery_parity runs compared new native code against
# a days-old bundle and reported phantom native/wasm divergence).
WASM=${CARGO_TARGET_DIR:-target}/wasm32-unknown-unknown/live/popkern_live.wasm

emit() { # target out-dir
  mkdir -p "$2"
  wasm-bindgen --target "$1" --out-dir "$2" --out-name popkern_live "$WASM"
  if [[ $OPT == 1 ]]; then
    # Rust targets the wasm32 baseline plus bulk memory, sign extension,
    # non-trapping float-to-int, mutable globals, reference types and
    # multi-value by default; binaryen must be told so or it rejects
    # the module ("Bulk memory operation (bulk memory is disabled)").
    wasm-opt "${WASM_OPT_LEVEL:--O3}" $WASM_OPT_FEATURES \
      -o "$2/popkern_live_bg.wasm" "$2/popkern_live_bg.wasm"
  fi
}
[[ $WEB == 1 ]] && emit web web/live/pkg
# Bakes (tools/bake-live.sh) answer for one bundle: drop any made
# for another, so the page never shows a model this kernel did not build.
if [[ $WEB == 1 && -f web/live/bakes/index.json ]]; then
  KERNEL=$(sha256sum web/live/pkg/popkern_live_bg.wasm | cut -c1-64)
  if ! grep -q "\"kernel\": \"$KERNEL\"" web/live/bakes/index.json; then
    rm -rf web/live/bakes
    echo "bakes: removed (made for another bundle; run tools/bake-live.sh)"
  fi
fi
[[ $NODE == 1 ]] && emit nodejs target/live/node

report() { # file
  local raw gz
  raw=$(stat -c %s "$1")
  gz=$(gzip -9 -c "$1" | wc -c)
  awk -v f="$1" -v raw="$raw" -v gz="$gz" \
    'BEGIN { printf "%s: %d bytes raw, %d bytes gzip (%.2f MB / %.2f MB)\n", f, raw, gz, raw/1048576, gz/1048576 }'
}
[[ $WEB == 1 ]] && report web/live/pkg/popkern_live_bg.wasm
[[ $NODE == 1 ]] && report target/live/node/popkern_live_bg.wasm
[[ $OPT == 1 ]] && echo "wasm-opt: ${WASM_OPT_LEVEL:--O3}" || echo "wasm-opt: not run (the shipped bundle; --opt measured slower)"
