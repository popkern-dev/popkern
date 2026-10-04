#!/usr/bin/env bash
# Pre-build popkern live's gallery examples (bakes): run every script
# in web/live/examples.json natively with no budgets, certify its
# continuity, and write web/live/bakes/ for the page to show
# instantly (crates/popkern-live/src/bake.rs has the format).
#
#   tools/bake-live.sh            # after tools/build-live.sh --web
#   tools/bake-live.sh --jobs 4   # parallel scripts (default: min(cores, 8))
#
# Bakes are keyed to the wasm bundle they ship beside (its SHA-256): a
# bake from the same bundle and the same script text is reused, so a
# rerun with nothing changed costs seconds; a new bundle rebakes all.
# build-live.sh deletes bakes made for another bundle, so a page
# never shows a model its own kernel did not build.
set -euo pipefail
cd "$(dirname "$0")/.."

JOBS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs) JOBS="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

WASM=web/live/pkg/popkern_live_bg.wasm
if [[ ! -f $WASM ]]; then
  echo "no wasm bundle at $WASM — run tools/build-live.sh --web first" >&2
  exit 1
fi
KERNEL=$(sha256sum "$WASM" | cut -c1-64)
cargo build --profile live -p popkern-live --example bake_gallery
target/live/examples/bake_gallery "$KERNEL" $JOBS
