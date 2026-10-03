#!/usr/bin/env bash
# Pre-build the playground's gallery examples (bakes): run every script
# in web/playground/examples.json natively with no budgets, certify its
# continuity, and write web/playground/bakes/ for the page to show
# instantly (crates/popkern-playground/src/bake.rs has the format).
#
#   tools/bake-playground.sh            # after tools/build-playground.sh --web
#   tools/bake-playground.sh --jobs 4   # parallel scripts (default: min(cores, 8))
#
# Bakes are keyed to the wasm bundle they ship beside (its SHA-256): a
# bake from the same bundle and the same script text is reused, so a
# rerun with nothing changed costs seconds; a new bundle rebakes all.
# build-playground.sh deletes bakes made for another bundle, so a page
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

WASM=web/playground/pkg/popkern_playground_bg.wasm
if [[ ! -f $WASM ]]; then
  echo "no wasm bundle at $WASM — run tools/build-playground.sh --web first" >&2
  exit 1
fi
KERNEL=$(sha256sum "$WASM" | cut -c1-64)
cargo build --profile playground -p popkern-playground --example bake_gallery
target/playground/examples/bake_gallery "$KERNEL" $JOBS
