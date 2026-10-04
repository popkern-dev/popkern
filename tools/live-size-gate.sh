#!/usr/bin/env bash
# Bundle-size ratchet for popkern live (PK-414).
#
#   tools/live-size-gate.sh              # check against the baseline
#   tools/live-size-gate.sh --update     # rewrite the baseline from the current bundle
#   (attribute first: tools/live-size-attribution.sh <last-baseline-commit>)
#
# The baseline (tools/live-size-baseline.txt) holds the measured
# raw and gzip byte counts of web/live/pkg/popkern_live_bg.wasm
# as shipped (unoptimised: wasm-opt measured ~20% slower for 1-3% less
# gzip, PK-414), plus the ceilings the gate enforces (baseline + 5%).
# A bundle that grows past a ceiling fails; a smaller bundle passes and
# the baseline is lowered only by an explicit --update, so every size
# change is a reviewed diff (the file-size ratchet's discipline).
set -euo pipefail
cd "$(dirname "$0")/.."
WASM=web/live/pkg/popkern_live_bg.wasm
BASELINE=tools/live-size-baseline.txt
if [[ ! -f "$WASM" ]]; then
  echo "no bundle at $WASM: run tools/build-live.sh --web first" >&2; exit 1
fi
raw=$(stat -c %s "$WASM")
gz=$(gzip -9 -c "$WASM" | wc -c)
if [[ "${1:-}" == "--update" ]]; then
  {
    echo "# popkern live bundle-size baseline (tools/live-size-gate.sh)."
    echo "# Measured on the shipped (unoptimised) bundle; ceilings are baseline + 5%. Update with --update."
    echo "raw_bytes=$raw"
    echo "gzip_bytes=$gz"
    echo "raw_ceiling=$(( raw + raw / 20 ))"
    echo "gzip_ceiling=$(( gz + gz / 20 ))"
  } > "$BASELINE"
  echo "baseline updated: raw $raw, gzip $gz"
  exit 0
fi
if [[ ! -f "$BASELINE" ]]; then
  echo "no baseline at $BASELINE: run tools/live-size-gate.sh --update" >&2; exit 1
fi
# shellcheck disable=SC1090
source <(grep -E '^[a-z_]+=' "$BASELINE")
status=0
printf 'bundle: raw %d bytes (baseline %d, ceiling %d); gzip %d bytes (baseline %d, ceiling %d)\n' \
  "$raw" "$raw_bytes" "$raw_ceiling" "$gz" "$gzip_bytes" "$gzip_ceiling"
if (( raw > raw_ceiling )); then echo "FAIL: raw bundle exceeds its ceiling" >&2; status=1; fi
if (( gz > gzip_ceiling )); then echo "FAIL: gzip bundle exceeds its ceiling" >&2; status=1; fi
if (( status == 0 )); then
  echo "bundle size within the ratchet"
else
  echo "attribute the growth before re-baselining: tools/live-size-attribution.sh <commit of the last baseline>" >&2
fi
exit $status
