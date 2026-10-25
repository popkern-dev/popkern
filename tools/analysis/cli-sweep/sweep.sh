#!/usr/bin/env bash
# PH-105: the CLI robustness sweep (sweep.py) over the frozen gallery.
# Extra arguments go to sweep.py
# (--cap SECONDS).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
work=${WORK:-target/cli-sweep}
cargo build --release -p popkern-analysis-cli
python3 tools/analysis/cli-sweep/sweep.py target/release/popkern-analysis "$work" "$@"
