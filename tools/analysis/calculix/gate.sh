#!/usr/bin/env bash
# PH-075: the CalculiX gate (gate.py, cases in gates.toml). Needs ccx
# 2.21, gmsh 4.12 with python3-gmsh and python3 >= 3.11. CalculiX's
# results are cached in the work dir (default target/calculix-gate);
# with them cached the gate takes a few minutes (ours only), and from
# scratch about 15 minutes more. Exits non-zero when a ratio leaves
# its tolerance. Extra arguments go to gate.py (--record,
# --tolerance-scale S).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
work=${WORK:-target/calculix-gate}
cargo build --release -p popkern-analysis-cli
python3 tools/analysis/calculix/gate.py target/release/popkern-analysis "$work" "$@"
