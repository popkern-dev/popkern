#!/usr/bin/env bash
# PH-050: popkern-analysis against CalculiX on the cases in cases.toml.
# Needs ccx (calculix-ccx 2.21), gmsh 4.12 with its Python module
# (python3-gmsh) and python3 >= 3.11. CalculiX results are cached in the
# work dir (default target/calculix); delete a mesh's dir to re-solve.
# The finest meshes take up to half an hour and 14 GB each.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
work=${WORK:-target/calculix}
cargo build --release -p popkern-analysis-cli
python3 tools/analysis/calculix/compare.py target/release/popkern-analysis "$work" "$@"
