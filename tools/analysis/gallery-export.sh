#!/usr/bin/env bash
# Exports every popkern live gallery script to STEP with this checkout's
# popkern CLI, into target/gallery-export. The 23 parts up to 1.3 MB
# are frozen in tests/fixtures/analysis/gallery (copied from here when
# the fixtures are regenerated, docs/architecture/analysis-integration.md);
# the larger ones run from this export in the ignored extended_corpus
# test.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
cargo build --release -p popkern-cli
mkdir -p target/gallery-export
for script in web/live/examples/*.rhai; do
  name=$(basename "$script" .rhai)
  target/release/popkern run "$script" --step "target/gallery-export/$name.step" \
    > "target/gallery-export/$name.log" 2>&1 || echo "$name: export failed"
done
