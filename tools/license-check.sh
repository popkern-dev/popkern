#!/usr/bin/env bash
# PK-373 licensing gate, local mirror of the CI `licenses` job.
#   tools/license-check.sh               # deny allowlist + notices drift check
#   tools/license-check.sh --regenerate  # rewrite THIRD-PARTY-NOTICES.md
# Needs: cargo install cargo-deny; cargo install cargo-about --features cli
set -euo pipefail
cd "$(dirname "$0")/.."
cargo deny check licenses bans sources
cargo about generate --all-features --workspace about.hbs -o THIRD-PARTY-NOTICES.md
if [[ "${1:-}" == "--regenerate" ]]; then
  echo "THIRD-PARTY-NOTICES.md regenerated"
elif ! git diff --quiet -- THIRD-PARTY-NOTICES.md; then
  echo "THIRD-PARTY-NOTICES.md is stale — run tools/license-check.sh --regenerate and commit" >&2
  exit 1
else
  echo "license policy ok, notices current"
fi
