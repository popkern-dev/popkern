#!/usr/bin/env bash
# Attribute popkern live's wasm bundle size by crate and module
# (PK-658): what grew, measured in shipped code bytes, not source lines.
#
#   tools/live-size-attribution.sh               # the working tree
#   tools/live-size-attribution.sh <git-ref>     # the working tree vs <git-ref>
#
# Builds the `live` profile with symbol names kept (the shipped bundle
# strips them) into its own target directory, so the main build cache
# and web/live/pkg are untouched, then sums twiggy's per-function
# shallow sizes by crate and by the top two module levels of
# popkern_core. A <git-ref> is built from a temporary worktree (its
# crate may still be named popkern-playground, before PK-657). The
# "function names" subsection only exists in these unstripped builds:
# ignore its row when reading shipped growth.
#
# Tools: twiggy (cargo install twiggy --locked), the wasm32 target.
set -euo pipefail
cd "$(dirname "$0")/.."
command -v twiggy >/dev/null || { echo "missing twiggy: run  cargo install twiggy --locked" >&2; exit 1; }
WORK="$PWD/target/live-size-attribution"
mkdir -p "$WORK"

build() { # source-dir out-csv
  local dir="$1" out="$2" crate=popkern-live profile=live
  if [[ ! -d "$dir/crates/popkern-live" ]]; then crate=popkern-playground; profile=playground; fi
  local upper; upper=$(echo "$profile" | tr '[:lower:]' '[:upper:]')
  (cd "$dir" && env "CARGO_PROFILE_${upper}_STRIP=false" CARGO_TARGET_DIR="$WORK/target" \
    cargo build --quiet --profile "$profile" --target wasm32-unknown-unknown -p "$crate")
  twiggy top -n 100000000 -f csv "$WORK/target/wasm32-unknown-unknown/$profile/${crate//-/_}.wasm" > "$out"
}

build . "$WORK/now.csv"
BASE_CSV=""
if [[ $# -ge 1 ]]; then
  TREE="$WORK/base-tree"
  git worktree remove --force "$TREE" 2>/dev/null || true
  git worktree add --quiet --detach "$TREE" "$1"
  trap 'git worktree remove --force "$TREE"' EXIT
  build "$TREE" "$WORK/base.csv"
  BASE_CSV="$WORK/base.csv"
fi

python3 - "$WORK/now.csv" "$BASE_CSV" <<'EOF'
import collections, csv, re, sys

CRATE = re.compile(r"\b(popkern_\w+|rhai|faer\w*|serde_json|serde|core|alloc|std|hashbrown|"
                   r"smallvec|ttf_parser|wasm_bindgen|js_sys|dlmalloc|pulp|gemm\w*|rayon\w*|indexmap|glam)\b")
MODULE = re.compile(r"popkern_core(?:\[\w+\])?((?:::\w+){1,2})")

def sums(path):
    crates, modules = collections.Counter(), collections.Counter()
    if not path:
        return crates, modules
    for row in csv.DictReader(open(path)):
        name, size = row["Name"], int(row["ShallowSize"])
        if name.startswith('"function names"'):
            crates["(names section, stripped when shipped)"] += size
            continue
        m = CRATE.search(name)
        crate = m.group(1) if m else "other"
        crates["popkern_live" if crate == "popkern_playground" else crate] += size
        m = MODULE.search(name)
        if m:
            modules[m.group(1)] += size
    return crates, modules

def table(title, now, base, rows):
    print(f"\n{title}")
    print(f"{'':48}{'base':>12}{'now':>12}{'delta':>12}")
    keys = sorted(set(now) | set(base), key=lambda k: -(now[k] - base[k]) if base else -now[k])
    for k in keys[:rows]:
        print(f"{k:48}{base[k]:>12}{now[k]:>12}{now[k] - base[k]:>12}")

now_c, now_m = sums(sys.argv[1])
base_c, base_m = sums(sys.argv[2])
table("by crate (shallow code + data bytes)", now_c, base_c, 20)
table("popkern_core by module", now_m, base_m, 15)
print(f"\ntotal: base {sum(base_c.values())}, now {sum(now_c.values())}")
EOF
