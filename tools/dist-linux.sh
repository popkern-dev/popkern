#!/usr/bin/env bash
# Build popkern's Linux release archives, and optionally publish them to
# the release bucket behind https://dl.popkern.dev (what install.sh reads).
#
#   tools/dist-linux.sh                         # build both targets
#   tools/dist-linux.sh --version 0.1.0-rc.1    # stamp a release candidate
#   tools/dist-linux.sh --target x86_64-unknown-linux-musl
#   tools/dist-linux.sh --no-page-build         # reuse web/live/pkg as built
#   tools/dist-linux.sh --upload [--latest]     # publish; --latest moves `latest`
#
# The binaries are static musl builds: no libc, loader or distribution
# dependency, so one archive per CPU runs on any Linux. Each archive is
#
#   popkern-<version>-<target>/bin/popkern
#   popkern-<version>-<target>/share/popkern/live/   (popkern live's page)
#
# and `popkern serve` finds share/popkern/live beside its bin/ folder
# (popkern_cli::installed_page). The page ships without bakes: examples
# build on first open instead of loading pre-built models.
#
# Output: target/dist/v<version>/{<archives>, SHA256SUMS}. The bucket
# holds the same tree, plus `latest` (the version install.sh installs by
# default). Uploading needs wrangler, logged in to the account that owns
# the popkern-releases bucket.
set -euo pipefail
cd "$(dirname "$0")/.."

BUCKET=popkern-releases
VERSION=$(awk -F'"' '/^version = /{print $2; exit}' crates/popkern-cli/Cargo.toml)
TARGETS=()
PAGE_BUILD=1; UPLOAD=0; LATEST=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --version) VERSION=${2#v}; shift ;;
    --target) TARGETS+=("$2"); shift ;;
    --no-page-build) PAGE_BUILD=0 ;;
    --upload) UPLOAD=1 ;;
    --latest) LATEST=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
[[ ${#TARGETS[@]} -eq 0 ]] && TARGETS=(x86_64-unknown-linux-musl aarch64-unknown-linux-musl)
[[ $LATEST == 1 && $UPLOAD == 0 ]] && { echo "--latest needs --upload" >&2; exit 2; }
if [[ ! $VERSION =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.]+)?$ ]]; then
  echo "version '$VERSION' is not semver" >&2; exit 2
fi

for target in "${TARGETS[@]}"; do
  if ! rustup target list --installed | grep -qx "$target"; then
    echo "missing target: run  rustup target add $target" >&2; exit 1
  fi
done

# The musl targets link with their self-contained CRT objects, so the
# toolchain's own lld links both CPUs. The one C dependency (mimalloc)
# compiles with zig cc, a cross C compiler for every target in one
# binary (https://ziglang.org/download/).
if ! command -v zig >/dev/null; then
  echo "missing zig: install it from https://ziglang.org/download/ (0.14 or later)" >&2; exit 1
fi
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=rust-lld
export CARGO_TARGET_AARCH64_UNKNOWN_LINUX_MUSL_LINKER=rust-lld
export CC_x86_64_unknown_linux_musl=$PWD/tools/zig-cc/x86_64-linux-musl-cc
export CC_aarch64_unknown_linux_musl=$PWD/tools/zig-cc/aarch64-linux-musl-cc
export POPKERN_VERSION=$VERSION

[[ $PAGE_BUILD == 1 ]] && tools/build-live.sh --web
if [[ ! -f web/live/pkg/popkern_live_bg.wasm ]]; then
  echo "no wasm bundle at web/live/pkg/ (run without --no-page-build)" >&2; exit 1
fi

TARGET_DIR=${CARGO_TARGET_DIR:-target}
OUT="$TARGET_DIR/dist/v$VERSION"
rm -rf "$OUT"
mkdir -p "$OUT"
for target in "${TARGETS[@]}"; do
  cargo build --profile dist --target "$target" -p popkern-cli --bin popkern
  name="popkern-$VERSION-$target"
  stage="$OUT/$name"
  mkdir -p "$stage/bin" "$stage/share/popkern/live"
  cp "$TARGET_DIR/$target/dist/popkern" "$stage/bin/popkern"
  # The runtime page only (the deploy's set, without bakes).
  tar -C web/live \
    --exclude './test' --exclude './bakes' --exclude './core.*' \
    --exclude './package.json' --exclude './pkg/*.d.ts' --exclude './pkg/package.json' \
    -cf - . | tar -C "$stage/share/popkern/live" -xf -
  cp LICENSE* "$stage/" 2>/dev/null || true
  # Reproducible-ish archives: sorted names, fixed owner, no source mtimes.
  tar -C "$OUT" --sort=name --owner=0 --group=0 --numeric-owner \
    --mtime='@0' -czf "$OUT/$name.tar.gz" "$name"
  rm -rf "$stage"
  if [[ $target == x86_64-* ]]; then
    tmp=$(mktemp -d)
    tar -C "$tmp" -xzf "$OUT/$name.tar.gz"
    got=$("$tmp/$name/bin/popkern" --version)
    rm -rf "$tmp"
    [[ $got == "popkern $VERSION" ]] || { echo "built binary says '$got', not 'popkern $VERSION'" >&2; exit 1; }
  fi
  echo "built $OUT/$name.tar.gz ($(du -h "$OUT/$name.tar.gz" | cut -f1))"
done
(cd "$OUT" && sha256sum -- *.tar.gz > SHA256SUMS)
cat "$OUT/SHA256SUMS"

[[ $UPLOAD == 0 ]] && exit 0
put() { # <key> <file> <content-type> <cache-control>
  wrangler r2 object put "$BUCKET/$1" --remote --file "$2" \
    --content-type "$3" --cache-control "$4" >/dev/null
  echo "uploaded https://dl.popkern.dev/$1"
}
for file in "$OUT"/*.tar.gz; do
  put "v$VERSION/$(basename "$file")" "$file" application/gzip "public, max-age=31536000, immutable"
done
put "v$VERSION/SHA256SUMS" "$OUT/SHA256SUMS" "text/plain; charset=utf-8" "public, max-age=31536000, immutable"
if [[ $LATEST == 1 ]]; then
  printf '%s\n' "$VERSION" > "$OUT/latest"
  put latest "$OUT/latest" "text/plain; charset=utf-8" "public, max-age=60"
fi
