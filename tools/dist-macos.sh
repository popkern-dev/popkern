#!/usr/bin/env bash
# Build popkern's macOS release archive, and optionally publish it to the
# release bucket behind https://dl.popkern.dev (what install.sh reads).
# The macOS counterpart of tools/dist-linux.sh and tools/dist-windows.sh;
# runs on a Mac (tools/setup-macos-build.sh installs what it needs).
#
#   tools/dist-macos.sh                         # build the universal archive
#   tools/dist-macos.sh --version 0.1.0-rc.1    # stamp a release candidate
#   tools/dist-macos.sh --no-page-build         # reuse web/live/pkg as built
#   tools/dist-macos.sh --upload [--latest]     # publish; --latest moves `latest`
#
# One universal binary: the x86_64-apple-darwin and aarch64-apple-darwin
# builds joined with lipo, so a single archive runs natively on Intel and
# Apple silicon Macs, macOS $MACOSX_DEPLOYMENT_TARGET or later. It links
# only system libraries (libSystem and frameworks under /System), which
# every macOS ships. The archive is
#
#   popkern-<version>-universal-apple-darwin/bin/popkern
#   popkern-<version>-universal-apple-darwin/share/popkern/live/
#
# the same layout as Linux, so `popkern serve` finds share/popkern/live
# beside its bin/ folder (popkern_cli::installed_page).
#
# The binary is ad-hoc signed (Apple silicon refuses to run unsigned
# code), not Developer ID signed or notarized.
# install.sh downloads with curl, which sets no quarantine attribute, so
# Gatekeeper does not check it; a browser download of the archive would be
# quarantined and blocked.
#
# Output: target/dist/v<version>/{<archive>, SHA256SUMS}. Uploading merges
# these checksums into the release's published SHA256SUMS
# (tools/dist-publish.sh), so Linux, Windows and macOS can publish one
# release from different machines.
set -euo pipefail
cd "$(dirname "$0")/.."

VERSION=$(awk -F'"' '/^version = /{print $2; exit}' crates/popkern-cli/Cargo.toml)
TARGETS=(x86_64-apple-darwin aarch64-apple-darwin)
UNIVERSAL=universal-apple-darwin
PAGE_BUILD=1; UPLOAD=0; LATEST=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --version) VERSION=${2#v}; shift ;;
    --no-page-build) PAGE_BUILD=0 ;;
    --upload) UPLOAD=1 ;;
    --latest) LATEST=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
[[ $LATEST == 1 && $UPLOAD == 0 ]] && { echo "--latest needs --upload" >&2; exit 2; }
if [[ ! $VERSION =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.]+)?$ ]]; then
  echo "version '$VERSION' is not semver" >&2; exit 2
fi

[[ $(uname -s) == Darwin ]] || { echo "dist-macos.sh links with Apple's toolchain: run it on a Mac" >&2; exit 1; }
for tool in lipo otool codesign; do
  command -v "$tool" >/dev/null || { echo "missing $tool: run  xcode-select --install" >&2; exit 1; }
done
for target in "${TARGETS[@]}"; do
  if ! rustup target list --installed | grep -qx "$target"; then
    echo "missing target: run  rustup target add $target" >&2; exit 1
  fi
done

# 11.0 (Big Sur) is the first macOS on Apple silicon, so both slices
# share one floor. rustc and the cc crate (mimalloc is musl-only, but any
# future C dependency) both read this variable.
export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
export POPKERN_VERSION=$VERSION
# No AppleDouble (._*) files for extended attributes in the archives.
export COPYFILE_DISABLE=1

[[ $PAGE_BUILD == 1 ]] && tools/build-live.sh --web
if [[ ! -f web/live/pkg/popkern_live_bg.wasm ]]; then
  echo "no wasm bundle at web/live/pkg/ (run without --no-page-build)" >&2; exit 1
fi

TARGET_DIR=${CARGO_TARGET_DIR:-target}
OUT="$TARGET_DIR/dist/v$VERSION"
rm -rf "$OUT"
mkdir -p "$OUT"
slices=()
for target in "${TARGETS[@]}"; do
  cargo build --profile dist --target "$target" -p popkern-cli --bin popkern
  slices+=("$TARGET_DIR/$target/dist/popkern")
done

name="popkern-$VERSION-$UNIVERSAL"
stage="$OUT/$name"
mkdir -p "$stage/bin" "$stage/share/popkern/live"
lipo -create -output "$stage/bin/popkern" "${slices[@]}"
bin="$stage/bin/popkern"
# ld signs only the arm64 slice (and leaves x86_64 unsigned), so sign the
# joined binary: an ad-hoc signature (identity "-") on every slice.
codesign --force --sign - "$bin"

# What a Mac without developer tools needs from the binary: both CPUs,
# the deployment floor, system libraries only, and a valid signature on
# every slice.
archs=$(lipo -archs "$bin")
[[ $archs == "x86_64 arm64" ]] || { echo "universal binary has slices '$archs', not 'x86_64 arm64'" >&2; exit 1; }
for arch in x86_64 arm64; do
  minos=$(otool -arch "$arch" -l "$bin" | awk '/LC_BUILD_VERSION/{f=1} f && $1 == "minos" {print $2; exit}')
  [[ $minos == "$MACOSX_DEPLOYMENT_TARGET" ]] ||
    { echo "$arch slice targets macOS '$minos', not $MACOSX_DEPLOYMENT_TARGET" >&2; exit 1; }
  if foreign=$(otool -arch "$arch" -L "$bin" | tail -n +2 | awk '{print $1}' |
      grep -vE '^/usr/lib/|^/System/Library/') && [[ -n $foreign ]]; then
    echo "$arch slice links libraries macOS does not ship: $(echo $foreign)" >&2; exit 1
  fi
done
codesign --verify --strict --all-architectures "$bin" ||
  { echo "universal binary's signature does not verify" >&2; exit 1; }

# The runtime page only (the deploy's set, without bakes), as on Linux.
tar -C web/live \
  --exclude './test' --exclude './bakes' --exclude './core.*' \
  --exclude './package.json' --exclude './pkg/*.d.ts' --exclude './pkg/package.json' \
  -cf - . | tar -C "$stage/share/popkern/live" -xf -
cp LICENSE* "$stage/" 2>/dev/null || true
# Reproducible-ish archive with macOS's bsdtar (no --sort): sorted names
# fed one by one, fixed owner and mtime, no extended attributes or ACLs,
# and gzip without its own name and timestamp.
find "$stage" -exec touch -h -t 197001010000 {} +
(cd "$OUT" && find "$name" | LC_ALL=C sort |
  TZ=UTC tar -cf - -n -T - --format ustar --uid 0 --gid 0 --uname '' --gname '' \
    --no-xattrs --no-acls --no-fflags --no-mac-metadata |
  gzip -n -9 > "$name.tar.gz")
rm -rf "$stage"

# Run what was packed, on every CPU this Mac can run: the version stamp
# has to survive the build, lipo and the archive. The x86_64 slice runs
# only where Rosetta 2 is installed.
tmp=$(mktemp -d)
tar -C "$tmp" -xzf "$OUT/$name.tar.gz"
for arch in arm64 x86_64; do
  if [[ $arch == x86_64 && $(uname -m) == arm64 ]] && ! arch -x86_64 /bin/bash -c : 2>/dev/null; then
    echo "  skip running the x86_64 slice (no Rosetta 2: softwareupdate --install-rosetta)"
    continue
  fi
  got=$(arch "-$arch" "$tmp/$name/bin/popkern" --version)
  [[ $got == "popkern $VERSION" ]] ||
    { rm -rf "$tmp"; echo "$arch slice says '$got', not 'popkern $VERSION'" >&2; exit 1; }
  echo "  ok $arch slice runs: $got"
done
rm -rf "$tmp"
echo "built $OUT/$name.tar.gz ($(du -h "$OUT/$name.tar.gz" | cut -f1))"

# Text-mode lines ("<hash>  <file>"), the form install.sh and install.ps1 parse.
(cd "$OUT" && shasum -a 256 -- *.tar.gz > SHA256SUMS)
cat "$OUT/SHA256SUMS"

[[ $UPLOAD == 0 ]] && exit 0
source tools/dist-publish.sh
dist_publish_release "$OUT" "$VERSION" "$LATEST" "$OUT"/*.tar.gz
