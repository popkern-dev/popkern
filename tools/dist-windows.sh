#!/usr/bin/env bash
# Build popkern's Windows release archives, and optionally publish them to
# the release bucket behind https://dl.popkern.dev (what install.ps1 reads).
# The Windows counterpart of tools/dist-linux.sh; runs in Git Bash on a
# Windows machine (tools/setup-windows-build.ps1 installs what it needs).
#
#   tools/dist-windows.sh                         # build both targets
#   tools/dist-windows.sh --version 0.1.0-rc.1    # stamp a release candidate
#   tools/dist-windows.sh --target x86_64-pc-windows-msvc
#   tools/dist-windows.sh --no-page-build         # reuse web/live/pkg as built
#   tools/dist-windows.sh --upload [--latest]     # publish; --latest moves `latest`
#
# The binaries are MSVC builds with the C runtime linked statically
# (+crt-static): they import only DLLs every Windows ships, so one archive
# per CPU runs on any Windows 10 or 11 without a Visual C++ redistributable.
# Each archive is
#
#   popkern-<version>-<target>/bin/popkern.exe
#   popkern-<version>-<target>/share/popkern/live/   (popkern live's page)
#
# the same layout as Linux, so `popkern serve` finds share/popkern/live
# beside its bin/ folder (popkern_cli::installed_page).
#
# Output: target/dist/v<version>/{<archives>, SHA256SUMS}. Uploading merges
# these checksums into the release's published SHA256SUMS
# (tools/dist-publish.sh), so Linux and Windows can publish one release
# from different machines.
set -euo pipefail
cd "$(dirname "$0")/.."

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
[[ ${#TARGETS[@]} -eq 0 ]] && TARGETS=(x86_64-pc-windows-msvc aarch64-pc-windows-msvc)
[[ $LATEST == 1 && $UPLOAD == 0 ]] && { echo "--latest needs --upload" >&2; exit 2; }
if [[ ! $VERSION =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.]+)?$ ]]; then
  echo "version '$VERSION' is not semver" >&2; exit 2
fi

case "$(uname -s)" in
  MINGW* | MSYS* | CYGWIN*) ;;
  *) echo "dist-windows.sh builds with MSVC: run it in Git Bash on Windows" >&2; exit 1 ;;
esac
for target in "${TARGETS[@]}"; do
  if ! rustup target list --installed | grep -qx "$target"; then
    echo "missing target: run  rustup target add $target" >&2; exit 1
  fi
done
# Windows' own bsdtar writes the zip archives (Git Bash's GNU tar cannot).
ZIP_TAR=$(cygpath -u "${SYSTEMROOT:-C:\\Windows}")/System32/tar.exe
[[ -x $ZIP_TAR ]] || { echo "missing $ZIP_TAR (Windows 10 1803 or later ships it)" >&2; exit 1; }

# Static CRT for the release targets only: build-live.sh's wasm build and
# the host's build scripts keep their defaults.
export CARGO_TARGET_X86_64_PC_WINDOWS_MSVC_RUSTFLAGS="-C target-feature=+crt-static"
export CARGO_TARGET_AARCH64_PC_WINDOWS_MSVC_RUSTFLAGS="-C target-feature=+crt-static"
export POPKERN_VERSION=$VERSION

[[ $PAGE_BUILD == 1 ]] && tools/build-live.sh --web
if [[ ! -f web/live/pkg/popkern_live_bg.wasm ]]; then
  echo "no wasm bundle at web/live/pkg/ (run without --no-page-build)" >&2; exit 1
fi

HOST=$(rustc -vV | awk '/^host: /{print $2}')
TARGET_DIR=${CARGO_TARGET_DIR:-target}
OUT="$TARGET_DIR/dist/v$VERSION"
rm -rf "$OUT"
mkdir -p "$OUT"
for target in "${TARGETS[@]}"; do
  cargo build --profile dist --target "$target" -p popkern-cli --bin popkern
  name="popkern-$VERSION-$target"
  stage="$OUT/$name"
  mkdir -p "$stage/bin" "$stage/share/popkern/live"
  cp "$TARGET_DIR/$target/dist/popkern.exe" "$stage/bin/popkern.exe"
  # The static CRT is what lets the binary run without a redistributable:
  # refuse a build that still imports the VC++ runtime or the UCRT DLLs
  # (their names appear verbatim in the PE import table).
  if crt=$(grep -aoiE 'vcruntime[0-9_]*\.dll|msvcp[0-9_]*\.dll|api-ms-win-crt-[a-z0-9-]*\.dll' "$stage/bin/popkern.exe" | sort -u) && [[ -n $crt ]]; then
    echo "$target binary imports the dynamic C runtime ($(echo $crt)); +crt-static did not apply" >&2; exit 1
  fi
  # The runtime page only (the deploy's set, without bakes), as on Linux.
  tar -C web/live \
    --exclude './test' --exclude './bakes' --exclude './core.*' \
    --exclude './package.json' --exclude './pkg/*.d.ts' --exclude './pkg/package.json' \
    -cf - . | tar -C "$stage/share/popkern/live" -xf -
  cp LICENSE* "$stage/" 2>/dev/null || true
  # Fixed timestamps so a rebuild of the same sources zips the same
  # entries (zip stores local DOS time from 1980, so not the epoch).
  find "$stage" -exec touch -d '1980-01-02T00:00:00' {} +
  (cd "$OUT" && "$ZIP_TAR" -a -cf "$name.zip" "$name")
  rm -rf "$stage"
  # Run what was packed, when this machine's CPU can: the version stamp
  # has to survive the build and the archive. (The installed-page lookup
  # is popkern-cli's version_and_installed_page_layout test.)
  if [[ $target == "$HOST" ]]; then
    tmp=$(mktemp -d)
    "$ZIP_TAR" -xf "$OUT/$name.zip" -C "$tmp"
    got=$("$tmp/$name/bin/popkern.exe" --version | tr -d '\r')
    rm -rf "$tmp"
    [[ $got == "popkern $VERSION" ]] || { echo "built binary says '$got', not 'popkern $VERSION'" >&2; exit 1; }
  fi
  echo "built $OUT/$name.zip ($(du -h "$OUT/$name.zip" | cut -f1))"
done
# Text-mode lines ("<hash>  <file>"), the form install.sh and install.ps1 parse.
(cd "$OUT" && sha256sum -- *.zip | sed 's/ \*/  /' > SHA256SUMS)
cat "$OUT/SHA256SUMS"

[[ $UPLOAD == 0 ]] && exit 0
source tools/dist-publish.sh
dist_publish_release "$OUT" "$VERSION" "$LATEST" "$OUT"/*.zip
