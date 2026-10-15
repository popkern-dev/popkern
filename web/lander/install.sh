#!/bin/sh
# Install popkern: the `popkern` command and popkern live's page.
#
#   curl -fsSL https://popkern.dev/install.sh | sh
#
# What it does, in order:
#   1. picks the build for this machine (Linux x86_64 or aarch64; static,
#      so any distribution and any libc);
#   2. downloads it and the release's SHA256SUMS, and refuses to install
#      anything whose checksum does not match;
#   3. unpacks it to $POPKERN_HOME/<version> (default ~/.local/share/popkern),
#      points $POPKERN_HOME/current at it, and links
#      $POPKERN_INSTALL_DIR/popkern (default ~/.local/bin) to the binary.
# Nothing outside those two folders is touched, and nothing runs as root.
#
# Settings (environment variables):
#   POPKERN_VERSION        a release to install, e.g. 0.1.0 (default: latest)
#   POPKERN_INSTALL_DIR    where the `popkern` link goes (default ~/.local/bin)
#   POPKERN_HOME           where releases unpack (default ~/.local/share/popkern)
#   POPKERN_DOWNLOAD_BASE  where releases are fetched from
#                          (default https://dl.popkern.dev)
#
# Uninstall: rm -rf ~/.local/share/popkern ~/.local/bin/popkern
set -eu

BASE=${POPKERN_DOWNLOAD_BASE:-https://dl.popkern.dev}
BASE=${BASE%/}
INSTALL_DIR=${POPKERN_INSTALL_DIR:-$HOME/.local/bin}
POPKERN_HOME=${POPKERN_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/popkern}

say() { printf 'popkern: %s\n' "$*"; }
fail() { printf 'popkern: %s\n' "$*" >&2; exit 1; }

# The body runs only once the whole script has arrived, so a download cut
# short mid-pipe cannot run half an installer.
main() {
  case "$(uname -s)" in
    Linux) ;;
    Darwin) fail "macOS builds are coming soon. For now, try popkern live in your browser: https://popkern.dev/live/" ;;
    *) fail "no build for $(uname -s) yet. Try popkern live in your browser: https://popkern.dev/live/" ;;
  esac
  case "$(uname -m)" in
    x86_64 | amd64) target=x86_64-unknown-linux-musl ;;
    aarch64 | arm64) target=aarch64-unknown-linux-musl ;;
    *) fail "no Linux build for $(uname -m) yet" ;;
  esac

  if command -v curl >/dev/null 2>&1; then
    fetch() { curl -fsSL --retry 3 --proto '=https' --tlsv1.2 -o "$2" "$1"; }
  elif command -v wget >/dev/null 2>&1; then
    fetch() { wget -q -O "$2" "$1"; }
  else
    fail "needs curl or wget"
  fi
  if command -v sha256sum >/dev/null 2>&1; then
    sha256() { sha256sum "$1" | cut -d' ' -f1; }
  elif command -v shasum >/dev/null 2>&1; then
    sha256() { shasum -a 256 "$1" | cut -d' ' -f1; }
  else
    fail "needs sha256sum or shasum to verify the download"
  fi
  command -v tar >/dev/null 2>&1 || fail "needs tar"

  tmp=$(mktemp -d 2>/dev/null || mktemp -d -t popkern)
  trap 'rm -rf "$tmp"' EXIT INT TERM

  version=${POPKERN_VERSION:-}
  if [ -z "$version" ]; then
    fetch "$BASE/latest" "$tmp/latest" || fail "cannot reach $BASE/latest"
    version=$(tr -d ' \r\n' <"$tmp/latest")
  fi
  version=${version#v}
  case "$version" in
    '' | *[!0-9A-Za-z.+-]*) fail "unexpected version '$version' from $BASE/latest" ;;
  esac

  archive="popkern-$version-$target.tar.gz"
  say "downloading popkern $version ($target)"
  fetch "$BASE/v$version/$archive" "$tmp/$archive" ||
    fail "no $archive in release $version ($BASE/v$version/)"
  fetch "$BASE/v$version/SHA256SUMS" "$tmp/SHA256SUMS" ||
    fail "release $version has no SHA256SUMS"

  want=$(awk -v f="$archive" '$2 == f || $2 == "*" f { print $1 }' "$tmp/SHA256SUMS")
  [ -n "$want" ] || fail "SHA256SUMS does not list $archive"
  got=$(sha256 "$tmp/$archive")
  [ "$got" = "$want" ] || fail "checksum mismatch for $archive (expected $want, got $got)"

  mkdir -p "$tmp/unpack"
  tar -xzf "$tmp/$archive" -C "$tmp/unpack"
  unpacked="$tmp/unpack/popkern-$version-$target"
  [ -x "$unpacked/bin/popkern" ] || fail "$archive has no bin/popkern"

  mkdir -p "$POPKERN_HOME" "$INSTALL_DIR"
  dest="$POPKERN_HOME/$version"
  rm -rf "$dest.new"
  mv "$unpacked" "$dest.new"
  rm -rf "$dest"
  mv "$dest.new" "$dest"
  ln -sfn "$version" "$POPKERN_HOME/current"
  ln -sf "$POPKERN_HOME/current/bin/popkern" "$INSTALL_DIR/popkern"

  installed=$("$INSTALL_DIR/popkern" --version 2>/dev/null) ||
    fail "installed, but $INSTALL_DIR/popkern does not run on this machine"
  say "installed $installed to $INSTALL_DIR/popkern"

  case ":$PATH:" in
    *":$INSTALL_DIR:"*) ;;
    *)
      say "$INSTALL_DIR is not on your PATH. Add it with:"
      printf '\n    export PATH="%s:$PATH"\n\n' "$INSTALL_DIR"
      say "(put that line in ~/.profile, ~/.bashrc or ~/.zshrc to keep it)"
      ;;
  esac
  cat <<'EOF'

Make your first part:

    popkern new my-parts/bracket
    popkern serve my-parts

then open http://127.0.0.1:7690 in your browser.
EOF
}

main "$@"
