#!/usr/bin/env bash
# Prepare a Mac to build and publish popkern's macOS release archive
# (tools/dist-macos.sh). The macOS counterpart of
# tools/setup-windows-build.ps1. Safe to re-run: every step checks first
# and only installs what is missing.
#
#   tools/setup-macos-build.sh
#
# Installs, per user (no sudo, no Homebrew):
#   - rustup, with the toolchain rust-toolchain.toml pins as the default and
#     the x86_64/aarch64 Apple and wasm32 targets
#   - wasm-bindgen-cli at the version Cargo.lock locks (tools/build-live.sh)
#   - Node.js LTS under ~/.local/share/node, and wrangler (uploads to the
#     popkern-releases R2 bucket), both linked into ~/.local/bin
# and checks for the Xcode Command Line Tools (clang, ld, lipo, codesign),
# starting their installer when they are missing.
#
# It does not log wrangler in: run `wrangler login` once (a browser OAuth
# grant), or for unattended nightly builds set CLOUDFLARE_API_TOKEN (an R2
# edit token) and CLOUDFLARE_ACCOUNT_ID in the build user's environment.
set -euo pipefail
cd "$(dirname "$0")/.."

NODE_MAJOR=24 # the Active LTS line; wrangler needs 20 or later
BIN=$HOME/.local/bin

step() { printf '== %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }
export PATH="$HOME/.cargo/bin:$BIN:$PATH"
mkdir -p "$BIN"

[[ $(uname -s) == Darwin ]] || { echo "setup-macos-build.sh is for macOS" >&2; exit 1; }

step 'Xcode Command Line Tools (clang, ld, lipo, codesign)'
if xcode-select -p >/dev/null 2>&1 && have clang && have lipo; then
  echo "   present: $(xcode-select -p)"
else
  xcode-select --install || true
  echo "   finish the Command Line Tools installer that just opened, then re-run this script" >&2
  exit 1
fi

step 'Rust (rustup, pinned toolchain, release targets)'
if ! have rustup; then
  tmp=$(mktemp -d)
  curl -fsSL --proto '=https' --tlsv1.2 -o "$tmp/rustup-init.sh" https://sh.rustup.rs
  # --no-modify-path: the rc-file edit is shown below instead of made silently.
  sh "$tmp/rustup-init.sh" -y --default-toolchain none --profile minimal --no-modify-path
  rm -rf "$tmp"
fi
pinned=$(awk -F'"' '/^channel = /{print $2; exit}' rust-toolchain.toml)
# The pin as the default too, so cargo also works outside the checkout
# (cargo install below) without fetching a second toolchain.
rustup toolchain install "$pinned" --profile minimal --component rustfmt --component clippy
rustup default "$pinned"
rustup target add --toolchain "$pinned" x86_64-apple-darwin aarch64-apple-darwin wasm32-unknown-unknown

step 'wasm-bindgen-cli (the version Cargo.lock locks)'
locked=$(awk '/^name = "wasm-bindgen"$/{getline; sub(/version = "/,""); sub(/"/,""); print; exit}' Cargo.lock)
[[ -n $locked ]] || { echo "cannot read the wasm-bindgen version from Cargo.lock" >&2; exit 1; }
current=$(have wasm-bindgen && wasm-bindgen --version | awk '{print $2}' || true)
if [[ $current == "$locked" ]]; then
  echo "   present: wasm-bindgen $current"
else
  (cd "$HOME" && cargo install wasm-bindgen-cli --version "$locked" --locked --force)
fi

step "Node.js $NODE_MAJOR LTS and wrangler"
if ! have node || (( $(node -p 'process.versions.node.split(".")[0]') < 20 )); then
  case "$(uname -m)" in arm64) arch=arm64 ;; *) arch=x64 ;; esac
  dist=https://nodejs.org/dist/latest-v$NODE_MAJOR.x
  tmp=$(mktemp -d)
  curl -fsSL --proto '=https' -o "$tmp/SHASUMS256.txt" "$dist/SHASUMS256.txt"
  read -r sum file < <(grep -E " node-v[0-9.]+-darwin-$arch\.tar\.gz$" "$tmp/SHASUMS256.txt")
  [[ -n ${file:-} ]] || { echo "no darwin-$arch tarball in $dist/SHASUMS256.txt" >&2; exit 1; }
  curl -fsSL --proto '=https' -o "$tmp/$file" "$dist/$file"
  [[ $(shasum -a 256 "$tmp/$file" | cut -d' ' -f1) == "$sum" ]] || { echo "checksum mismatch for $file" >&2; exit 1; }
  rm -rf "$HOME/.local/share/node.new"
  mkdir -p "$HOME/.local/share/node.new"
  tar -xzf "$tmp/$file" -C "$HOME/.local/share/node.new" --strip-components 1
  rm -rf "$HOME/.local/share/node" "$tmp"
  mv "$HOME/.local/share/node.new" "$HOME/.local/share/node"
  for tool in node npm npx; do ln -sfn "$HOME/.local/share/node/bin/$tool" "$BIN/$tool"; done
fi
echo "   node $(node --version)"
if ! have wrangler; then
  # --prefix ~/.local: the wrangler link lands in ~/.local/bin, beside node.
  # npm 11 skips dependencies' install scripts unless allowed; wrangler's
  # esbuild and workerd fetch their platform binaries in theirs.
  npm install --global --prefix "$HOME/.local" --allow-scripts=esbuild,workerd wrangler
fi
echo "   wrangler $(wrangler --version)"

cat <<EOF

Ready. Make sure your shell finds the tools (once, in ~/.zshrc):
    export PATH="\$HOME/.cargo/bin:\$HOME/.local/bin:\$PATH"
Next, once:  wrangler login
Then, from the checkout:
    tools/dist-macos.sh                       # build + verify the archive
    tools/dist-macos.sh --upload [--latest]   # publish to dl.popkern.dev
EOF
