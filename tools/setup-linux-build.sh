#!/usr/bin/env bash
# Prepare a Linux machine to build and publish popkern's Linux release
# archives (tools/dist-linux.sh). The Linux counterpart of
# tools/setup-macos-build.sh. Safe to re-run: every step checks first and
# only installs what is missing.
#
#   tools/setup-linux-build.sh
#
# Installs, per user (no sudo, no package manager):
#   - rustup, with the toolchain rust-toolchain.toml pins as the default and
#     the x86_64/aarch64 musl and wasm32 targets
#   - wasm-bindgen-cli at the version Cargo.lock locks (tools/build-live.sh)
#   - zig $ZIG_VERSION under ~/.local/share/zig, the cross C/C++ compiler
#     and linker for both musl targets (tools/zig-cc), checksum-pinned
#   - Node.js LTS under ~/.local/share/node, and wrangler (uploads to the
#     popkern-releases R2 bucket), all linked into ~/.local/bin
# and checks for the host tools it cannot install without root: a C
# compiler (cargo links build scripts with cc), git, curl, xz, and
# readelf (dist-linux.sh's static-PIE gate), naming the apt package.
#
# It does not log wrangler in: run `wrangler login` once (a browser OAuth
# grant; over SSH, forward its callback port first:
# ssh -L 8976:localhost:8976 <this machine>), or for unattended nightly
# builds set CLOUDFLARE_API_TOKEN (an R2 edit token) and
# CLOUDFLARE_ACCOUNT_ID in the build user's environment.
set -euo pipefail
cd "$(dirname "$0")/.."

NODE_MAJOR=24 # the LTS line setup-macos-build.sh installs; wrangler needs 20 or later
ZIG_VERSION=0.17.0
BIN=$HOME/.local/bin

step() { printf '== %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }
export PATH="$HOME/.cargo/bin:$BIN:$PATH"
mkdir -p "$BIN"

[[ $(uname -s) == Linux ]] || { echo "setup-linux-build.sh is for Linux" >&2; exit 1; }
case "$(uname -m)" in
  x86_64) ARCH=x86_64; NODE_ARCH=x64
    ZIG_SHA256=1cbe9df9f27e6b78d14ccbca43b6703a404ef79ef1c463de901d7f088d4e2026 ;;
  aarch64 | arm64) ARCH=aarch64; NODE_ARCH=arm64
    ZIG_SHA256=9e8d11661d4ae3bd57702a3832781e23ad151dde5798e16a5ccd503f65234ff8 ;;
  *) echo "no Linux build setup for $(uname -m)" >&2; exit 1 ;;
esac

# install_tree <dir> <tarball>: unpack into <dir> by way of <dir>.new, so
# an interrupted run never leaves a half-replaced tool behind.
install_tree() {
  rm -rf "$1.new"
  mkdir -p "$1.new"
  tar -xf "$2" -C "$1.new" --strip-components 1
  rm -rf "$1"
  mv "$1.new" "$1"
}

step 'Host tools (cc, git, curl, xz, readelf)'
missing=()
for tool in cc:build-essential git:git curl:curl xz:xz-utils readelf:binutils; do
  have "${tool%%:*}" || missing+=("${tool#*:}")
done
if (( ${#missing[@]} )); then
  echo "   missing: install them with  sudo apt install ${missing[*]}" >&2; exit 1
fi
echo "   present"

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
rustup target add --toolchain "$pinned" \
  x86_64-unknown-linux-musl aarch64-unknown-linux-musl wasm32-unknown-unknown

step 'wasm-bindgen-cli (the version Cargo.lock locks)'
locked=$(awk '/^name = "wasm-bindgen"$/{getline; sub(/version = "/,""); sub(/"/,""); print; exit}' Cargo.lock)
[[ -n $locked ]] || { echo "cannot read the wasm-bindgen version from Cargo.lock" >&2; exit 1; }
current=$(have wasm-bindgen && wasm-bindgen --version | awk '{print $2}' || true)
if [[ $current == "$locked" ]]; then
  echo "   present: wasm-bindgen $current"
else
  (cd "$HOME" && cargo install wasm-bindgen-cli --version "$locked" --locked --force)
fi

step "zig $ZIG_VERSION (cross C/C++ compiler and linker)"
if [[ $(have zig && zig version || true) != "$ZIG_VERSION" ]]; then
  file=zig-$ARCH-linux-$ZIG_VERSION.tar.xz
  tmp=$(mktemp -d)
  curl -fsSL --proto '=https' -o "$tmp/$file" "https://ziglang.org/download/$ZIG_VERSION/$file"
  [[ $(sha256sum "$tmp/$file" | cut -d' ' -f1) == "$ZIG_SHA256" ]] || { echo "checksum mismatch for $file" >&2; exit 1; }
  install_tree "$HOME/.local/share/zig" "$tmp/$file"
  rm -rf "$tmp"
  ln -sfn "$HOME/.local/share/zig/zig" "$BIN/zig"
fi
echo "   zig $(zig version)"

step "Node.js $NODE_MAJOR LTS and wrangler"
if ! have node || (( $(node -p 'process.versions.node.split(".")[0]') < 20 )); then
  dist=https://nodejs.org/dist/latest-v$NODE_MAJOR.x
  tmp=$(mktemp -d)
  curl -fsSL --proto '=https' -o "$tmp/SHASUMS256.txt" "$dist/SHASUMS256.txt"
  read -r sum file < <(grep -E " node-v[0-9.]+-linux-$NODE_ARCH\.tar\.xz$" "$tmp/SHASUMS256.txt")
  [[ -n ${file:-} ]] || { echo "no linux-$NODE_ARCH tarball in $dist/SHASUMS256.txt" >&2; exit 1; }
  curl -fsSL --proto '=https' -o "$tmp/$file" "$dist/$file"
  [[ $(sha256sum "$tmp/$file" | cut -d' ' -f1) == "$sum" ]] || { echo "checksum mismatch for $file" >&2; exit 1; }
  install_tree "$HOME/.local/share/node" "$tmp/$file"
  rm -rf "$tmp"
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

Ready. Make sure your shell finds the tools (once, in ~/.bashrc):
    export PATH="\$HOME/.cargo/bin:\$HOME/.local/bin:\$PATH"
Next, once:  wrangler login
Then, from the checkout:
    tools/dist-linux.sh                       # build + verify both archives
    tools/dist-linux.sh --upload [--latest]   # publish to dl.popkern.dev
    tools/register-nightly-linux.sh           # nightly builds (tools/nightly.sh)
EOF
