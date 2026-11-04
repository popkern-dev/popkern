#!/usr/bin/env bash
# Nightly build: build this machine's platform at the night's agreed
# commit with its tools/dist-<platform>.sh, and publish it to
# dl.popkern.dev as
#
#   <next patch>-nightly.<YYYYMMDD>.<commit>    e.g. 0.1.1-nightly.20261102.9852157f
#
# <next patch> is popkern-cli's version with the patch number raised:
# Cargo.toml carries the last release, and a nightly previews the next.
# <YYYYMMDD> is the UTC date, so build machines in different time zones
# name the same night alike.
#
#   tools/nightly.sh              # build + publish, unless already published
#   tools/nightly.sh --force      # rebuild and republish this night's version
#   tools/nightly.sh --no-upload  # build and verify only; bucket untouched
#
# Agreed commit: each platform builds on its own machine at its own time,
# and main can move in between. The first machine to run on a UTC date
# records origin/main in the bucket at nightly-commits/<YYYYMMDD>; every
# machine that night builds the recorded commit, so one release folder
# never mixes builds of different commits. The record is written, then read
# back and the read value adopted, so two machines that race on the write
# still agree unless both checked for the record before either wrote it:
# schedule the platforms apart (Windows runs at 02:00 local, macOS at
# 02:30, Linux at 03:00) rather than rely on that.
#
# The build runs in its own worktree of this repository
# ($POPKERN_NIGHTLY_HOME/checkout, default %LOCALAPPDATA%\popkern-nightly on
# Windows, ~/.local/state/popkern-nightly elsewhere), detached at the agreed
# commit, so it never touches the branch or working tree you are editing,
# and it keeps its own target/ between nights.
#
# A nightly never moves `latest` (what a plain install gets). It moves the
# `nightly` pointer, so others install the newest nightly with the
# installers as they are:
#
#   $env:POPKERN_VERSION = (irm https://dl.popkern.dev/nightly).Trim()
#   irm https://popkern.dev/install.ps1 | iex
#
#   curl -fsSL https://popkern.dev/install.sh |
#     POPKERN_VERSION=$(curl -fsSL https://dl.popkern.dev/nightly) sh
#
# Open item for the second platform: `nightly` moves when any platform
# publishes, so until the others publish the same version, an install on
# those platforms finds no build for itself.
#
# Logs: $POPKERN_NIGHTLY_HOME/logs/<date>.log (kept 30 days). Publishing
# needs wrangler logged in for this user (`wrangler login`), or
# CLOUDFLARE_API_TOKEN + CLOUDFLARE_ACCOUNT_ID in the user environment.
# Scheduling: tools/register-nightly-windows.ps1 on Windows,
# tools/register-nightly-macos.sh on macOS, tools/register-nightly-linux.sh
# on Linux.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FORCE=0; UPLOAD=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE=1 ;;
    --no-upload) UPLOAD=0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done

# The platform's dist script, and the archive-name suffix its builds share
# (to tell whether this platform already published a version).
case "$(uname -s)" in
  MINGW* | MSYS* | CYGWIN*) PLATFORM=windows; ARCHIVES='-pc-windows-msvc\.zip$' ;;
  Linux) PLATFORM=linux; ARCHIVES='-unknown-linux-musl\.tar\.gz$' ;;
  Darwin) PLATFORM=macos; ARCHIVES='-apple-darwin\.tar\.gz$' ;;
  *) echo "no nightly build for $(uname -s)" >&2; exit 1 ;;
esac

# A scheduled job starts with a bare environment: put rustup's and npm's
# shims (cargo, wasm-bindgen, wrangler) and node (which wrangler runs on)
# on PATH explicitly.
if [[ $PLATFORM == windows ]]; then
  export PATH="$HOME/.cargo/bin:$(cygpath -u "$APPDATA")/npm:$(cygpath -u "$PROGRAMFILES")/nodejs:$PATH"
  DEFAULT_HOME=$(cygpath -u "$LOCALAPPDATA")/popkern-nightly
else
  # ~/.local/bin: node, wrangler and zig as tools/setup-macos-build.sh and
  # tools/setup-linux-build.sh install them.
  export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
  DEFAULT_HOME=${XDG_STATE_HOME:-$HOME/.local/state}/popkern-nightly
fi
NIGHTLY_HOME=${POPKERN_NIGHTLY_HOME:-$DEFAULT_HOME}
CHECKOUT=$NIGHTLY_HOME/checkout
mkdir -p "$NIGHTLY_HOME/logs"
NIGHT=$(date -u +%Y%m%d)
exec > >(tee -a "$NIGHTLY_HOME/logs/$NIGHT.log") 2>&1
find "$NIGHTLY_HOME/logs" -name '*.log' -mtime +30 -delete
echo "== nightly $PLATFORM $(date -Iseconds) (night $NIGHT UTC)"

# wrangler from the checkout, never from wherever the scheduler starts:
# run in $HOME, it creates ~/.wrangler (its project cache), which it then
# takes for its legacy global config folder, and every wrangler command
# for this user, scheduled or not, stops finding the login.
cd "$REPO"
source "$REPO/tools/dist-publish.sh"
SCRATCH=$(mktemp -d)
trap 'rm -rf "$SCRATCH"' EXIT

# bucket_get <key> <file>: 0 when fetched, 1 when the key does not exist;
# any other failure stops the run rather than be read as "missing".
bucket_get() {
  if wrangler r2 object get "$DIST_BUCKET/$1" --remote --pipe > "$2" 2> "$SCRATCH/err"; then
    return 0
  fi
  grep -q 'specified key does not exist' "$SCRATCH/err" && return 1
  echo "cannot read $DIST_BUCKET/$1:" >&2; cat "$SCRATCH/err" >&2
  exit 1
}

git -C "$REPO" fetch --quiet origin main
MAIN=$(git -C "$REPO" rev-parse origin/main)
RECORD=nightly-commits/$NIGHT
if [[ $UPLOAD == 0 ]]; then
  # A dry run reads the night's record but never writes one.
  if bucket_get "$RECORD" "$SCRATCH/commit"; then SHA=$(tr -d '[:space:]' < "$SCRATCH/commit"); else SHA=$MAIN; fi
else
  if ! bucket_get "$RECORD" "$SCRATCH/commit"; then
    printf '%s\n' "$MAIN" > "$SCRATCH/commit"
    dist_put "$RECORD" "$SCRATCH/commit" "text/plain; charset=utf-8" "$DIST_SHORT"
    bucket_get "$RECORD" "$SCRATCH/commit" || { echo "$RECORD vanished after writing it" >&2; exit 1; }
  fi
  SHA=$(tr -d '[:space:]' < "$SCRATCH/commit")
fi
[[ $SHA =~ ^[0-9a-f]{40}$ ]] || { echo "$RECORD holds '$SHA', not a commit" >&2; exit 1; }
if [[ $SHA != "$MAIN" ]]; then
  echo "night $NIGHT is agreed at ${SHA:0:8}; origin/main has since moved to ${MAIN:0:8}"
  git -C "$REPO" cat-file -e "$SHA^{commit}" 2>/dev/null || git -C "$REPO" fetch --quiet origin "$SHA"
fi

BASE=$(git -C "$REPO" show "$SHA:crates/popkern-cli/Cargo.toml" | awk -F'"' '/^version = /{print $2; exit}')
IFS=. read -r MAJOR MINOR PATCH <<< "${BASE%%-*}"
VERSION="$MAJOR.$MINOR.$((PATCH + 1))-nightly.$NIGHT.${SHA:0:8}"

if [[ $UPLOAD == 1 && $FORCE == 0 ]] && bucket_get "v$VERSION/SHA256SUMS" "$SCRATCH/sums" &&
   grep -qE "  popkern-$VERSION-.*$ARCHIVES" "$SCRATCH/sums"; then
  echo "$PLATFORM build of $VERSION is already published: nothing to do"
  exit 0
fi

if [[ ! -d $CHECKOUT ]]; then
  git -C "$REPO" worktree prune
  git -C "$REPO" worktree add --detach "$CHECKOUT" "$SHA"
else
  git -C "$CHECKOUT" checkout --quiet --detach --force "$SHA"
  # Everything but the build cache: a stale web/live/pkg or leftover
  # dist folder must not leak into tonight's archives.
  git -C "$CHECKOUT" clean -ffdxq -e /target
fi
cd "$CHECKOUT"
echo "building $VERSION: $(git log -1 --format=%s)"

if [[ $UPLOAD == 0 ]]; then
  "tools/dist-$PLATFORM.sh" --version "$VERSION"
  exit 0
fi
"tools/dist-$PLATFORM.sh" --version "$VERSION" --upload

# Point `nightly` at this build only once its archives and the merged
# SHA256SUMS are up, so the pointer never names a half-published release.
printf '%s\n' "$VERSION" > "$SCRATCH/nightly"
dist_put nightly "$SCRATCH/nightly" "text/plain; charset=utf-8" "public, max-age=60"
echo "== published popkern $VERSION ($PLATFORM)"
