#!/usr/bin/env bash
# What the nightly schedulers run, so a night runs tools/nightly.sh as
# origin/main has it, never as the registering checkout happens to have it
# (a branch, an old main, a rebase half done).
#
#   <nightly home>/bootstrap.sh <repository>
#
# tools/register-nightly-{windows.ps1,macos.sh,linux.sh} install a copy of
# this file as bootstrap.sh in the nightly home (the folder tools/nightly.sh
# keeps its logs and build worktree in) and schedule that copy, so nothing
# the scheduler starts lives in the checkout. Each run it
#
#   1. fetches origin/main into <repository>,
#   2. moves the runner worktree (<nightly home>/runner) to it, detached,
#   3. hands over to the runner's tools/nightly.sh.
#
# The runner is a worktree of <repository>, separate from the build
# worktree tools/nightly.sh moves to the night's agreed commit: bash reads a
# script as it runs it, so the folder holding the running script is never
# the one being checked out.
#
# This file changes rarely, and the installed copy does not follow main:
# when main's copy differs, the run says so in its log; re-run the
# register script to install it.
set -euo pipefail

[[ $# -ge 1 ]] || { echo "usage: bootstrap.sh <repository> [tools/nightly.sh arguments]" >&2; exit 2; }
REPO=$1; shift
NIGHTLY_HOME=$(cd "$(dirname "$0")" && pwd)
RUNNER=$NIGHTLY_HOME/runner
export POPKERN_NIGHTLY_HOME=$NIGHTLY_HOME

# tools/nightly.sh logs from its start on; what fails before that goes to
# the same night's log, so a failed bootstrap is not silent.
mkdir -p "$NIGHTLY_HOME/logs"
LOG=$NIGHTLY_HOME/logs/$(date -u +%Y%m%d).log
fail() { echo "== nightly bootstrap $(date -Iseconds): $*" >> "$LOG"; echo "$*" >&2; exit 1; }

[[ -d $REPO ]] || fail "repository $REPO is missing (re-run the register script from the checkout)"
git -C "$REPO" fetch --quiet origin main >> "$LOG" 2>&1 || fail "git fetch origin main failed in $REPO"

if ! git -C "$RUNNER" rev-parse --git-dir >/dev/null 2>&1; then
  # Missing, or left behind by a repository that moved: start it again.
  rm -rf "$RUNNER"
  git -C "$REPO" worktree prune
  git -C "$REPO" worktree add --quiet --detach "$RUNNER" origin/main >> "$LOG" 2>&1 ||
    fail "cannot add the runner worktree at $RUNNER"
else
  git -C "$RUNNER" checkout --quiet --detach --force origin/main >> "$LOG" 2>&1 ||
    fail "cannot move the runner worktree to origin/main"
fi

if ! cmp -s "$0" "$RUNNER/tools/nightly-bootstrap.sh"; then
  echo "== nightly bootstrap: origin/main's tools/nightly-bootstrap.sh differs from this installed copy; re-run the register script to install it" >> "$LOG"
fi
exec "$RUNNER/tools/nightly.sh" "$@"
