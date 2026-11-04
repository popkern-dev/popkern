#!/usr/bin/env bash
# Register the cron job that runs tools/nightly.sh at 03:00 every day, as
# the current user. The Linux counterpart of
# tools/register-nightly-windows.ps1 and tools/register-nightly-macos.sh.
# Safe to re-run: it replaces the job.
#
#   tools/register-nightly-linux.sh                 # daily 03:00 local time
#   tools/register-nightly-linux.sh --at 04:15      # another local start time
#   tools/register-nightly-linux.sh --unregister
#
# 03:00, after Windows (02:00) and macOS (02:30): tools/nightly.sh has
# build machines start apart, because the first one to run on a night
# records the commit every platform builds that night.
#
# cron, not a systemd timer: it is on every Linux, with or without
# systemd as init, and a user crontab runs whether or not the user is
# logged in, without root. The job is two lines of your crontab tagged
# "# popkern-nightly"; other lines are left as they are. The second line
# (@reboot) catches up a start missed while the machine was off, which
# cron itself never does; it is a no-op once tonight's Linux build is
# published (tools/nightly.sh checks the bucket). flock keeps a second
# copy from starting while one runs, and a run is stopped after 6 hours,
# so a hang cannot hold the next night up. The job waits up to 5 minutes
# for a default route, which a boot brings up after cron starts.
#
# Not covered, unlike macOS and Windows: waking a sleeping machine (that
# needs root, e.g. rtcwake) and holding sleep off during a build (logind
# refuses a cron job's inhibitor lock without interactive authentication).
# Build on a machine that stays awake: a night slept through is skipped.
#
# Prerequisites: tools/setup-linux-build.sh, then `wrangler login` once.
set -euo pipefail

TAG='# popkern-nightly'
AT=03:00
LIMIT_HOURS=6
UNREGISTER=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --at) AT=${2:-}; shift ;;
    --unregister) UNREGISTER=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
[[ $(uname -s) == Linux ]] || { echo "register-nightly-linux.sh is for Linux" >&2; exit 1; }
[[ $AT =~ ^([01][0-9]|2[0-3]):([0-5][0-9])$ ]] || { echo "--at takes a 24-hour HH:MM, not '$AT'" >&2; exit 2; }
HOUR=$((10#${BASH_REMATCH[1]})); MINUTE=$((10#${BASH_REMATCH[2]}))
for tool in crontab flock timeout; do
  command -v "$tool" >/dev/null || { echo "missing $tool (install cron / util-linux / coreutils)" >&2; exit 1; }
done

SCRIPT=$(cd "$(dirname "$0")" && pwd)/nightly.sh
[[ -x $SCRIPT ]] || { echo "missing $SCRIPT" >&2; exit 1; }
# Where tools/nightly.sh keeps its worktree and logs (its own default).
NIGHTLY_HOME=${POPKERN_NIGHTLY_HOME:-${XDG_STATE_HOME:-$HOME/.local/state}/popkern-nightly}

# The crontab without our lines; an empty or missing crontab is no error.
others=$(crontab -l 2>/dev/null | grep -vF "$TAG" || true)

if [[ $UNREGISTER == 1 ]]; then
  if [[ -n $others ]]; then printf '%s\n' "$others" | crontab -; else crontab -r 2>/dev/null || true; fi
  echo "removed the popkern-nightly cron job"
  exit 0
fi

# cron runs the line with /bin/sh, PATH=/usr/bin:/bin and HOME set;
# tools/nightly.sh puts the build tools on PATH itself. Arguments are
# single-quoted for sh, and cron reads a bare % as a newline.
q() { printf "'%s'" "${1//\'/\'\\\'\'}"; }
# The default route from /proc, not ip(8), which cron's PATH may not reach.
WAIT_FOR_NETWORK='for _ in $(seq 60); do awk '"'"'$2 == "00000000" { up = 1 } END { exit !up }'"'"' /proc/net/route && break; sleep 5; done; exec "$0"'
COMMAND="POPKERN_NIGHTLY_HOME=$(q "$NIGHTLY_HOME") flock -n $(q "$NIGHTLY_HOME/lock")"
COMMAND+=" timeout --kill-after=5m ${LIMIT_HOURS}h /bin/bash -c $(q "$WAIT_FOR_NETWORK") $(q "$SCRIPT")"
# tools/nightly.sh tees its own log into logs/; cron.log catches only
# what fails before that starts.
COMMAND+=" >/dev/null 2>>$(q "$NIGHTLY_HOME/logs/cron.log")"
[[ $COMMAND != *%* ]] || { echo "a path contains %, which cron cannot run" >&2; exit 1; }

mkdir -p "$NIGHTLY_HOME/logs"
{
  [[ -n $others ]] && printf '%s\n' "$others"
  printf '%d %d * * * %s %s\n' "$MINUTE" "$HOUR" "$COMMAND" "$TAG"
  printf '@reboot %s %s\n' "$COMMAND" "$TAG"
} | crontab -
pgrep -x cron >/dev/null || pgrep -x crond >/dev/null ||
  echo "warning: no cron daemon is running, so the job will not start" >&2
echo "registered popkern-nightly: daily $AT ($(date +%Z)) and at boot, as $(id -un), logged in or not"
echo "run now:  $SCRIPT"
echo "logs:     $NIGHTLY_HOME/logs"
