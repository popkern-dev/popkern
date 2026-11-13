#!/usr/bin/env bash
# Register the launchd job that runs tools/nightly.sh at 02:30 every day,
# as the current user. The macOS counterpart of
# tools/register-nightly-windows.ps1. Safe to re-run: it replaces the job.
#
# The job runs an installed copy of tools/nightly-bootstrap.sh
# (<nightly home>/bootstrap.sh), which runs origin/main's tools/nightly.sh
# from a runner worktree: a night never depends on what this checkout has
# checked out. Re-run this script when the bootstrap changes on main (a
# night's log says so); the checkout must stay where it is, since the
# runner and build worktrees share its repository.
#
#   tools/register-nightly-macos.sh                 # runs while you are logged in
#   sudo tools/register-nightly-macos.sh            # runs whether or not you are logged in
#   tools/register-nightly-macos.sh --at 03:15      # another local start time
#   [sudo] tools/register-nightly-macos.sh --unregister
#
# 02:30 rather than the Windows build's 02:00: tools/nightly.sh has build
# machines start apart, because the first one to run on a night records
# the commit every platform builds that night.
#
# The job survives restarts (launchd reloads it from its plist). A start
# missed while the Mac slept runs when it wakes; one missed while it was
# off runs at the next boot (with sudo) or login (without), and the job
# also runs once right after registering. Each of those is a no-op once
# tonight's macOS build is published (tools/nightly.sh checks the bucket).
# caffeinate keeps the Mac awake while it builds. A run is stopped after
# 6 hours, and launchd then ends its process group: launchd never starts a
# second copy of a running job, so a hang must not hold the next night up.
#
# With sudo it is a LaunchDaemon that runs as you (never as root), logged
# in or not, and a `pmset repeat` event wakes a sleeping Mac 2 minutes
# before the start. pmset keeps one repeating event per machine: a
# different one already set is left alone, and reported. Without sudo it
# is a LaunchAgent, which runs only while you are logged in, and says so.
#
# Prerequisites: tools/setup-macos-build.sh, then `wrangler login` once.
set -euo pipefail

LABEL=dev.popkern.nightly-macos
AT=02:30
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
[[ $(uname -s) == Darwin ]] || { echo "register-nightly-macos.sh is for macOS" >&2; exit 1; }
[[ $AT =~ ^([01][0-9]|2[0-3]):([0-5][0-9])$ ]] || { echo "--at takes a 24-hour HH:MM, not '$AT'" >&2; exit 2; }
HOUR=$((10#${BASH_REMATCH[1]})); MINUTE=$((10#${BASH_REMATCH[2]}))

TOOLS=$(cd "$(dirname "$0")" && pwd)
REPO=$(dirname "$TOOLS")
[[ -f $TOOLS/nightly-bootstrap.sh ]] || { echo "missing $TOOLS/nightly-bootstrap.sh" >&2; exit 1; }

# Who the job runs as: you, also under sudo.
if [[ $EUID == 0 ]]; then
  [[ -n ${SUDO_USER:-} && $SUDO_USER != root ]] ||
    { echo "run with sudo from your own account, not as root: the job runs as that account" >&2; exit 1; }
  RUN_USER=$SUDO_USER; MODE=daemon
else
  RUN_USER=$(id -un); MODE=agent
fi
RUN_UID=$(id -u "$RUN_USER")
RUN_HOME=$(dscl . -read "/Users/$RUN_USER" NFSHomeDirectory | awk '{print $2}')
AGENT_PLIST=$RUN_HOME/Library/LaunchAgents/$LABEL.plist
DAEMON_PLIST=/Library/LaunchDaemons/$LABEL.plist
# Where tools/nightly.sh keeps its worktrees and logs (its own default).
NIGHTLY_HOME=${POPKERN_NIGHTLY_HOME:-$RUN_HOME/.local/state/popkern-nightly}
BOOTSTRAP=$NIGHTLY_HOME/bootstrap.sh

# The wake event: 2 minutes before the start, as pmset sets and prints it.
WAKE=$(( (HOUR * 60 + MINUTE + 1440 - 2) % 1440 ))
WAKE_SET=$(printf '%02d:%02d:00' $((WAKE / 60)) $((WAKE % 60)))
WAKE_SHOWN="wakepoweron at $(date -j -f %H:%M:%S "$WAKE_SET" +%-I:%M%p) every day"
our_wake_is_set() { pmset -g sched | grep -qF "$WAKE_SHOWN"; }

remove_agent() {
  launchctl bootout "gui/$RUN_UID/$LABEL" 2>/dev/null || true
  rm -f "$AGENT_PLIST"
}

if [[ $UNREGISTER == 1 ]]; then
  if [[ $MODE == daemon ]]; then
    launchctl bootout "system/$LABEL" 2>/dev/null || true
    rm -f "$DAEMON_PLIST"
    if our_wake_is_set; then pmset repeat cancel; echo "cancelled the pmset wake ($WAKE_SHOWN)"; fi
  elif [[ -e $DAEMON_PLIST ]]; then
    echo "the job is a LaunchDaemon: remove it with  sudo $0 --unregister" >&2; exit 1
  fi
  remove_agent
  echo "removed $LABEL"
  exit 0
fi

if [[ $MODE == agent && -e $DAEMON_PLIST ]]; then
  echo "$LABEL is registered as a LaunchDaemon: re-run with sudo to replace it" >&2; exit 1
fi

# Text for the plist's <string> elements.
xml() { printf '%s' "$1" | sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g'; }

# perl's alarm outlives its exec, so the time limit covers the whole run;
# the route check waits up to 5 minutes for the network a boot or wake
# brings up only after launchd starts the job.
WAIT_FOR_NETWORK='for _ in $(seq 60); do /sbin/route -n get default >/dev/null 2>&1 && break; sleep 5; done; exec "$0" "$@"'
ARGS=(/usr/bin/perl -e 'alarm shift; exec @ARGV or die "exec: $!\n"' $((LIMIT_HOURS * 3600))
      /usr/bin/caffeinate -i /bin/bash -c "$WAIT_FOR_NETWORK" "$BOOTSTRAP" "$REPO")

plist() {
  local arg user_keys=''
  [[ $MODE == daemon ]] && user_keys="  <key>UserName</key><string>$(xml "$RUN_USER")</string>"
  cat <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>$LABEL</string>
$user_keys
  <key>ProgramArguments</key>
  <array>
$(for arg in "${ARGS[@]}"; do printf '    <string>%s</string>\n' "$(xml "$arg")"; done)
  </array>
  <key>StartCalendarInterval</key>
  <dict><key>Hour</key><integer>$HOUR</integer><key>Minute</key><integer>$MINUTE</integer></dict>
  <key>RunAtLoad</key><true/>
  <key>EnvironmentVariables</key>
  <dict>
    <key>HOME</key><string>$(xml "$RUN_HOME")</string>
    <key>USER</key><string>$(xml "$RUN_USER")</string>
    <key>PATH</key><string>/usr/bin:/bin:/usr/sbin:/sbin</string>
    <key>POPKERN_NIGHTLY_HOME</key><string>$(xml "$NIGHTLY_HOME")</string>
  </dict>
  <key>WorkingDirectory</key><string>$(xml "$REPO")</string>
  <key>StandardOutPath</key><string>/dev/null</string>
  <key>StandardErrorPath</key><string>$(xml "$NIGHTLY_HOME/logs/launchd.log")</string>
</dict>
</plist>
EOF
}
# tools/nightly.sh tees its own log into logs/; launchd.log catches only
# what fails before that starts.

if [[ $MODE == daemon ]]; then
  sudo -u "$RUN_USER" mkdir -p "$NIGHTLY_HOME/logs"
  sudo -u "$RUN_USER" install -m 755 "$TOOLS/nightly-bootstrap.sh" "$BOOTSTRAP"
  remove_agent
  tmp=$(mktemp)
  plist > "$tmp"
  plutil -lint -s "$tmp"
  install -m 644 -o root -g wheel "$tmp" "$DAEMON_PLIST"
  rm -f "$tmp"
  launchctl bootout "system/$LABEL" 2>/dev/null || true
  launchctl bootstrap system "$DAEMON_PLIST"
  echo "registered $LABEL: daily $AT as $RUN_USER, whether or not logged in"
  if our_wake_is_set; then
    echo "wake event already set: $WAKE_SHOWN"
  elif pmset -g sched | grep -q 'Repeating power events'; then
    echo "NOT setting the wake event: pmset already has a repeating event (pmset -g sched)."
    echo "Replace it with  sudo pmset repeat wakeorpoweron MTWRFSU $WAKE_SET"
  else
    pmset repeat wakeorpoweron MTWRFSU "$WAKE_SET"
    echo "wakes the Mac at $WAKE_SET (pmset repeat)"
  fi
  echo "run now:  sudo launchctl kickstart system/$LABEL"
else
  mkdir -p "$NIGHTLY_HOME/logs" "$(dirname "$AGENT_PLIST")"
  install -m 755 "$TOOLS/nightly-bootstrap.sh" "$BOOTSTRAP"
  plist > "$AGENT_PLIST"
  plutil -lint -s "$AGENT_PLIST"
  launchctl bootout "gui/$RUN_UID/$LABEL" 2>/dev/null || true
  launchctl bootstrap "gui/$RUN_UID" "$AGENT_PLIST"
  echo "registered $LABEL: daily $AT as $RUN_USER, ONLY while logged in"
  echo "(re-run with sudo to run while logged out and wake the Mac from sleep)"
  echo "run now:  launchctl kickstart gui/$RUN_UID/$LABEL"
fi
echo "logs:     $NIGHTLY_HOME/logs"
