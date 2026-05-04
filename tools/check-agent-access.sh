#!/usr/bin/env bash
# Verify local tools and repo files needed by Claude/Codex workflows.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

missing=0

need_cmd() {
    local cmd="$1"
    if command -v "$cmd" >/dev/null 2>&1; then
        printf 'PASS command %-8s %s\n' "$cmd" "$(command -v "$cmd")"
    else
        printf 'FAIL command %-8s missing\n' "$cmd"
        missing=1
    fi
}

need_file() {
    local path="$1"
    if [[ -r "$path" ]]; then
        printf 'PASS file    %s\n' "$path"
    else
        printf 'FAIL file    %s not readable\n' "$path"
        missing=1
    fi
}

need_exec() {
    local path="$1"
    if [[ -x "$path" ]]; then
        printf 'PASS exec    %s\n' "$path"
    else
        printf 'FAIL exec    %s not executable\n' "$path"
        missing=1
    fi
}

for cmd in cargo rustc rustfmt jq sqlite3 cloc f3d rg codex; do
    need_cmd "$cmd"
done

need_file CLAUDE.md
need_file AGENTS.md
need_file .claude/health-config.json
need_file .claude/settings.local.json
need_file .codex/skills/popkern-workflow/SKILL.md

need_exec tools/task
need_exec tools/blast-radius.sh
need_exec tools/generate-atlas.sh
need_exec tools/check-semantic-names.sh
need_exec tools/health-check.sh
need_exec tools/quality-check.sh
need_exec .claude/hooks/log-prompt.sh
need_exec .claude/hooks/session-start.sh

echo
if [[ "$missing" -eq 0 ]]; then
    echo "Agent access check passed."
else
    echo "Agent access check failed."
fi

exit "$missing"
