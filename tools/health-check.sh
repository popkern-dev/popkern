#!/usr/bin/env bash
# Read-only project health/preflight checks for popkern.
set -euo pipefail

# The checks pipe their search tools through `|| true` (no match is not
# a failure), which also swallowed command-not-found: without the tool
# the gate scanned nothing and passed (PK-830). Require it up front.
for tool in rg jq; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$(basename "$0"): '$tool' is required but not on PATH" >&2
        exit 127
    fi
done

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

CONFIG=".claude/health-config.json"
pre_commit=false
metric=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pre-commit)
            pre_commit=true
            shift
            ;;
        --metric)
            metric="${2:-}"
            [[ -n "$metric" ]] || { echo "--metric requires a value" >&2; exit 2; }
            shift 2
            ;;
        -h|--help)
            cat <<'USAGE'
Usage: tools/health-check.sh [--pre-commit] [--metric name]

Metrics:
  file_lines
  blast_radius
  pub_api_changes
  atlas_staleness
  flowchart_staleness
  ip_clearance_gaps
  claude_md_budget
USAGE
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

cfg_num() {
    local key="$1" fallback="$2"
    jq -r --arg key "$key" '.thresholds[$key] // empty' "$CONFIG" 2>/dev/null | {
        read -r value
        if [[ -n "$value" && "$value" != "null" ]]; then
            echo "$value"
        else
            echo "$fallback"
        fi
    }
}

enabled() {
    local name="$1"
    if [[ -n "$metric" ]]; then
        [[ "$metric" == "$name" ]]
        return
    fi
    jq -e --arg name "$name" '.metrics[$name].enabled == true' "$CONFIG" >/dev/null 2>&1
}

base_ref() {
    if git rev-parse --verify origin/main >/dev/null 2>&1; then
        echo "origin/main...HEAD"
    elif git rev-parse --verify main >/dev/null 2>&1; then
        echo "main...HEAD"
    else
        echo "HEAD"
    fi
}

status_pass=0
status_warn=0
status_fail=0

pass() {
    echo "PASS: $*"
    status_pass=$((status_pass + 1))
}

warn() {
    echo "WARN: $*"
    status_warn=$((status_warn + 1))
}

fail() {
    echo "FAIL: $*"
    status_fail=$((status_fail + 1))
}

check_file_lines() {
    local warn_limit hard_limit violations hard
    warn_limit="$(cfg_num file_line_limit 450)"
    hard_limit="$(cfg_num file_line_hard_limit 500)"
    violations="$(find crates -name '*.rs' -exec wc -l {} + | awk -v warn="$warn_limit" '$2 != "total" && $1 >= warn {print $1 " " $2}' | sort -rn || true)"
    hard="$(printf '%s\n' "$violations" | awk -v hard="$hard_limit" '$1 > hard {print}' || true)"
    if [[ -n "$hard" ]]; then
        fail "file lines over hard limit $hard_limit"
        printf '%s\n' "$hard"
    elif [[ -n "$violations" ]]; then
        warn "file lines at or above warning limit $warn_limit"
        printf '%s\n' "$violations"
    else
        pass "file lines under $warn_limit"
    fi
}

check_blast_radius() {
    tools/blast-radius.sh
}

check_pub_api_changes() {
    local ref changes
    ref="$(base_ref)"
    if [[ "$ref" == "HEAD" ]]; then
        changes="$(git diff -- 'crates/*/src/**/*.rs' | rg '^\+.*pub (fn|struct|trait|enum|type) ' || true)"
    else
        changes="$(git diff "$ref" -- 'crates/*/src/**/*.rs' | rg '^\+.*pub (fn|struct|trait|enum|type) ' || true)"
    fi
    if [[ -n "$changes" ]]; then
        warn "public API additions/changes detected"
        printf '%s\n' "$changes"
    else
        pass "no public API additions detected"
    fi
}

check_atlas_staleness() {
    if [[ ! -f docs/atlas/index.md ]]; then
        warn "atlas missing docs/atlas/index.md"
        return
    fi
    local atlas_commit latest_pub_commit
    atlas_commit="$(git log -1 --format=%ct -- docs/atlas/index.md 2>/dev/null || echo 0)"
    latest_pub_commit="$(git log -1 --format=%ct -- crates 2>/dev/null || echo 0)"
    if [[ "${latest_pub_commit:-0}" -gt "${atlas_commit:-0}" ]]; then
        warn "atlas may be stale; run tools/generate-atlas.sh"
    else
        pass "atlas is not older than latest committed crate change"
    fi
}

check_flowchart_staleness() {
    shopt -s nullglob
    local dbs=(flowcharts/*.db)
    if [[ ${#dbs[@]} -eq 0 ]]; then
        pass "flowcharts: no built databases"
        return
    fi
    local stale=0
    for db in "${dbs[@]}"; do
        local subsystem seed_time refs
        subsystem="$(basename "$db" .db)"
        seed_time="$(git log -1 --format=%ct -- "flowcharts/seeds/$subsystem.sql" 2>/dev/null || echo 0)"
        refs="$(sqlite3 "$db" "SELECT DISTINCT substr(code_ref, 1, instr(code_ref, ':') - 1) FROM nodes WHERE code_ref IS NOT NULL AND instr(code_ref, ':') > 1;" 2>/dev/null || true)"
        while IFS= read -r ref; do
            [[ -n "$ref" && -f "$ref" ]] || continue
            local changed
            changed="$(git log -1 --format=%ct -- "$ref" 2>/dev/null || echo 0)"
            if [[ "${changed:-0}" -gt "${seed_time:-0}" ]]; then
                warn "flowchart $subsystem may be stale: $ref changed after seed"
                stale=1
            fi
        done <<< "$refs"
    done
    [[ "$stale" == "1" ]] || pass "flowcharts not stale by committed timestamps"
}

check_ip_clearance_gaps() {
    shopt -s nullglob
    local dbs=(flowcharts/*.db)
    if [[ ${#dbs[@]} -eq 0 ]]; then
        pass "IP clearance: no built flowchart databases"
        return
    fi
    local total=0
    for db in "${dbs[@]}"; do
        local count
        count="$(sqlite3 "$db" "SELECT COUNT(*) FROM annotations a JOIN nodes n ON a.node_id = n.id WHERE a.kind = 'comparison' AND (n.metadata IS NULL OR json_extract(n.metadata, '$.ip_status') IS NULL);" 2>/dev/null || echo 0)"
        if [[ "$count" -gt 0 ]]; then
            warn "IP clearance gaps in $(basename "$db"): $count comparison annotations"
            total=$((total + count))
        fi
    done
    [[ "$total" -gt 0 ]] || pass "IP clearance: no unchecked comparison annotations"
}

check_claude_md_budget() {
    local warn_limit hard_limit lines
    warn_limit="$(cfg_num claude_md_line_warn 145)"
    hard_limit="$(cfg_num claude_md_line_limit 160)"
    lines="$(wc -l < CLAUDE.md)"
    if [[ "$lines" -gt "$hard_limit" ]]; then
        fail "CLAUDE.md over budget: $lines/$hard_limit lines"
    elif [[ "$lines" -gt "$warn_limit" ]]; then
        warn "CLAUDE.md approaching budget: $lines/$hard_limit lines"
    else
        pass "CLAUDE.md within budget: $lines/$hard_limit lines"
    fi
}

run_metric() {
    local name="$1"
    enabled "$name" || return 0
    if [[ "$pre_commit" == true ]]; then
        case "$name" in
            file_lines|blast_radius|pub_api_changes|claude_md_budget) ;;
            *) return 0 ;;
        esac
    fi
    echo "== $name =="
    "check_$name"
    echo
}

echo "=== Project Health Check ==="
run_metric file_lines
run_metric blast_radius
run_metric pub_api_changes
run_metric atlas_staleness
run_metric flowchart_staleness
run_metric ip_clearance_gaps
run_metric claude_md_budget

echo "Summary: $status_pass pass, $status_warn warn, $status_fail fail"
[[ "$status_fail" -eq 0 ]]
