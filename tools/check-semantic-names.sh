#!/usr/bin/env bash
# Flag planning labels that should not become implementation names.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

scope="${1:---changed}"

case "$scope" in
    --changed)
        mapfile -t files < <(
            {
                git diff --name-only --diff-filter=ACMR
                git diff --cached --name-only --diff-filter=ACMR
                git ls-files --others --exclude-standard
            } | sort -u
        )
        ;;
    --all)
        mapfile -t files < <(rg --files crates docs | sort)
        ;;
    *)
        echo "Usage: tools/check-semantic-names.sh [--changed|--all]" >&2
        exit 2
        ;;
esac

check_files=()
for file in "${files[@]}"; do
    case "$file" in
        crates/*.rs|docs/*.md|docs/**/*.md)
            ;;
        *)
            continue
            ;;
    esac

    case "$file" in
        docs/superpowers/plans/*|docs/journal/*|docs/decisions/*)
            continue
            ;;
    esac

    [[ -f "$file" ]] && check_files+=("$file")
done

if [[ ${#check_files[@]} -eq 0 ]]; then
    echo "PASS: no changed long-lived code/docs to scan"
    exit 0
fi

patterns=(
    '\bPath[ _-]?[ABC]\b'
    '\bpath_[abc]\b'
    '\bPlan[ _-]?[ABC]\b'
    '\bplan_[abc]\b'
    '\bApproach[ _-]?[0-9]+\b'
    '\bapproach_[0-9]+\b'
    '\bPhase[0-9]+Impl\b'
    '\bphase_[0-9]+_impl\b'
    '\btemporary_builder\b'
    '\bnew_new_[A-Za-z0-9_]*\b'
)

regex="$(IFS='|'; echo "${patterns[*]}")"

matches="$(
    rg -n "$regex" "${check_files[@]}" || true
)"

if [[ -n "$matches" ]]; then
    echo "Planning labels found in long-lived code/docs. Use semantic domain names instead."
    echo "$matches"
    exit 1
fi

echo "PASS: no planning-label implementation names found"
