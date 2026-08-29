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
        docs/superpowers/*|docs/journal/*|docs/decisions/*|docs/reports/*|docs/atlas/*)
            # journals, ADRs, generated reports/atlas: historical or
            # regenerated content.
            continue
            ;;
        crates/popkern-fillet/src/unfillet_plan_contract_tests.rs)
            # Contains the banned literals as NEGATIVE assertions (error
            # messages must not leak plan coordinates).
            continue
            ;;
        docs/architecture/kp20-*|docs/architecture/shared-corner-repeated-apex-substrate.md|docs/reference/kp20-bi6-deep-interior.md|docs/architecture/internal-innovations.md)
            # Historical campaign knowledge records; filename renames are
            # tracked by the KP semantic-rename task (plan section 5.2
            # item 7).
            continue
            ;;
        docs/architecture/public-release-readiness-plan.md|docs/reference/kp-naming-map.md)
            # These docs quote the banned labels as their own
            # rename/disposition inventory.
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
    '\bpath_[abc]_[a-z0-9_]+\b'
    '\bPlan[ _-]?[ABC]\b'
    '\bplan_[abc]\b'
    '\bApproach[ _-]?[0-9]+\b'
    '\bapproach_[0-9]+\b'
    '\bPhase[0-9]+Impl\b'
    '\bphase_[0-9]+_impl\b'
    '\bphase_[0-9]+[a-z]?(_[a-z0-9]+)+\b'
    '\bPhase [0-9]+[a-z]?-[A-Za-z][A-Za-z0-9.]*\b'
    '\bPhase 6c\b'
    '\bPhase [DEGHI]\b'
    '\bPhase I[IV]*\b'
    'Round [0-9]+[a-z]?'
    '\btemporary_builder\b'
    '\bnew_new_[A-Za-z0-9_]*\b'
    '_tmp\.rs\b'
)

# New identifiers must not carry ticket or paper-codename prefixes as
# bare test-file names (pkNNN_*.rs) — semantic names with the ticket
# cited in a //! line instead.
file_patterns=(
    'pk[0-9]{3}_[a-z0-9_]+\.rs$'
    'phase_[0-9]+[a-z]?[a-z0-9_]*\.rs$'
)

regex="$(IFS='|'; echo "${patterns[*]}")"

matches="$(
    rg -n "$regex" "${check_files[@]}" || true
)"

# Bare "legacy" is banned in code (name the superseded mechanism);
# docs keep their historical vocabulary until the export-step pruning.
crate_files=()
for file in "${check_files[@]}"; do
    [[ "$file" == crates/* ]] && crate_files+=("$file")
done
if [[ ${#crate_files[@]} -gt 0 ]]; then
    legacy_matches="$(rg -n '[Ll]egacy|LEGACY' "${crate_files[@]}" || true)"
else
    legacy_matches=""
fi
matches="${matches}${legacy_matches:+
$legacy_matches}"

file_regex="$(IFS='|'; echo "${file_patterns[*]}")"
bad_names="$(
    printf '%s\n' "${check_files[@]}" | rg "$file_regex" || true
)"

fail=0
if [[ -n "$matches" ]]; then
    echo "Planning labels found in long-lived code/docs. Use semantic domain names instead."
    echo "$matches"
    fail=1
fi
if [[ -n "$bad_names" ]]; then
    echo "Ticket/plan-coordinate file names found. Rename semantically (cite the ticket in a //! line)."
    echo "$bad_names"
    fail=1
fi
[[ $fail -eq 1 ]] && exit 1

echo "PASS: no planning-label implementation names found"
