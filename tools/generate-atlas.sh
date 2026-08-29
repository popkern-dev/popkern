#!/usr/bin/env bash
set -euo pipefail

# Move to repo root (parent of the tools/ directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

ATLAS_DIR="docs/atlas"

# Clean and recreate atlas directory
rm -rf "$ATLAS_DIR"
mkdir -p "$ATLAS_DIR"

# ──────────────────────────────────────────────────────────────────────────────
# Helper: get description for a crate
# Try Cargo.toml description first, then //! doc comment in lib.rs
# ──────────────────────────────────────────────────────────────────────────────
get_description() {
    local crate_dir="$1"
    local cargo_toml="$crate_dir/Cargo.toml"
    local lib_rs="$crate_dir/src/lib.rs"

    # Try Cargo.toml description field
    local desc
    desc=$({ grep '^description' "$cargo_toml" 2>/dev/null || true; } \
           | sed 's/^description\s*=\s*"//' \
           | sed 's/"$//' \
           | head -1)
    if [[ -n "$desc" ]]; then
        echo "$desc"
        return
    fi

    # Fall back to first //! line in lib.rs
    if [[ -f "$lib_rs" ]]; then
        local crate_name
        crate_name=$(basename "$crate_dir")
        desc=$({ grep '^//!' "$lib_rs" || true; } \
               | sed 's|^//!\s*||' \
               | { grep -v '^$' || true; } \
               | head -1 \
               | sed "s|^${crate_name}:[[:space:]]*||")
        if [[ -n "$desc" ]]; then
            echo "$desc"
            return
        fi
    fi

    echo "(no description)"
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: count source files (excluding test files and tests/ dirs)
# ──────────────────────────────────────────────────────────────────────────────
count_source_files() {
    local src_dir="$1"
    find "$src_dir" -name "*.rs" \
        ! -name "*_tests.rs" \
        ! -name "test_*.rs" \
        ! -path "*/tests/*" \
        2>/dev/null \
        | wc -l | tr -d ' '
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: count total lines (excluding test files)
# ──────────────────────────────────────────────────────────────────────────────
count_lines() {
    local src_dir="$1"
    find "$src_dir" -name "*.rs" \
        ! -name "*_tests.rs" \
        ! -name "test_*.rs" \
        ! -path "*/tests/*" \
        2>/dev/null \
        -exec cat {} + \
        | wc -l | tr -d ' '
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: count #[test] annotations (in ALL files, including test files)
# ──────────────────────────────────────────────────────────────────────────────
count_tests() {
    local src_dir="$1"
    { grep -rh '#\[test\]' "$src_dir" 2>/dev/null || true; } | wc -l | tr -d ' '
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: get internal (popkern-*) dependencies from Cargo.toml
# ──────────────────────────────────────────────────────────────────────────────
get_internal_deps() {
    local cargo_toml="$1"
    { grep '^popkern-' "$cargo_toml" 2>/dev/null || true; } \
        | sed 's/\s*=.*//' \
        | sort \
        | tr '\n' ',' \
        | sed 's/,$//' \
        | sed 's/,/, /g'
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: get reverse dependencies (which crates depend on this one)
# ──────────────────────────────────────────────────────────────────────────────
get_reverse_deps() {
    local crate_name="$1"
    local result=""
    for other_dir in crates/popkern-*/; do
        local other_crate
        other_crate=$(basename "$other_dir")
        if [[ "$other_crate" == "$crate_name" ]]; then
            continue
        fi
        local other_toml="$other_dir/Cargo.toml"
        if grep -q "^${crate_name}\s*=" "$other_toml" 2>/dev/null; then
            result="${result}${other_crate}, "
        fi
    done
    # Remove trailing ", "
    echo "${result%, }"
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: extract public API items from source files (not test files)
# Returns lines in format: KIND|NAME|DOC
# KIND is one of: type, trait, fn
# ──────────────────────────────────────────────────────────────────────────────
extract_public_api() {
    local src_dir="$1"

    # Find all non-test source files
    local files
    mapfile -t files < <(find "$src_dir" -name "*.rs" \
        ! -name "*_tests.rs" \
        ! -name "test_*.rs" \
        ! -path "*/tests/*" \
        2>/dev/null \
        | sort)

    for file in "${files[@]}"; do
        local prev_doc=""
        local in_doc=0
        while IFS= read -r line; do
            # Accumulate leading /// lines; reset on non-/// non-blank
            if [[ "$line" =~ ^[[:space:]]*///' ' ]]; then
                # First /// line: capture as summary; subsequent lines ignored
                if [[ $in_doc -eq 0 ]]; then
                    prev_doc="${line#*/// }"
                    prev_doc="${prev_doc## }"
                    in_doc=1
                fi
            elif [[ "$line" =~ ^[[:space:]]*'///' ]]; then
                # Empty /// line
                : # keep current doc state
            else
                # Check for pub item declarations
                # Match: pub struct, pub enum, pub type (-> "type" kind)
                if [[ "$line" =~ ^pub[[:space:]]+(struct|enum|type)[[:space:]]+([A-Za-z_][A-Za-z0-9_]*) ]]; then
                    local name="${BASH_REMATCH[2]}"
                    echo "type|${name}|${prev_doc}"
                elif [[ "$line" =~ ^pub[[:space:]]+trait[[:space:]]+([A-Za-z_][A-Za-z0-9_]*) ]]; then
                    local name="${BASH_REMATCH[1]}"
                    echo "trait|${name}|${prev_doc}"
                elif [[ "$line" =~ ^pub[[:space:]]+fn[[:space:]]+([A-Za-z_][A-Za-z0-9_]*) ]]; then
                    local name="${BASH_REMATCH[1]}"
                    echo "fn|${name}|${prev_doc}"
                elif [[ "$line" =~ ^pub[[:space:]]+(async[[:space:]]+fn|unsafe[[:space:]]+fn)[[:space:]]+([A-Za-z_][A-Za-z0-9_]*) ]]; then
                    local name="${BASH_REMATCH[2]}"
                    echo "fn|${name}|${prev_doc}"
                else
                    # Non-pub, non-doc line: reset doc accumulation
                    # But don't reset on attributes (#[...]) or blank lines
                    if [[ ! "$line" =~ ^[[:space:]]*'#' ]] && \
                       [[ ! "$line" =~ ^[[:space:]]*$ ]]; then
                        prev_doc=""
                        in_doc=0
                    fi
                    continue
                fi
                # After capturing an item, reset doc state
                prev_doc=""
                in_doc=0
                continue
            fi
            # If we got here we were on a doc line — don't reset
        done < "$file"
    done
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper: generate module structure listing (flat, relative to crate root)
# Excludes test files
# ──────────────────────────────────────────────────────────────────────────────
get_module_structure() {
    local src_dir="$1"
    # Strip the crate root prefix (everything before "src/") from paths
    local crate_root
    crate_root="$(dirname "$src_dir")/"
    find "$src_dir" -name "*.rs" \
        ! -name "*_tests.rs" \
        ! -name "test_*.rs" \
        ! -path "*/tests/*" \
        2>/dev/null \
        | sed "s|${crate_root}||" \
        | sort
}

# ──────────────────────────────────────────────────────────────────────────────
# Collect all crate data first (needed for index table)
# ──────────────────────────────────────────────────────────────────────────────
declare -a ALL_CRATES=()
declare -A CRATE_DESC=()
declare -A CRATE_FILES=()
declare -A CRATE_LINES=()
declare -A CRATE_TESTS=()
declare -A CRATE_DEPS=()

for crate_dir in crates/popkern-*/; do
    crate=$(basename "$crate_dir")
    src_dir="${crate_dir}src"

    # Skip crates with no src directory
    if [[ ! -d "$src_dir" ]]; then
        continue
    fi

    ALL_CRATES+=("$crate")
    CRATE_DESC["$crate"]=$(get_description "$crate_dir")
    CRATE_FILES["$crate"]=$(count_source_files "$src_dir")
    CRATE_LINES["$crate"]=$(count_lines "$src_dir")
    CRATE_TESTS["$crate"]=$(count_tests "$src_dir")
    CRATE_DEPS["$crate"]=$(get_internal_deps "${crate_dir}Cargo.toml")
done

# ──────────────────────────────────────────────────────────────────────────────
# Generate per-crate pages
# ──────────────────────────────────────────────────────────────────────────────
for crate in "${ALL_CRATES[@]}"; do
    echo "Generating atlas for ${crate}..."

    crate_dir="crates/${crate}"
    src_dir="${crate_dir}/src"
    out_file="${ATLAS_DIR}/${crate}.md"

    desc="${CRATE_DESC[$crate]}"
    files="${CRATE_FILES[$crate]}"
    lines="${CRATE_LINES[$crate]}"
    tests="${CRATE_TESTS[$crate]}"
    deps="${CRATE_DEPS[$crate]}"
    rev_deps=$(get_reverse_deps "$crate")

    # Format deps display
    deps_display="${deps:-(none)}"
    rev_deps_display="${rev_deps:-(none)}"

    # Extract public API
    mapfile -t api_items < <(extract_public_api "$src_dir" | sort -u)

    # Separate into types, traits, functions
    declare -a types=()
    declare -a traits=()
    declare -a fns=()
    for item in "${api_items[@]}"; do
        IFS='|' read -r kind name doc <<< "$item"
        case "$kind" in
            type)  types+=("$name|$doc") ;;
            trait) traits+=("$name|$doc") ;;
            fn)    fns+=("$name|$doc") ;;
        esac
    done

    # Module structure
    mapfile -t modules < <(get_module_structure "$src_dir")

    # Determine entry point
    if [[ -f "${src_dir}/lib.rs" ]]; then
        entry_point="src/lib.rs"
    elif [[ -f "${src_dir}/main.rs" ]]; then
        entry_point="src/main.rs"
    else
        entry_point="(none)"
    fi

    # ── Write markdown page ────────────────────────────────────────────────────
    {
        echo "# ${crate}"
        echo ""
        echo "> ${desc}"
        echo ""
        echo "## Metrics"
        echo ""
        echo "| Metric | Value |"
        echo "|--------|-------|"
        echo "| Source files | ${files} |"
        echo "| Lines of code | ${lines} |"
        echo "| Tests | ${tests} |"
        echo ""
        echo "## Dependencies"
        echo ""
        echo "**Uses:** ${deps_display}"
        echo ""
        echo "**Used by:** ${rev_deps_display}"
        echo ""
        echo "## Public API"
        echo ""

        # Types section
        echo "### Types"
        echo ""
        if [[ ${#types[@]} -gt 0 ]]; then
            for item in "${types[@]}"; do
                IFS='|' read -r name doc <<< "$item"
                if [[ -n "$doc" ]]; then
                    echo "- \`${name}\` — ${doc}"
                else
                    echo "- \`${name}\`"
                fi
            done
        else
            echo "(none)"
        fi
        echo ""

        # Traits section
        echo "### Traits"
        echo ""
        if [[ ${#traits[@]} -gt 0 ]]; then
            for item in "${traits[@]}"; do
                IFS='|' read -r name doc <<< "$item"
                if [[ -n "$doc" ]]; then
                    echo "- \`${name}\` — ${doc}"
                else
                    echo "- \`${name}\`"
                fi
            done
        else
            echo "(none)"
        fi
        echo ""

        # Functions section
        echo "### Functions"
        echo ""
        if [[ ${#fns[@]} -gt 0 ]]; then
            for item in "${fns[@]}"; do
                IFS='|' read -r name doc <<< "$item"
                if [[ -n "$doc" ]]; then
                    echo "- \`fn ${name}\` — ${doc}"
                else
                    echo "- \`fn ${name}\`"
                fi
            done
        else
            echo "(none)"
        fi
        echo ""

        # Module structure
        echo "## Module Structure"
        echo ""
        echo '```'
        for mod in "${modules[@]}"; do
            echo "$mod"
        done
        echo '```'
        echo ""

        # Entry point
        echo "## Entry Point"
        echo ""
        echo "Start reading: \`${entry_point}\`"

    } > "$out_file"

    # Reset arrays for next iteration
    types=()
    traits=()
    fns=()
done

# ──────────────────────────────────────────────────────────────────────────────
# Generate index page
# ──────────────────────────────────────────────────────────────────────────────
echo "Generating atlas index..."

{
    echo "# Crate Atlas"
    echo ""
    echo "Auto-generated overview of all popkern crates. Regenerate with \`tools/generate-atlas.sh\`."
    echo ""
    echo "| Crate | Description | Files | Lines | Tests | Depends on |"
    echo "|-------|-------------|-------|-------|-------|------------|"

    for crate in "${ALL_CRATES[@]}"; do
        desc="${CRATE_DESC[$crate]}"
        files="${CRATE_FILES[$crate]}"
        lines="${CRATE_LINES[$crate]}"
        tests="${CRATE_TESTS[$crate]}"
        deps="${CRATE_DEPS[$crate]}"

        # Shorten deps: strip "popkern-" prefix for brevity
        deps_short=$(echo "$deps" | sed 's/popkern-//g')
        deps_short="${deps_short:-(none)}"

        echo "| [${crate}](${crate}.md) | ${desc} | ${files} | ${lines} | ${tests} | ${deps_short} |"
    done

} > "${ATLAS_DIR}/index.md"

# ──────────────────────────────────────────────────────────────────────────────
# Generate cross-crate type flow map
# ──────────────────────────────────────────────────────────────────────────────
echo "Generating cross-crate type flow map..."

# Core types to track (defined in popkern-core)
CORE_TYPES=(
    BRep
    FaceHandle
    EdgeHandle
    VertexHandle
    ShellHandle
    LoopHandle
    HalfEdgeHandle
    SolidHandle
    SurfaceGeometry
    CurveGeometry
    NurbsSurface
    NurbsCurve
    Tolerance
    Point3
    Vec3
)

# Collect per-type usage data across all crates (excluding popkern-core itself)
declare -A TYPE_CRATES=()   # type -> space-separated list of short crate names
declare -A TYPE_COUNT=()    # type -> count of crates using it

for type_name in "${CORE_TYPES[@]}"; do
    crate_list=""
    count=0
    for crate_dir in crates/popkern-*/; do
        crate=$(basename "$crate_dir")
        if [[ "$crate" == "popkern-core" ]]; then
            continue
        fi
        src_dir="${crate_dir}src"
        if [[ ! -d "$src_dir" ]]; then
            continue
        fi
        # Count files (not lines) that reference this type, excluding test files
        file_count=$(find "$src_dir" -name "*.rs" \
            ! -name "*_tests.rs" \
            ! -name "test_*.rs" \
            ! -path "*/tests/*" \
            2>/dev/null \
            -exec grep -l "\b${type_name}\b" {} + 2>/dev/null \
            | wc -l | tr -d ' ' || true)
        if [[ "${file_count:-0}" -gt 0 ]]; then
            short="${crate#popkern-}"
            crate_list="${crate_list}${short}, "
            (( count++ )) || true
        fi
    done
    TYPE_CRATES["$type_name"]="${crate_list%, }"
    TYPE_COUNT["$type_name"]="$count"
done

# Collect error type data — which crate defines each error, and re-exports
declare -a ERROR_TYPES_LIST=()
declare -A ERROR_DEFINED_IN=()
declare -A ERROR_REEXPORTED_BY=()

# Known error types per crate
declare -A CRATE_ERRORS=(
    [popkern-core]="ValidationError CoreError"
    [popkern-bool]="BoolError"
    [popkern-draft]="DraftError"
    [popkern-blend]="BlendError"
    [popkern-sweep]="SweepError"
    [popkern-shell]="ShellError"
    [popkern-io]="IoError"
    [popkern-mesh]="MeshError"
    [popkern-sketch]="SketchError"
    [popkern-history]="HistoryError ReplayError"
    [popkern-measure]="MeasureError"
    [popkern-constraint]="ConstraintError"
    [popkern-viz]="VizError"
)

for crate in "${!CRATE_ERRORS[@]}"; do
    for err in ${CRATE_ERRORS[$crate]}; do
        ERROR_TYPES_LIST+=("$err")
        ERROR_DEFINED_IN["$err"]="${crate#popkern-}"
        # Check if any other crate re-exports this error
        reexported=""
        for other_dir in crates/popkern-*/; do
            other=$(basename "$other_dir")
            if [[ "$other" == "$crate" ]]; then
                continue
            fi
            if grep -r "pub use.*\b${err}\b" "${other_dir}src" --include="*.rs" -q 2>/dev/null; then
                reexported="${reexported}${other#popkern-}, "
            fi
        done
        ERROR_REEXPORTED_BY["$err"]="${reexported%, }"
    done
done

# Sort ERROR_TYPES_LIST for stable output
IFS=$'\n' sorted_errors=($(printf '%s\n' "${ERROR_TYPES_LIST[@]}" | sort)) ; unset IFS

# Sort CORE_TYPES by usage count descending for the high-impact list
declare -a sorted_types=()
for type_name in "${CORE_TYPES[@]}"; do
    sorted_types+=("${TYPE_COUNT[$type_name]} $type_name")
done
IFS=$'\n' sorted_types=($(printf '%s\n' "${sorted_types[@]}" | sort -rn)) ; unset IFS

# Total crate count (excluding core, for context)
non_core_count=$(( ${#ALL_CRATES[@]} - 1 ))

{
    echo "# Cross-Crate Type Flow"
    echo ""
    echo "Auto-generated map of which core types are used across crate boundaries."
    echo "Regenerate with \`tools/generate-atlas.sh\`."
    echo ""
    echo "Shows: \"If you change type X, these crates are affected.\""
    echo ""
    echo "## Core Types"
    echo ""
    echo "| Type | Defined in | Used by (N crates) | Crates |"
    echo "|------|-----------|---------------------|--------|"
    for type_name in "${CORE_TYPES[@]}"; do
        count="${TYPE_COUNT[$type_name]}"
        crates_str="${TYPE_CRATES[$type_name]:-(none)}"
        echo "| \`${type_name}\` | popkern-core | ${count} | ${crates_str} |"
    done
    echo ""
    echo "## Error Types"
    echo ""
    echo "| Type | Defined in | Re-exported by |"
    echo "|------|-----------|----------------|"
    for err in "${sorted_errors[@]}"; do
        defined="${ERROR_DEFINED_IN[$err]:-(unknown)}"
        reexported="${ERROR_REEXPORTED_BY[$err]:-(none)}"
        echo "| \`${err}\` | popkern-${defined} | ${reexported} |"
    done
    echo ""
    echo "## High-Impact Types"
    echo ""
    echo "Types used by 8+ crates — changes here have the widest blast radius:"
    echo ""
    rank=1
    any_high_impact=0
    for entry in "${sorted_types[@]}"; do
        count="${entry%% *}"
        type_name="${entry#* }"
        if [[ "$count" -ge 8 ]]; then
            echo "${rank}. \`${type_name}\` — ${count} crates"
            (( rank++ )) || true
            any_high_impact=1
        fi
    done
    if [[ "$any_high_impact" -eq 0 ]]; then
        echo "(no types used by 8+ crates)"
    fi
    echo ""
    echo "## All Types by Usage"
    echo ""
    echo "Sorted by number of crates using the type (descending):"
    echo ""
    for entry in "${sorted_types[@]}"; do
        count="${entry%% *}"
        type_name="${entry#* }"
        crates_str="${TYPE_CRATES[$type_name]:-(none)}"
        echo "- \`${type_name}\` — ${count}/${non_core_count} crates: ${crates_str}"
    done

} > "${ATLAS_DIR}/cross-crate.md"

echo "Done. Atlas written to ${ATLAS_DIR}/"
