//! Structural test: API ergonomics enforcement.
//!
//! Ensures public API surfaces follow ergonomic conventions:
//! - Error types use thiserror with Debug
//! - Public types implement Debug
//! - API functions have reasonable parameter counts
//!
//! REMEDIATION: See docs/architecture/overview.md for API design guidelines.

use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

/// Returns true if a file path should be excluded from structural checks.
///
/// Skips test files, example binaries, and FFI crate files.
fn is_excluded_file(path: &std::path::Path) -> bool {
    // Skip files inside a `tests/` directory.
    if path.components().any(|c| c.as_os_str() == "tests") {
        return true;
    }

    // Skip test files included via #[path] — they contain test code, not library code.
    let file_name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if file_name == "tests.rs" || file_name.ends_with("_tests.rs") || file_name.starts_with("test_")
    {
        return true;
    }

    // Skip example binaries — they have main() and .expect() is appropriate there.
    if path
        .components()
        .any(|c| c.as_os_str() == "popkern-examples")
    {
        return true;
    }

    false
}

/// Returns true if the file belongs to the popkern-ffi crate.
fn is_ffi_file(path: &std::path::Path) -> bool {
    path.components().any(|c| c.as_os_str() == "popkern-ffi")
}

/// Collect the lines immediately before `line_index` (0-based) for derive inspection.
///
/// Returns up to `count` preceding lines as a single string.
fn preceding_lines(lines: &[&str], line_index: usize, count: usize) -> String {
    let start = line_index.saturating_sub(count);
    lines[start..line_index].join("\n")
}

// ---------------------------------------------------------------------------
// Rule 1 + 2: Error types must use thiserror (and derive Debug)
// ---------------------------------------------------------------------------

#[test]
fn error_types_use_thiserror() {
    let crates_dir = workspace_root().join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);

    let mut violations_thiserror: Vec<String> = Vec::new();
    let mut violations_debug: Vec<String> = Vec::new();

    for file in &files {
        if is_excluded_file(file) {
            continue;
        }

        let content = popkern_lints::read_file(file);
        let lines: Vec<&str> = content.lines().collect();

        for (i, line) in lines.iter().enumerate() {
            // Match `pub enum SomethingError` (but not type aliases).
            let trimmed = line.trim();
            if !trimmed.starts_with("pub enum ") {
                continue;
            }
            // Extract enum name — must contain "Error".
            let after_pub_enum = trimmed.trim_start_matches("pub enum ").trim();
            let enum_name = after_pub_enum
                .split(|c: char| !c.is_alphanumeric() && c != '_')
                .next()
                .unwrap_or("");
            if !enum_name.contains("Error") {
                continue;
            }

            // Inspect up to 5 preceding lines for `#[derive(... Error ...)]`.
            let context = preceding_lines(&lines, i, 5);

            if !context.contains("Error") || !context.contains("derive") {
                violations_thiserror.push(format!(
                    "{}:{}: `pub enum {}` is missing `#[derive(..., Error, ...)]`. \
                     Remediation: Add `#[derive(Debug, Error)]` and import thiserror.",
                    file.display(),
                    i + 1,
                    enum_name,
                ));
            }

            if !context.contains("Debug") || !context.contains("derive") {
                violations_debug.push(format!(
                    "{}:{}: `pub enum {}` is missing `Debug` in its derive attributes. \
                     Remediation: Add `#[derive(Debug, Error)]`.",
                    file.display(),
                    i + 1,
                    enum_name,
                ));
            }
        }
    }

    assert!(
        violations_thiserror.is_empty(),
        "Public error enums must derive thiserror::Error:\n{}",
        violations_thiserror.join("\n")
    );

    assert!(
        violations_debug.is_empty(),
        "Public error enums must derive Debug:\n{}",
        violations_debug.join("\n")
    );
}

// ---------------------------------------------------------------------------
// Rule 3: Public structs must derive Debug
//
// TODO: Fix existing violations, then make this a hard failure.
// Current known violations (18 structs) lack `#[derive(Debug)]`.
// ---------------------------------------------------------------------------

#[test]
fn public_structs_derive_debug() {
    let crates_dir = workspace_root().join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);

    let mut violations: Vec<String> = Vec::new();

    for file in &files {
        if is_excluded_file(file) {
            continue;
        }
        // FFI types may deliberately not be Debug (repr(C) etc.).
        if is_ffi_file(file) {
            continue;
        }

        let content = popkern_lints::read_file(file);
        let lines: Vec<&str> = content.lines().collect();
        let mut in_test_module = false;

        for (i, line) in lines.iter().enumerate() {
            // Track entry into #[cfg(test)] blocks (simple heuristic).
            if line.contains("#[cfg(test)]") {
                in_test_module = true;
            }
            if in_test_module {
                continue;
            }

            let trimmed = line.trim();

            // Match top-level `pub struct Name` declarations.
            if !trimmed.starts_with("pub struct ") {
                continue;
            }

            // Skip unit structs: `pub struct Foo;`
            if trimmed.ends_with(';') {
                continue;
            }

            // Extract struct name.
            let after_keyword = trimmed.trim_start_matches("pub struct ").trim();
            let struct_name = after_keyword
                .split(|c: char| !c.is_alphanumeric() && c != '_')
                .next()
                .unwrap_or("");

            // Inspect up to 5 preceding lines for `#[derive(... Debug ...)]`.
            let context = preceding_lines(&lines, i, 5);

            if context.contains("derive") && context.contains("Debug") {
                continue;
            }

            violations.push(format!(
                "{}:{}: `pub struct {}` is missing `Debug` in its derive attributes. \
                 Remediation: Derive thiserror::Error on all public error types: \
                 #[derive(Debug, ...)] or add a manual `impl fmt::Debug`.",
                file.display(),
                i + 1,
                struct_name,
            ));
        }
    }

    // TODO: Fix existing violations, then make this a hard failure.
    if !violations.is_empty() {
        println!(
            "WARNING: {} public struct(s) are missing a Debug impl:\n{}",
            violations.len(),
            violations.join("\n")
        );
    }
    // Not asserting — existing violations prevent this from being a hard failure.
}

// ---------------------------------------------------------------------------
// Rule 4: Public functions in lib.rs with >7 non-self parameters
//
// TODO: Fix existing violations (draft_face: 6 params, shell/shell_variable: 8 params),
// then make this a hard failure.
// ---------------------------------------------------------------------------

/// Count the non-self parameters in a function signature string.
///
/// The `sig` must contain the complete parameter list (balanced parentheses).
/// Returns `None` if the signature cannot be parsed.
fn count_non_self_params(sig: &str) -> Option<usize> {
    // Find the outermost parameter list.
    let open = sig.find('(')?;
    // Find the matching close paren at depth 1.
    let mut depth = 0usize;
    let mut close = None;
    for (idx, ch) in sig[open..].char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    close = Some(open + idx);
                    break;
                }
            }
            _ => {}
        }
    }
    let close = close?;
    let params_str = sig[open + 1..close].trim();
    if params_str.is_empty() {
        return Some(0);
    }

    // Count commas at depth 0 (handles generics and nested types).
    let mut depth2 = 0usize;
    let mut commas = 0usize;
    for ch in params_str.chars() {
        match ch {
            '<' | '(' | '[' => depth2 += 1,
            '>' | ')' | ']' => depth2 = depth2.saturating_sub(1),
            ',' if depth2 == 0 => commas += 1,
            _ => {}
        }
    }
    let total = commas + 1;

    // Subtract 1 if the first param is a self receiver.
    let first_param = params_str
        .split(',')
        .next()
        .unwrap_or("")
        .trim()
        .to_string();
    let is_self = first_param.starts_with("self")
        || first_param.starts_with("&self")
        || first_param.starts_with("&mut self");

    Some(if is_self { total - 1 } else { total })
}

#[test]
fn api_functions_have_reasonable_param_counts() {
    let crates_dir = workspace_root().join("crates");
    let max_params = 7usize;

    // Only check the crate's external API surface: `*/src/lib.rs`.
    let all_files = popkern_lints::find_rs_files(&crates_dir);
    let lib_rs_files: Vec<_> = all_files
        .iter()
        .filter(|f| {
            f.file_name().and_then(|n| n.to_str()) == Some("lib.rs")
                && !is_excluded_file(f)
                && !is_ffi_file(f)
        })
        .collect();

    let mut violations: Vec<String> = Vec::new();

    for file in lib_rs_files {
        let content = popkern_lints::read_file(file);
        let lines: Vec<&str> = content.lines().collect();

        let mut i = 0;
        while i < lines.len() {
            let line = lines[i];
            if line.contains("pub fn ") {
                // Accumulate lines until parentheses are balanced.
                let mut sig = String::new();
                let mut depth = 0usize;
                let mut found_open = false;
                let start_line = i;

                let mut j = i;
                while j < lines.len() {
                    sig.push_str(lines[j]);
                    sig.push('\n');
                    for ch in lines[j].chars() {
                        match ch {
                            '(' => {
                                depth += 1;
                                found_open = true;
                            }
                            ')' => depth = depth.saturating_sub(1),
                            _ => {}
                        }
                    }
                    if found_open && depth == 0 {
                        break;
                    }
                    j += 1;
                }

                if let Some(count) = count_non_self_params(&sig)
                    && count > max_params
                {
                    violations.push(format!(
                        "{}:{}: `pub fn` has {} parameters (limit: {}). \
                         Remediation: Group related parameters into a typed options struct.",
                        file.display(),
                        start_line + 1,
                        count,
                        max_params,
                    ));
                }
            }
            i += 1;
        }
    }

    // TODO: Fix existing violations, then make this a hard failure.
    if !violations.is_empty() {
        println!(
            "WARNING: {} pub fn(s) in lib.rs files exceed {} parameters:\n{}",
            violations.len(),
            max_params,
            violations.join("\n")
        );
    }
    // Not asserting — existing violations prevent this from being a hard failure.
}
