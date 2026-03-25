//! Structural test: no .unwrap() or .expect() in library code.
//!
//! Library code in `crates/` must use `?` or typed errors instead of
//! `.unwrap()` or `.expect()`. Test code (#[cfg(test)] modules) is exempt.
//!
//! REMEDIATION: Use `?` operator or return a typed error via thiserror.

use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

#[test]
fn no_unwrap_or_expect_in_library_code() {
    let crates_dir = workspace_root().join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);

    let mut violations = Vec::new();
    let forbidden_patterns = [".unwrap()", ".expect("];

    for file in &files {
        // Skip files inside a `tests/` directory — those are integration tests, not library code.
        let is_test_file = file.components().any(|c| c.as_os_str() == "tests");
        if is_test_file {
            continue;
        }

        // Skip test files included via #[path = "tests.rs"] — they contain test code, not library code.
        let file_name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let is_test_module = file_name == "tests.rs"
            || file_name.ends_with("_tests.rs")
            || file_name.starts_with("test_");
        if is_test_module {
            continue;
        }

        // Skip example binaries — they have main() and .expect() is appropriate there.
        let is_example = file.components().any(|c| c.as_os_str() == "popkern-examples");
        if is_example {
            continue;
        }

        let content = popkern_lints::read_file(file);
        let mut in_test_module = false;

        for (line_num, line) in content.lines().enumerate() {
            // Track #[cfg(test)] module boundaries (simple heuristic)
            if line.contains("#[cfg(test)]") {
                in_test_module = true;
            }

            if in_test_module {
                continue;
            }

            let trimmed = line.trim();
            // Skip comments
            if trimmed.starts_with("//") || trimmed.starts_with("///") {
                continue;
            }

            for pattern in &forbidden_patterns {
                if line.contains(pattern) {
                    violations.push(format!(
                        "{}:{}: {} found. Use `?` operator or return a typed error via thiserror.",
                        file.display(),
                        line_num + 1,
                        pattern
                    ));
                }
            }
        }
    }

    assert!(
        violations.is_empty(),
        "Found .unwrap()/.expect() in library code:\n{}",
        violations.join("\n")
    );
}
