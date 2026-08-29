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
        // Skip files inside a `tests/` or `integration_tests/`
        // directory — those are test code, not library code (the
        // latter is mounted `#[cfg(test)] mod integration_tests;`).
        let is_test_file = file
            .components()
            .any(|c| c.as_os_str() == "tests" || c.as_os_str() == "integration_tests");
        if is_test_file {
            continue;
        }

        // Skip files inside a `benches/` directory — Criterion benchmarks
        // are not library code, and `.unwrap()` is the standard bench
        // setup pattern (panic-fail is acceptable for a perf harness).
        let is_bench_file = file.components().any(|c| c.as_os_str() == "benches");
        if is_bench_file {
            continue;
        }

        // Skip test files included via `#[cfg(test)] mod <name>;` from a
        // parent module. Several file-name patterns all denote test code:
        //   `tests.rs`                  — inline tests (included via #[path])
        //   `*_tests.rs`                — popkern convention (e.g. `cache_tests.rs`)
        //   `test_*.rs`                 — legacy / single-suite files
        //   `tests_*.rs`                — alternate suffix order (e.g. `tests_dispatch.rs`)
        //   `integration_tests_*.rs`    — popkern-bool crate uses this convention
        let file_name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let is_test_module = file_name == "tests.rs"
            || file_name.contains("_tests")
            || file_name.ends_with("_test_helpers.rs")
            || file_name.starts_with("test_")
            || file_name.starts_with("tests_")
            || file_name.starts_with("integration_tests_");
        if is_test_module {
            continue;
        }

        // Skip example binaries — they have main() and .expect() is appropriate there.
        let is_example = file
            .components()
            .any(|c| c.as_os_str() == "popkern-examples");
        if is_example {
            continue;
        }

        let content = popkern_lints::read_file(file);
        if content.lines().any(|line| line.trim() == "#![cfg(test)]") {
            continue;
        }
        // Brace-depth tracking for #[cfg(test)] items: a pending
        // marker arms on the attribute; the item that follows is
        // skipped — a `mod x;`/`use`/fn-signature line terminated by
        // `;` clears it, a braced item skips until its braces
        // balance. The old latch never reset, leaving everything
        // after a file's FIRST #[cfg(test)] unchecked (this codebase
        // routinely mounts `#[cfg(test)] mod x_tests;` near the top
        // of files, so whole library files went unscanned).
        let mut cfg_test_pending = false;
        let mut skip_depth: i32 = 0;

        for (line_num, line) in content.lines().enumerate() {
            if skip_depth > 0 {
                skip_depth += line.matches('{').count() as i32;
                skip_depth -= line.matches('}').count() as i32;
                continue;
            }
            if line.contains("#[cfg(test)]") {
                cfg_test_pending = true;
                continue;
            }
            if cfg_test_pending {
                // Everything from the attribute to the item's body
                // opening (attribute lines, doc lines, multi-line fn
                // signatures) keeps the marker armed; the item either
                // opens a brace block (skip until balanced) or ends
                // with `;` (a `mod x;` mount — skip that line only).
                let opens = line.matches('{').count() as i32;
                let closes = line.matches('}').count() as i32;
                if opens > closes {
                    skip_depth = opens - closes;
                    cfg_test_pending = false;
                } else if line.trim_end().ends_with(';') {
                    cfg_test_pending = false;
                }
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
