//! Structural test: no source file exceeds 500 lines.
//!
//! REMEDIATION: Split into submodules. See docs/architecture/ for module boundaries.

use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

#[test]
fn no_file_exceeds_500_lines() {
    let crates_dir = workspace_root().join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);

    let mut violations = Vec::new();
    let limit = 500;

    for file in &files {
        let count = popkern_lints::line_count(file);
        if count > limit {
            violations.push(format!(
                "{}: {} lines (limit: {}). Split into submodules. See docs/architecture/ for module boundaries.",
                file.display(),
                count,
                limit
            ));
        }
    }

    assert!(
        violations.is_empty(),
        "Files exceeding {} lines:\n{}",
        limit,
        violations.join("\n")
    );
}
