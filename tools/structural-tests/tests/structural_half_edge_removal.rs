//! Structural test: no raw `half_edges.remove(` in library code.
//!
//! `Vertex::outgoing` is derived navigational state; removing a
//! half-edge without re-anchoring the vertices that point at it makes
//! their fan walks silently skip them forever (measured on
//! `BRep::torus`: one dangling anchor exempted the seam vertex from
//! the shell vertex recompute, PK-324). Library code must remove
//! half-edges through `BRep::remove_half_edge_reanchored` (or an
//! anchor-aware Euler operator with its own local repair).
//!
//! REMEDIATION: replace `x.half_edges.remove(h)` with
//! `x.remove_half_edge_reanchored(h)`, or add the anchor repair
//! locally and whitelist the file here with a justification.

use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

#[test]
fn half_edge_removal_goes_through_the_reanchoring_path() {
    let crates_dir = workspace_root().join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);

    // Files with their own anchor-aware repair around removal.
    let whitelist = [
        "reanchor.rs", // the sanctioned implementation itself
        "kve.rs",      // local surviving-vertex repair
        "mev.rs",      // rollback path removes freshly created pairs
        "mef.rs",      // kef arm re-points endpoint anchors explicitly
    ];

    let mut violations = Vec::new();
    for file in &files {
        let file_name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let is_test_file = file.components().any(|c| c.as_os_str() == "tests")
            || file_name == "tests.rs"
            || file_name.ends_with("_tests.rs")
            || file_name.starts_with("test_")
            || file_name.starts_with("tests_")
            || file_name.starts_with("integration_tests_");
        if is_test_file || whitelist.contains(&file_name) {
            continue;
        }
        let Ok(content) = std::fs::read_to_string(file) else {
            continue;
        };
        for (idx, line) in content.lines().enumerate() {
            if line.contains("half_edges.remove(") {
                violations.push(format!("{}:{}: {}", file.display(), idx + 1, line.trim()));
            }
        }
    }

    assert!(
        violations.is_empty(),
        "raw half_edges.remove() in library code (use remove_half_edge_reanchored):\n{}",
        violations.join("\n")
    );
}
