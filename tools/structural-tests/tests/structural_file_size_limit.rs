//! Structural test: no source file exceeds 500 lines.
//!
//! Enforcement is a RATCHET: `file_size_baseline.txt` freezes the
//! pre-existing violators at their current line counts (PK-118
//! containment — the flat rule had been red for the whole legacy set,
//! and the binary sat beyond the sweep's fail-fast horizon so it
//! enforced nothing). A file fails when it exceeds 500 lines and is
//! either absent from the baseline or LARGER than its frozen count.
//! Shrinking a file below its entry is always legal (lower or delete
//! the entry); raising an entry is not — split into submodules
//! instead. See docs/architecture/ for module boundaries.

use std::collections::HashMap;
use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

/// Frozen line counts for pre-existing violators, keyed by
/// workspace-relative path with forward slashes.
fn ratchet_baseline() -> HashMap<String, usize> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("file_size_baseline.txt");
    let text = std::fs::read_to_string(&path).expect("file_size_baseline.txt must be readable");
    let mut map = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.rsplitn(2, ' ');
        let count: usize = parts
            .next()
            .and_then(|c| c.parse().ok())
            .expect("baseline line must end with a line count");
        let file = parts.next().expect("baseline line must name a file");
        map.insert(file.to_string(), count);
    }
    map
}

#[test]
fn no_file_exceeds_500_lines() {
    let root = workspace_root();
    let crates_dir = root.join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);
    let baseline = ratchet_baseline();

    let mut violations = Vec::new();
    let limit = 500;

    for file in &files {
        let count = popkern_lints::line_count(file);
        if count <= limit {
            continue;
        }
        let rel = file
            .strip_prefix(&root)
            .unwrap_or(file)
            .to_string_lossy()
            .replace(std::path::MAIN_SEPARATOR, "/");
        match baseline.get(rel.as_str()) {
            Some(&frozen) if count <= frozen => {}
            Some(&frozen) => violations.push(format!(
                "{rel}: {count} lines, GREW past its frozen ratchet entry ({frozen}). \
                 Split into submodules instead of raising the baseline."
            )),
            None => violations.push(format!(
                "{rel}: {count} lines (limit: {limit}). New files must stay under the \
                 limit — split into submodules. See docs/architecture/ for module boundaries."
            )),
        }
    }

    assert!(
        violations.is_empty(),
        "Files exceeding {} lines:\n{}",
        limit,
        violations.join("\n")
    );
}
