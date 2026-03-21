#![warn(missing_docs)]

//! Structural lint helpers for popkern.
//!
//! Provides utilities used by structural tests (`tests/structural_*.rs`)
//! to enforce architectural rules. These are deterministic checks,
//! not LLM-based.

use std::fs;
use std::path::{Path, PathBuf};

/// Recursively find all `.rs` files under a directory.
pub fn find_rs_files(dir: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    if dir.is_dir()
        && let Ok(entries) = fs::read_dir(dir)
    {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                files.extend(find_rs_files(&path));
            } else if path.extension().and_then(|e| e.to_str()) == Some("rs") {
                files.push(path);
            }
        }
    }
    files
}

/// Count the number of lines in a file.
pub fn line_count(path: &Path) -> usize {
    fs::read_to_string(path)
        .map(|s| s.lines().count())
        .unwrap_or(0)
}

/// Read a file's contents as a string.
pub fn read_file(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_default()
}
