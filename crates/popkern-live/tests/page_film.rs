//! The film state rules (`web/live/film-state.js`;
//! docs/architecture/presentation-films.md §10): holds, eases, composed
//! moves, the camera, captions and section blending, checked in node.
//! Skips when node is not installed.

use std::path::{Path, PathBuf};
use std::process::Command;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

#[test]
fn film_state_checks_pass() {
    let output = match Command::new("node")
        .arg(root().join("web/live/test/film.test.mjs"))
        .output()
    {
        Ok(output) => output,
        Err(_) => {
            eprintln!("[page_film] SKIP: node is not installed");
            return;
        }
    };
    assert!(
        output.status.success(),
        "film.test.mjs failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let row: serde_json::Value = serde_json::from_slice(&output.stdout).expect("json");
    assert!(row["checks"].as_u64().unwrap_or(0) >= 9, "{row}");
}
