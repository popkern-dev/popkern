//! PH-056: a bore's resolution status. bored_block's 0.18 bore in
//! tension (the PH-050 case) reads 0.803, 0.875 and 0.941 of CalculiX's
//! peak at 3, 4 and 6.9 cells per radius; the report must say how low
//! to expect it. The bore table is a worst case over grid shifts and
//! sources (C6 sets rows 2-6, PH-058), so it covers this one-shift
//! shortfall without overstating it by more than 10 points.

use std::path::{Path, PathBuf};

use popkern_analysis_cli::run::run;
use serde_json::Value;

fn repo(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

/// The bored_block case's bore checks at cell size `h` (default budget
/// when `None`).
fn bores(h: Option<f64>) -> Vec<Value> {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("bores");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let mut setup =
        std::fs::read_to_string(repo("tools/analysis/calculix/cases/bored_block.setup.toml"))
            .unwrap_or_else(|e| panic!("{e}"));
    if let Some(h) = h {
        setup += &format!("\n[options]\ntarget_cell_size = {h}\nmax_active_cells = 2000000\n");
    }
    let path = dir.join(format!(
        "{}.toml",
        h.map_or("default".into(), |h| h.to_string())
    ));
    std::fs::write(&path, setup).unwrap_or_else(|e| panic!("{e}"));
    let out = run(
        &repo("tests/fixtures/analysis/gallery/bored_block.step"),
        &path,
    );
    assert_eq!(out.report["status"], "solved");
    out.report["result"]["assessment"]["bores"]
        .as_array()
        .cloned()
        .unwrap_or_default()
}

#[test]
fn the_bore_reports_how_low_its_peak_reads() {
    // (h, CalculiX shortfall measured in PH-050).
    for (h, shortfall) in [(0.06, 0.197), (0.045, 0.125)] {
        let b = bores(Some(h));
        assert_eq!(b.len(), 1, "{b:?}");
        let low = b[0]["status"]["UnderResolved"]["expected_low_by"]
            .as_f64()
            .unwrap_or(f64::NAN);
        assert!(
            low >= shortfall && low <= shortfall + 0.10,
            "h {h}: {low} vs {shortfall}"
        );
    }
    // At 6.9 cells (0.941 of CalculiX) the bore is resolved; the block's
    // convex rounded corners are not bores.
    let b = bores(None);
    assert_eq!(b.len(), 1, "{b:?}");
    assert_eq!(b[0]["status"], "Resolved");
    assert_eq!(
        b[0]["radius"].as_f64().map(|r| (r * 100.0).round()),
        Some(18.0)
    );
}
