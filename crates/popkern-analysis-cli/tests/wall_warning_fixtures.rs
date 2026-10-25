//! PH-055: the thin-wall warning on gallery parts at the default
//! budget. Before: boss_plate warned of a 5.3e-5 mm wall (its sphere
//! corners were tessellated as most of each sphere, PK-787) and
//! tapered_fillet of a 0.004 mm wall (rays beside the acute edge where
//! the tapered fillet meets an end face cross the wedge). Real walls
//! still warn (`tw2_walls_under_one_and_a_half_cells_are_under_resolved`,
//! down to 0.3 cells).
//!
//! The frozen boss_plate was exported before PK-787 fixed the sphere
//! faces' same_sense; on 2026-10-25 its four corner faces took that
//! correction (docs/architecture/analysis-integration.md fixture log),
//! so they mesh as the patches they are and nothing is mismatched. The
//! mismatch warning itself is covered in `model/tessellation_tests.rs`.

use std::path::{Path, PathBuf};

use popkern_analysis_cli::run::run;
use serde_json::Value;

fn repo(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

/// The warnings of the part's PH-050 case at the default budget.
fn warnings(part: &str) -> Vec<Value> {
    let out = run(
        &repo(&format!("tests/fixtures/analysis/gallery/{part}.step")),
        &repo(&format!("tools/analysis/calculix/cases/{part}.setup.toml")),
    );
    assert_eq!(out.report["status"], "solved", "{part}");
    out.report["solve"]["warnings"]
        .as_array()
        .cloned()
        .unwrap_or_default()
}

#[test]
fn no_false_thin_walls_on_wedges_or_mismatched_tessellation() {
    let tapered = warnings("tapered_fillet");
    assert!(tapered.is_empty(), "{tapered:?}");
    let boss = warnings("boss_plate");
    assert!(
        boss.iter()
            .all(|w| w.get("UnderResolved").is_none() && w.get("TessellationMismatch").is_none()),
        "{boss:?}"
    );
}

/// PH-094: centrifugal_impeller's thin curved blades tessellate 1.35%
/// short of the kernel's volume (1.44% of OpenCASCADE's) and the run
/// says so with both volumes;
/// propeller (0.18%, the next largest) does not warn.
#[test]
fn a_tessellation_short_of_the_exact_volume_is_named() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("volume_warning");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("heat.toml");
    let text = "analysis = \"thermal\"\nunits = \"mm\"\nmaterial = \"Aluminium 6061-T6\"\n\
        [[heat_source]]\npower_w = 1\n\
        [[convection]]\nfaces = \"all_other\"\nh = 10\nambient_c = 20\n\
        [options]\nmax_active_cells = 3000\n";
    std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
    let volume_warnings = |part: &str| -> Vec<Value> {
        let out = run(
            &repo(&format!("tests/fixtures/analysis/gallery/{part}.step")),
            &setup,
        );
        assert_eq!(out.report["status"], "solved", "{part}: {}", out.report);
        out.report["solve"]["warnings"]
            .as_array()
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|w| w.get("TessellationVolume").cloned())
            .collect()
    };
    let impeller = volume_warnings("centrifugal_impeller");
    assert_eq!(impeller.len(), 1, "{impeller:?}");
    let (t, e) = (
        impeller[0]["tessellated"].as_f64().unwrap_or(0.0),
        impeller[0]["kernel"].as_f64().unwrap_or(1.0),
    );
    assert!((1.0 - t / e - 0.0135).abs() < 0.001, "{t} of {e}");
    assert!(volume_warnings("propeller").is_empty());
}
