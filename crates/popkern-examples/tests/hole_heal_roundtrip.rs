//! PK-199 acceptance 315: a cube with a face removed round-trips —
//! detect hole → fill → commit → `validate()` → watertight →
//! OCC-clean STEP.
//!
//! `BRep::open_box` is exactly the fixture the acceptance names: a
//! unit cube whose first face was removed, leaving four registered
//! boundary edges with null-sentinel twins. `heal_holes` detects the
//! opening, records the G¹ PnS decline (its sector complex is not
//! fan-committable — measured, see the PK-199 resolution notes) and
//! the Coons decline, and commits the ruled fan cap through the
//! `cap_boundary_ring_fan` Euler operator.

use popkern_core::healing::repair::{HoleFillerKind, heal_holes};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_part::Part;

#[test]
fn healed_open_box_is_watertight_and_occ_clean() {
    let tolerance = Tolerance::default();
    let mut brep = BRep::open_box(tolerance).expect("open box fixture");
    assert_eq!(brep.boundary_edges.len(), 4, "fixture starts open");

    let report = heal_holes(&mut brep, tolerance.linear());
    assert_eq!(report.holes_detected, 1, "one opening detected");
    assert_eq!(report.holes_committed, 1, "the opening commits");
    let attempt = &report.attempts[0];
    assert_eq!(attempt.committed_by, Some(HoleFillerKind::RuledFan));
    assert!(
        attempt.pns_decline.is_some() && attempt.coons_decline.is_some(),
        "G1 declines are recorded, never silent"
    );

    // Strict validation re-arms once boundary_edges empties.
    assert!(brep.boundary_edges.is_empty());
    brep.validate_euler_poincare()
        .expect("strict genus-0 Euler");
    brep.validate().expect("fast validation");

    // Part::from_brep adds the serving-boundary curve-surface gate.
    let part = Part::from_brep(brep, tolerance, "healed open box").expect("serving gates");
    let open_edges = part.open_edge_count().expect("welded mesh");
    assert_eq!(open_edges, 0, "healed body is watertight");

    // STEP + OCC oracle leg (skips with a notice when the OCP venv is
    // absent — the corpus convention).
    let out_dir = std::path::Path::new("../../examples/output/fixtures");
    std::fs::create_dir_all(out_dir).expect("fixture dir");
    let path = out_dir.join("healed_open_box.step");
    part.export_step(&path).expect("STEP export");

    let oracle_python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let oracle_script = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !oracle_python.exists() || !oracle_script.exists() {
        eprintln!("[hole-heal-occ] OCP oracle unavailable — leg skipped");
        return;
    }
    let output = std::process::Command::new(&oracle_python)
        .arg(oracle_script)
        .arg(&path)
        .output()
        .expect("oracle invocation");
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("\"passed\": true"),
        "OCC rejected the healed STEP: {}",
        report.chars().take(400).collect::<String>()
    );
}
