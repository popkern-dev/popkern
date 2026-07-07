//! Pentagon-pyramid apex blends, written against the `popkern-part`
//! modeling API — the acceptance gate for the API's phase 3 (see
//! `docs/architecture/part-modeling-api-design.md`).
//!
//! The N = 5 apex is the classic high-valence vertex blend test
//! vehicle. Both continuity registers run: the G1 request (the
//! shipped production path) and the G2 request (curvature-continuous
//! constructions). Blend reports come back through the API instead
//! of being dropped.

use popkern_core::topology::BRep;
use popkern_part::dressup::Continuity;
use popkern_part::prelude::*;

fn blended(continuity: Continuity) -> Result<Part, PartError> {
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, Tolerance::default())
        .map_err(|e| PartError::NotApplicable(format!("pyramid constructor: {e}")))?;
    let part = Part::from_brep(
        pyramid,
        Tolerance::default(),
        "pentagonal pyramid r=1 h=1.2",
    )?;

    // The five slant edges meet at the apex; they are the straight
    // edges that leave the base plane. Radius 0.16 is the PK-097
    // parity-gate value: the KP20 setback cap (production for both
    // continuity levels) accepts it and produces the watertight
    // 20-patch apex fan. Smaller radii (≤ ~0.1) currently make the
    // setback boundary surgery decline ("expected 2 arc/boundary
    // crossings, found 1") and the apex degrades to a mitred hole —
    // that radius-robustness gap is tracked in the PK-100 ledger
    // notes.
    part.fillet_with(
        |f| f.radius(0.16).continuity(continuity),
        |e| e.straight().filter(|i| i.midpoint.z > 1e-6),
    )
}

fn main() {
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("out dir");

    for (label, continuity) in [("g1", Continuity::G1), ("g2", Continuity::G2)] {
        let part =
            blended(continuity).unwrap_or_else(|e| panic!("{label} pentagon blend failed: {e}"));
        let report = part.last_blend_report().expect("blend report");
        println!(
            "pentagon {label}: {} faces, {} edges blended, {} corner patches, {} warnings",
            part.face_count(),
            report.edges_blended,
            report.corner_patches,
            report.warnings.len()
        );
        for w in &report.warnings {
            println!("  warning: {w}");
        }
        let part = part
            .assert_watertight()
            .unwrap_or_else(|e| panic!("{label} watertight gate: {e}"));
        let path = out_dir.join(format!("pentagon_apex_{label}_api.step"));
        part.export_step(&path)
            .unwrap_or_else(|e| panic!("{label} export: {e}"));
        println!("  watertight ✓  wrote {}", path.display());
    }
}
