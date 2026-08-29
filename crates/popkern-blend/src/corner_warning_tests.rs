//! §5.3 C5 — lib-level tests for
//! [`crate::error::BlendWarning::CornerContinuityDowngraded`]
//! emission. Kept separate from `lib_tests.rs` to stay under the
//! 500-line structural limit.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Trihedral cube corners (N=3) at G1 are inside the Gregory /
/// RollingBallSphere comfort zone — no corner downgrade warning.
#[test]
fn full_cube_g1_emits_no_corner_continuity_downgraded_warning() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let count = r
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CornerContinuityDowngraded { .. }))
        .count();
    assert_eq!(
        count, 0,
        "cube corner at G1 is trihedral and matched — no corner downgrade expected"
    );
}

/// §5.3 C6 — trihedral cube corners with all-convex edges classify
/// as `SameSideConvex` on `CornerPatchInfo::classification`. UI
/// consumers use this to colour corners by configuration-hardness.
#[test]
fn full_cube_corner_patches_classify_as_same_side_convex() {
    use crate::vertex_classifier::VertexBlendClassification;
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    assert_eq!(r.corner_patches.len(), 8, "cube has 8 corners");
    for p in &r.corner_patches {
        assert_eq!(
            p.classification,
            VertexBlendClassification::SameSideConvex,
            "all cube corner edges are convex (outward-rounded)"
        );
    }
}

/// At G2 request every trihedral cube corner now serves the KP20 G2
/// setback fan (PK-225 — the N=3 setback rung is production for
/// corners the sphere does not own; G2 strips are NURBS). The
/// historical G1 downgrade class is gone for this fixture: achieved
/// continuity is G2 and NO downgrade warnings fire.
#[test]
fn full_cube_g2_request_serves_setback_g2_corners_without_downgrade() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let opts = FilletOptions {
        continuity: crate::ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let r = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let corner_warns = r
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CornerContinuityDowngraded { .. }))
        .count();
    assert_eq!(
        corner_warns, 0,
        "all 8 corners serve achieved-G2 setback fans; no downgrade may fire",
    );
    let setback_g2 = r
        .corner_patches
        .iter()
        .filter(|cp| cp.setback_cap && cp.achieved_continuity == crate::ContinuityLevel::G2)
        .count();
    assert!(
        setback_g2 > 0 && setback_g2 == r.corner_patches.len(),
        "every corner patch must be a setback cap at achieved G2, got {}/{}",
        setback_g2,
        r.corner_patches.len(),
    );
}
