use super::*;
use crate::healing::repair::hole_fill::detect_holes;
use crate::tolerance::Tolerance;

fn open_box_with_hole() -> (BRep, DetectedHole) {
    let brep = BRep::open_box(Tolerance::default()).expect("open box");
    let holes = detect_holes(&brep);
    assert_eq!(holes.len(), 1, "open box has exactly one hole");
    let hole = holes[0].clone();
    (brep, hole)
}

/// The acceptance fixture (PK-224): the open box's 4-sided hole
/// commits the single-patch G1 cap. Planar walls + straight rims
/// make the ribbon construction exact, so the measured mid-band
/// cap↔wall G1 must be near machine zero — this doubles as the
/// calibration pin for the Gregory inner-pair slot mapping (a
/// flipped mapping cannot pass the 5e-2 gate, let alone this pin).
#[test]
fn open_box_commits_single_patch_cap_with_exact_mid_band_g1() {
    let (mut brep, hole) = open_box_with_hole();
    let tolerance = Tolerance::default().linear();
    let commit =
        commit_single_patch_hole_cap(&mut brep, &hole, tolerance).expect("single-patch commit");
    assert!(
        commit.worst_mid_band_g1 <= 1.0e-9,
        "planar walls must measure near-exact mid-band G1, got {:.3e}",
        commit.worst_mid_band_g1
    );
    assert!(
        commit.worst_l2_gap <= tolerance * 10.0,
        "rim containment L2 {:.3e} above budget",
        commit.worst_l2_gap
    );
    assert!(brep.boundary_edges.is_empty(), "the hole must be closed");
    let errors = brep.validate_all();
    assert!(
        errors.is_empty(),
        "capped box must validate clean: {errors:?}"
    );
}

/// The cap must decline transactionally on non-4-sided holes and
/// leave the body untouched.
#[test]
fn non_quad_holes_decline_by_side_count() {
    let (brep, hole) = open_box_with_hole();
    let mut narrowed = hole.clone();
    narrowed.boundary_edges.pop();
    let before = serde_json::to_string(&brep).expect("serializes");
    let mut candidate = brep.clone();
    let error =
        commit_single_patch_hole_cap(&mut candidate, &narrowed, Tolerance::default().linear())
            .expect_err("3 edges must decline");
    assert!(matches!(
        error,
        SinglePatchCapError::UnsupportedSideCount { sides: 3 }
    ));
    assert_eq!(
        serde_json::to_string(&candidate).expect("serializes"),
        before
    );
}
