use super::*;
use crate::tolerance::Tolerance;

#[test]
fn heal_holes_closes_the_open_box_with_the_single_patch_g1_cap() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let report = heal_holes(&mut brep, Tolerance::default().linear());

    assert_eq!(report.holes_detected, 1);
    assert_eq!(report.holes_committed, 1);
    assert_eq!(report.attempts.len(), 1);
    let attempt = &report.attempts[0];
    // PK-224: the 4-sided hole serves the single-patch G1 cap — no
    // internal seams, ribbon-tangent to every wall. The PnS/Coons/
    // ruled rungs are the recorded fallbacks and never fire here.
    assert_eq!(attempt.committed_by, Some(HoleFillerKind::SinglePatchG1));
    assert_eq!(attempt.num_sides, 4);
    assert!(attempt.single_patch_decline.is_none());
    assert!(attempt.pns_decline.is_none());
    assert!(attempt.coons_decline.is_none());
    assert!(attempt.ruled_decline.is_none());

    assert!(brep.boundary_edges.is_empty());
    brep.validate_euler_poincare().expect("strict Euler holds");
    brep.validate().expect("full fast validation holds");
}

#[test]
fn heal_holes_on_a_closed_body_is_a_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let before = serde_json::to_string(&brep).expect("serializes");
    let report = heal_holes(&mut brep, Tolerance::default().linear());
    assert_eq!(report.holes_detected, 0);
    assert_eq!(report.holes_committed, 0);
    assert!(report.attempts.is_empty());
    assert_eq!(serde_json::to_string(&brep).expect("serializes"), before);
}
