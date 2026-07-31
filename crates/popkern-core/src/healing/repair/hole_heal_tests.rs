use super::*;
use crate::tolerance::Tolerance;

#[test]
fn heal_holes_closes_the_open_box_with_the_ruled_fan() {
    let mut brep = BRep::open_box(Tolerance::default()).expect("open box");
    let report = heal_holes(&mut brep, Tolerance::default().linear());

    assert_eq!(report.holes_detected, 1);
    assert_eq!(report.holes_committed, 1);
    assert_eq!(report.attempts.len(), 1);
    let attempt = &report.attempts[0];
    // The G1 PnS complex declines with its measured internal seam gap
    // (recorded, never silent); the ruled fan closes the body.
    assert_eq!(attempt.committed_by, Some(HoleFillerKind::RuledFan));
    assert_eq!(attempt.num_sides, 4);
    assert!(
        attempt
            .pns_decline
            .as_deref()
            .is_some_and(|reason| reason.contains("coherence")),
        "PnS decline must be recorded with its reason: {:?}",
        attempt.pns_decline
    );
    assert!(attempt.coons_decline.is_some());
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
