//! The general path serves what the touching shortcut answers: a ring
//! whose bore sits on a shaft, its top flush with the shaft's end, in
//! any rotation (PK-1043).

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;
use crate::options::NonManifoldPolicy;

/// The shaft and the ring on it, turned by `turn`.
fn ring_on_shaft(turn: &Transform3, tol: Tolerance) -> (BRep, BRep) {
    let shaft = BRep::cylinder(5.0, 10.0, tol).expect("shaft");
    let outer = BRep::cylinder(9.0, 6.0, tol).expect("outer");
    let bore = BRep::cylinder(5.0, 8.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
        .expect("bore");
    let ring = crate::boolean(BoolOperation::Subtraction, outer, bore, tol, None)
        .and_then(|r| {
            r.transform(&Transform3::translation(0.0, 0.0, 4.0))
                .map_err(|e| crate::BoolError::MergeFailed {
                    stage: "test".into(),
                    reason: e.to_string(),
                })
        })
        .expect("ring");
    (
        shaft.transform(turn).expect("shaft turned"),
        ring.transform(turn).expect("ring turned"),
    )
}

/// The general path (no shortcut) serves the intersection empty and
/// the subtraction as the shaft, and the union validates, unrotated
/// and turned about each axis. Turned, the rim records closed off the
/// seam and the ring's top walked a disc over its own bore (a hole
/// filler bounded by the records the common blocks put in place of the
/// rim), which the complex declined as a same-domain mismatch and the
/// cellular fallback served non-manifold.
#[test]
fn a_flush_ring_on_a_shaft_serves_on_the_general_path_in_any_rotation() {
    let tol = Tolerance::default();
    for turn in [
        Transform3::identity(),
        Transform3::rotation_z(1.0),
        Transform3::rotation_x(std::f64::consts::PI),
        Transform3::rotation_y(0.4).compose(&Transform3::rotation_x(1.3)),
    ] {
        let (shaft, ring) = ring_on_shaft(&turn, tol);
        let ix = crate::intersect::find_intersections(&shaft, &ring, &tol).expect("intersections");
        for op in [
            BoolOperation::Intersection,
            BoolOperation::Subtraction,
            BoolOperation::Union,
        ] {
            let body = super::super::serve::select_and_serve(
                &shaft,
                &ring,
                &ix,
                op,
                &tol,
                &mut Vec::new(),
                &mut None,
                NonManifoldPolicy::Decline,
            );
            match (op, body) {
                (BoolOperation::Intersection, Ok(b)) => {
                    assert!(b.faces.is_empty(), "{turn:?}: {} faces", b.faces.len())
                }
                (BoolOperation::Subtraction, Ok(b)) => {
                    assert_eq!(b.faces.len(), shaft.faces.len(), "{turn:?}");
                    assert!(
                        b.validate_all().is_empty(),
                        "{turn:?}: {:?}",
                        b.validate_all()
                    );
                }
                (BoolOperation::Union, Ok(b)) => {
                    assert!(
                        b.validate_all().is_empty(),
                        "{turn:?}: {:?}",
                        b.validate_all()
                    );
                }
                (op, r) => panic!(
                    "{op:?} {turn:?}: {:?}",
                    r.map(|b| b.faces.len()).map_err(|e| e.to_string())
                ),
            }
        }
    }
}
