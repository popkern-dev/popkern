//! Two faces crossing at their tangency (a V configuration, tangency
//! design §4.8, W3c): the intersection's branches cross at the point.

use popkern_bool::{BoolError, BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The cylinder of radius 1 on the z axis and the ball of radius 1.5
/// centred at (-0.5, 0, 0) are tangent at (1, 0, 0), crossing both ways
/// there (across z the cylinder is flat and the ball bends away, round
/// the axis the cylinder bends more).
/// The same faces crossing at their tangency (the cylinder turned so its
/// seam is off the point, the ball about its axis): the exact
/// cylinder x ball arm keeps the two branches through the tangent
/// crossing (it declined loops winding the ball's longitude to the
/// marching path, which emitted fragments round the singular point and
/// every operation's cellular result was invalid). The union and the
/// common lie below the higher and the lower face and serve, their
/// volumes summing to the operands'; the difference is two lobes meeting
/// at the point and declines there.
#[test]
fn a_cylinder_and_ball_crossing_at_their_tangency_serve_and_decline() {
    let tol = Tolerance::default();
    let volume = |b: &BRep| popkern_measure::body_volume(b, &tol).expect("volume");
    for turn in [0.9, 0.0, 2.3] {
        let cylinder = BRep::cylinder(1.0, 4.0, tol)
            .and_then(|c| c.transform(&Transform3::rotation_z(1.0)))
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -2.0)))
            .expect("cylinder");
        let ball = BRep::sphere(1.5, tol)
            .and_then(|s| s.transform(&Transform3::rotation_z(turn)))
            .and_then(|s| s.transform(&Transform3::translation(-0.5, 0.0, 0.0)))
            .expect("ball");
        let union = boolean(
            BoolOperation::Union,
            cylinder.clone(),
            ball.clone(),
            tol,
            None,
        )
        .unwrap_or_else(|e| panic!("turn {turn} union: {e}"));
        let common = boolean(
            BoolOperation::Intersection,
            cylinder.clone(),
            ball.clone(),
            tol,
            None,
        )
        .unwrap_or_else(|e| panic!("turn {turn} common: {e}"));
        for (name, b) in [("union", &union), ("common", &common)] {
            assert!(
                b.validate_all().is_empty(),
                "turn {turn} {name}: {:?}",
                b.validate_all().first()
            );
        }
        let want = volume(&cylinder) + volume(&ball);
        let got = volume(&union) + volume(&common);
        assert!(
            (got - want).abs() < 1e-4 * want,
            "turn {turn}: {got} vs {want}"
        );
        match boolean(BoolOperation::Subtraction, cylinder, ball, tol, None) {
            Err(BoolError::NonManifoldResult { contacts }) => assert!(
                contacts
                    .iter()
                    .any(|c| c.location.distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-6),
                "turn {turn}: {contacts:?}"
            ),
            other => panic!("turn {turn}: {:?}", other.map(|b| b.faces.len())),
        }
    }
}
