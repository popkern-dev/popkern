//! PK-1153: a ball beside a box's top edge, the plane x sphere circle
//! tangent to that edge from outside, is no crossing; sunk into the
//! edge it is.

use crate::error::ValidationError;
use crate::geometry::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A unit cube and a ball of radius 0.4 beside its edge x = z = 1,
/// `clear` from it along the edge's outward bisector, in one body of
/// two shells, scaled by `scale`.
fn cube_and_ball_beside_edge(clear: f64, scale: f64) -> BRep {
    let tol = Tolerance::default();
    let r = 0.4;
    let d = (r + clear / scale) / std::f64::consts::SQRT_2;
    let mut body = BRep::unit_cube(tol).expect("cube");
    let ball = BRep::sphere(r, tol)
        .expect("sphere")
        .transform(&Transform3::translation(1.0 + d, 0.5, 1.0 + d))
        .expect("moves");
    body.append_body(&ball);
    body.transform(&Transform3::scale_uniform(scale))
        .expect("scales")
}

#[test]
fn a_circle_touching_a_face_edge_from_outside_is_no_crossing() {
    for scale in [1.0, 1.0e-3, 1.0e3] {
        let body = cube_and_ball_beside_edge(2.0e-6, scale);
        assert!(
            body.validate_self_intersection().is_ok(),
            "scale {scale}: {:?}",
            body.validate_self_intersection()
        );
    }
}

#[test]
fn a_ball_sunk_into_the_edge_still_crosses() {
    for scale in [1.0, 1.0e-3, 1.0e3] {
        let body = cube_and_ball_beside_edge(-1.0e-2 * scale, scale);
        assert!(
            matches!(
                body.validate_self_intersection(),
                Err(ValidationError::SelfIntersection { .. })
            ),
            "scale {scale}"
        );
    }
}
