//! A ball resting on a box's top beside the rim of a rod fused to it
//! (PK-1177). The ball (radius 3, its poles turned onto y) touches the
//! box's top at P = (5, 5, 10); a rod of radius 1 fused to the ball runs
//! down into the box, its rim on the top `delta` from P. Within
//! sqrt(2 R tol) of P the ball's curve with the rod and the rod's rim
//! meet under the tolerance, rightly one vertex V, and the ball's face
//! and the box's top share V: the self-contact search took faces sharing
//! a vertex to meet only there, and the union served the ball touching
//! the top at P, 1e-3 from V, as one valid solid with no contact.

use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean_with_options,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Box [0, 10]^3 united with the ball and rod (`delta` the rim's offset
/// from P, beside it when positive), the scene turned `turn` about an
/// axis through P.
fn united(delta: f64, turn: f64) -> (Result<BRep, BoolError>, Transform3) {
    let tol = Tolerance::default();
    let spin = Transform3::rotation_around_axis(
        Point3::new(5.0, 5.0, 10.0),
        Vec3::new(0.3, 1.0, 0.2),
        turn,
    )
    .expect("axis");
    let cube = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        tol,
    )
    .expect("box");
    let ball = BRep::sphere(3.0, tol)
        .and_then(|b| b.transform(&Transform3::rotation_x(std::f64::consts::FRAC_PI_2)))
        .and_then(|b| b.transform(&Transform3::translation(5.0, 5.0, 13.0)))
        .expect("ball");
    let rod = BRep::cylinder(1.0, 5.5, tol)
        .and_then(|b| b.transform(&Transform3::translation(6.0 + delta, 5.0, 8.0)))
        .expect("rod");
    let fused =
        popkern_bool::boolean(BoolOperation::Union, ball, rod, tol, None).expect("ball and rod");
    let options = BoolOptions {
        non_manifold: NonManifoldPolicy::Decline,
    };
    let (cube, fused) = (
        cube.transform(&spin).expect("turn"),
        fused.transform(&spin).expect("turn"),
    );
    (
        boolean_with_options(BoolOperation::Union, cube, fused, tol, None, &options),
        spin,
    )
}

/// Beside the rim by 3e-3 down to 3e-5 (inside sqrt(2 R tol) from 1e-3
/// on), as built and turned: the union declines with its one contact at
/// P (it serves a valid body with no contact without the search of
/// faces facing each other at a shared vertex, from 1e-3 on).
#[test]
fn a_ball_touching_the_top_beside_the_rim_declines_at_its_touch() {
    let tol = Tolerance::default();
    for turn in [0.0, 0.4] {
        for delta in [3e-3, 1e-3, 1e-4, 3e-5] {
            let (outcome, spin) = united(delta, turn);
            let p = spin.apply_point(&Point3::new(5.0, 5.0, 10.0));
            match outcome {
                Err(BoolError::NonManifoldResult { contacts }) => {
                    assert_eq!(contacts.len(), 1, "delta {delta} turn {turn}: {contacts:?}");
                    assert!(
                        contacts[0].location.distance_to(&p) < 10.0 * tol.linear(),
                        "delta {delta} turn {turn}: at {:?}",
                        contacts[0].location
                    );
                }
                other => panic!(
                    "delta {delta} turn {turn}: {}",
                    other.map_or_else(
                        |e| e.to_string(),
                        |b| format!("served {} faces", b.faces.len())
                    )
                ),
            }
        }
    }
}

/// With P inside the rod (the rim beyond it) the ball rests on the rod
/// where it meets the top, and the union is one valid solid of the box,
/// the ball and the rod's part above the box.
#[test]
fn a_ball_over_the_rod_unites_as_one_solid() {
    for turn in [0.0, 0.4] {
        for delta in [-1e-3, -1e-4] {
            let (outcome, _) = united(delta, turn);
            let body = outcome.unwrap_or_else(|e| panic!("delta {delta} turn {turn}: {e}"));
            assert!(body.validate_all().is_empty(), "delta {delta} turn {turn}");
            assert!(body.contacts().is_empty(), "delta {delta} turn {turn}");
        }
    }
}
