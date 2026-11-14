//! A ball resting on a face anywhere is a contact (PK-1160 emergent).
//!
//! Separate shells of a served result touching are found by
//! alternating projection from seeds on both faces. Between a plane and
//! a ball of radius r that iteration moves a seed s off the touch by
//! only s^3 / 2r^2 a round, so the touch was found only when a seed lay
//! on it: the ball's pole vertex (an unturned ball) or the face's
//! interior point (a ball set at the middle). A ball turned and set off
//! the middle of a cube's top united with the cube served as one solid
//! of two shells touching, which every manifold check passes. The
//! touch is now also located as a critical point of the signed gap.

use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean, boolean_with_options,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

fn scene(scale: f64, r: f64, at: (f64, f64), lift: f64, place: &Transform3) -> (BRep, BRep) {
    let tol = Tolerance::default();
    let turn =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), 1.0)
            .expect("axis");
    let s = Transform3::scale(scale, scale, scale);
    let cube = BRep::unit_cube(tol)
        .and_then(|b| b.transform(&s))
        .and_then(|b| b.transform(place))
        .expect("cube");
    let ball = BRep::sphere(r * scale, tol)
        .and_then(|b| b.transform(&turn))
        .and_then(|b| {
            b.transform(&Transform3::translation(
                at.0 * scale,
                at.1 * scale,
                (1.0 + r) * scale + lift,
            ))
        })
        .and_then(|b| b.transform(place))
        .expect("ball");
    (cube, ball)
}

#[test]
fn a_turned_ball_off_the_middle_of_a_face_touches_it() {
    let tol = Tolerance::default();
    let placed =
        Transform3::rotation_around_axis(Point3::new(0.3, 0.1, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    let mut problems = Vec::new();
    for (place_name, place) in [("identity", Transform3::identity()), ("turned", placed)] {
        for scale in [1.0, 40.0] {
            for r in [0.25, 0.5] {
                for at in [(0.5, 0.5), (0.3, 0.6), (0.7, 0.35), (0.42, 0.58)] {
                    let name = format!("{place_name} scale {scale} r {r} at {at:?}");
                    let (cube, ball) = scene(scale, r, at, 0.0, &place);
                    match boolean(BoolOperation::Union, cube.clone(), ball.clone(), tol, None) {
                        Err(BoolError::NonManifoldResult { contacts }) if contacts.len() == 1 => {}
                        Err(e) => problems.push(format!("{name}: declined {e}")),
                        Ok(b) => problems.push(format!(
                            "{name}: served {:?} of {} shells, touching",
                            b.body_kind(),
                            b.shells.len()
                        )),
                    }
                    let general = BoolOptions {
                        non_manifold: NonManifoldPolicy::General,
                    };
                    match boolean_with_options(
                        BoolOperation::Union,
                        cube,
                        ball,
                        tol,
                        None,
                        &general,
                    ) {
                        Ok(b)
                            if b.body_kind() == BodyKind::General
                                && b.contacts().len() == 1
                                && b.validate_all().is_empty() => {}
                        Ok(b) => problems.push(format!(
                            "{name}: general served {:?}, {} contacts",
                            b.body_kind(),
                            b.contacts().len()
                        )),
                        Err(e) => problems.push(format!("{name}: general declined {e}")),
                    }
                }
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}

#[test]
fn a_ball_ten_tolerances_above_the_face_is_clear_of_it() {
    // The near miss: two shells, no contact, at every place.
    let tol = Tolerance::default();
    for at in [(0.5, 0.5), (0.3, 0.6)] {
        let (cube, ball) = scene(1.0, 0.25, at, 10.0 * tol.linear(), &Transform3::identity());
        let body = boolean(BoolOperation::Union, cube, ball, tol, None).expect("served apart");
        assert_eq!(body.body_kind(), BodyKind::Solid, "{at:?}");
        assert!(body.validate_all().is_empty(), "{at:?}");
        assert_eq!(body.solid_regions().len(), 2, "{at:?}");
    }
}
