use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval as _;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{LoopSide, loop_side_exact, point_in_loop_exact};

/// A box [-10, 10]² × [-5, 0] united with a ball of radius 3 at the
/// origin: its top face holds the ball's equator as a hole.
fn holed_top() -> (BRep, popkern_core::topology::FaceHandle) {
    let tol = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(-10.0, -10.0, z),
            Point3::new(10.0, -10.0, z),
            Point3::new(10.0, 10.0, z),
            Point3::new(-10.0, 10.0, z),
        ]
    };
    let bx = BRep::polygon_frustum(&rect(-5.0), &rect(0.0), tol).expect("box");
    let ball = BRep::sphere(3.0, tol).expect("ball");
    let body = crate::boolean(crate::BoolOperation::Union, bx, ball, tol, None).expect("union");
    let face = body
        .faces
        .iter()
        .find(|(_, f)| !f.inner_loops.is_empty())
        .map(|(h, _)| h)
        .expect("holed face");
    (body, face)
}

#[test]
fn a_hole_rim_tells_inside_beyond_the_tolerance_from_on_the_rim() {
    let (body, face) = holed_top();
    let f = &body.faces[face];
    let hole = f.inner_loops[0];
    let surface = &f.surface.as_ref().expect("surface").geometry;
    let side_at = |r: f64| {
        let angle: f64 = 1.0;
        let p = Point3::new(r * angle.cos(), r * angle.sin(), 0.0);
        let (u, v) = surface.closest_point(&p).expect("uv");
        (
            loop_side_exact(&body, hole, u, v),
            point_in_loop_exact(&body, hole, u, v),
        )
    };
    // Two tolerances inside the rim: in the hole, though the rim's
    // polygon chords (sagitta 3.6e-3 at 64 stations) pass outside it.
    assert_eq!(side_at(3.0 - 2e-6), (Some(LoopSide::Inside), Some(true)));
    assert_eq!(side_at(3.0 - 0.5e-6), (Some(LoopSide::OnRim), Some(true)));
    assert_eq!(side_at(3.0 + 0.5e-6), (Some(LoopSide::OnRim), Some(true)));
    assert_eq!(side_at(3.0 + 2e-6), (Some(LoopSide::Outside), Some(false)));
    assert_eq!(side_at(1.0), (Some(LoopSide::Inside), Some(true)));
    assert_eq!(side_at(5.0), (Some(LoopSide::Outside), Some(false)));
}

#[test]
fn a_straight_hole_rim_is_on_the_rim_within_the_tolerance() {
    // A block [0, 2]³ less a unit cube from its top: the top face holds
    // the square [0.5, 1.5]² as a hole. A point on a side of the square
    // is on the rim whichever way the crossing count rounds it (the
    // hole filter dropped the refilled pocket's rim pieces as inside).
    let tol = Tolerance::default();
    let cube = |at: f64, s: f64, z: f64| {
        let rect = |h: f64| {
            vec![
                Point3::new(at, at, h),
                Point3::new(at + s, at, h),
                Point3::new(at + s, at + s, h),
                Point3::new(at, at + s, h),
            ]
        };
        BRep::polygon_frustum(&rect(z), &rect(z + s), tol).expect("cube")
    };
    let body = crate::boolean(
        crate::BoolOperation::Subtraction,
        cube(0.0, 2.0, 0.0),
        cube(0.5, 1.0, 1.0),
        tol,
        None,
    )
    .expect("pocket");
    let (_, f) = body
        .faces
        .iter()
        .find(|(_, f)| !f.inner_loops.is_empty())
        .expect("holed face");
    let hole = f.inner_loops[0];
    let surface = &f.surface.as_ref().expect("surface").geometry;
    let side_at = |x: f64, y: f64| {
        let (u, v) = surface.closest_point(&Point3::new(x, y, 2.0)).expect("uv");
        loop_side_exact(&body, hole, u, v)
    };
    for (x, y) in [(0.5, 1.0), (1.0, 0.5), (1.5, 1.2), (0.7, 1.5), (0.5, 0.5)] {
        assert_eq!(side_at(x, y), Some(LoopSide::OnRim), "({x}, {y})");
    }
    assert_eq!(side_at(0.5 + 0.5e-6, 1.0), Some(LoopSide::OnRim));
    assert_eq!(side_at(0.5 + 2e-6, 1.0), Some(LoopSide::Inside));
    assert_eq!(side_at(0.5 - 2e-6, 1.0), Some(LoopSide::Outside));
    assert_eq!(side_at(1.0, 1.0), Some(LoopSide::Inside));
}
