//! Touches inside one shell (PK-1155): a box united with a body that
//! overlaps it at its foot and, elsewhere, rests a ball on the box's top
//! edge or top face. The union is one shell touching itself at the
//! point, where neither operand's topology stands: it served one valid
//! solid, under every policy, the touch unseen. Each case runs with
//! the ball as built (its pole at the touch on the face) and turned
//! about its centre (no vertex near the touch).

use std::f64::consts::{FRAC_1_SQRT_2, PI};

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{Case, Expect, Status, ball, box_at, moved};

/// A prism of the polygon `profile` in the xz plane, 10 deep in y.
fn prism_xz(profile: &[(f64, f64)], tol: Tolerance) -> BRep {
    let points: Vec<Point3> = profile
        .iter()
        .map(|&(x, z)| Point3::new(x, 0.0, z))
        .collect();
    BRep::prism(&points, Vec3::new(0.0, 1.0, 0.0), 10.0, tol).expect("prism")
}

/// The ball of radius 3 at `c`, turned about its centre when `turned`.
fn resting_ball(c: (f64, f64, f64), turned: bool, tol: Tolerance) -> BRep {
    let b = ball(3.0, c, tol);
    if !turned {
        return b;
    }
    moved(
        b,
        &Transform3::rotation_around_axis(
            Point3::new(c.0, c.1, c.2),
            Vec3::new(1.0, 2.0, 3.0),
            0.7,
        )
        .expect("axis"),
    )
}

/// The box [0, 10]³ and, united with it, an L prism (a plate under the
/// box's foot, x 8..12, z 0..2, and a column x 12..20, z 0..20) holding
/// a ball of radius 3 that rests on the box's top edge x = z = 10 at
/// (10, 5, 10) from 45° outside, its centre 3/√2 past both faces.
fn on_edge(turned: bool, tol: Tolerance) -> (BRep, BRep) {
    let d = 3.0 * FRAC_1_SQRT_2;
    let l = prism_xz(
        &[
            (8.0, 0.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (12.0, 20.0),
            (12.0, 2.0),
            (8.0, 2.0),
        ],
        tol,
    );
    let b = boolean(
        BoolOperation::Union,
        l,
        resting_ball((10.0 + d, 5.0, 10.0 + d), turned, tol),
        tol,
        None,
    )
    .expect("the column holds the ball");
    (box_at((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), tol), b)
}

/// The box [0, 10]³ and a C prism (the plate, the column, and a roof x
/// 2..12, z 15..20) holding a ball of radius 3 at (6, 5, 13) that rests
/// on the box's top face at (6, 5, 10), its top cap (above z = 15, one
/// high) in the roof.
fn on_face(turned: bool, tol: Tolerance) -> (BRep, BRep) {
    let c = prism_xz(
        &[
            (8.0, 0.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (2.0, 20.0),
            (2.0, 15.0),
            (12.0, 15.0),
            (12.0, 2.0),
            (8.0, 2.0),
        ],
        tol,
    );
    let b = boolean(
        BoolOperation::Union,
        c,
        resting_ball((6.0, 5.0, 13.0), turned, tol),
        tol,
        None,
    )
    .expect("the roof holds the ball");
    (box_at((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), tol), b)
}

/// The edge case's union. The ball's part outside the column is its cap
/// x < 12, of height h = 3 − (10 + 3/√2 − 12) = 5 − 3/√2, outside the
/// box too (the ball is beyond both faces wherever x < 10). Volume:
/// box 1000 + L 168 × 10 − their overlap 2 × 10 × 2 + the cap
/// πh²(9 − h)/3. Area: the floor 20 × 10; each side the box's 100 and
/// the profile's 168 less their overlap 4; the box's left 100, top 100
/// and right 100 less its foot in the plate 20; the plate's top beyond
/// the box 2 × 10; the column's left 18 × 10 less the disc the ball
/// cuts in x = 12, π(9 − δ²) with δ = 3/√2 − 2; its top 80, right 200;
/// the cap's zone 6πh.
fn on_edge_measures() -> (f64, f64) {
    let h = 5.0 - 3.0 * FRAC_1_SQRT_2;
    let delta = 3.0 * FRAC_1_SQRT_2 - 2.0;
    let volume = 2640.0 + PI * h * h * (9.0 - h) / 3.0;
    let area = 1488.0 + 6.0 * PI * h - PI * (9.0 - delta * delta);
    (volume, area)
}

/// The face case's union. Volume: box 1000 + C 218 × 10 − overlap 40 +
/// the ball less its cap in the roof (height 1: π(9 − 1)/3). Area: the
/// floor 200; each side 100 + 218 − 4; the box's left 100, top 100,
/// right 80; the plate's top beyond the box 20; the column's left
/// (z 2..15) 130, right 200; the top 180; the roof's underside 100 less
/// the disc the ball cuts at z = 15 (radius² 9 − 4) and its left 50; the
/// ball less its cap's zone 6π.
fn on_face_measures() -> (f64, f64) {
    (3140.0 + 100.0 * PI / 3.0, 1788.0 + 25.0 * PI)
}

pub(super) fn one_shell_cases() -> Vec<Case> {
    use BoolOperation::Union;
    use Expect::*;
    let edge = Point3::new(10.0, 5.0, 10.0);
    let face = Point3::new(6.0, 5.0, 10.0);
    let (ve, ae) = on_edge_measures();
    let (vf, af) = on_face_measures();
    let row = |name, build: fn(Tolerance) -> (BRep, BRep), at: Point3, (volume, area)| Case {
        name,
        op: Union,
        build,
        decline: NonManifold(vec![at]),
        general: Some(Body { volume, area }),
        turns: Vec::new(),
        status: Status::Gate,
        general_status: Status::Gate,
    };
    vec![
        row(
            "ball-on-an-edge-inside-one-shell",
            |tol| on_edge(false, tol),
            edge,
            (ve, ae),
        ),
        row(
            "ball-turned-on-an-edge-inside-one-shell",
            |tol| on_edge(true, tol),
            edge,
            (ve, ae),
        ),
        row(
            "ball-on-a-face-inside-one-shell",
            |tol| on_face(false, tol),
            face,
            (vf, af),
        ),
        row(
            "ball-turned-on-a-face-inside-one-shell",
            |tol| on_face(true, tol),
            face,
            (vf, af),
        ),
    ]
}
