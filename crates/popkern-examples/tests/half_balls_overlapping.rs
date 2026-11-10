//! Two half balls of one sphere, charted about perpendicular axes,
//! overlapping on a quarter ball (PK-1138).
//!
//! Each half's sphere face is charted about its own pole axis. The
//! union's assembly kept one copy of the common quarter, and the
//! post-processing merged co-surfacial faces across the two charts: the
//! merged face enclosed the survivor's pole and crossed its seam with
//! no seam edge, its pcurves 0.1 off its edges, and the boolean
//! declined (and then misreported the operands as only touching). Faces
//! merge only within one chart. The half ball's hemisphere, bounded by
//! its equator alone, measured no area: a loop along one chart row read
//! as a uv rectangle of no height.

use std::f64::consts::{FRAC_PI_2, PI};

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn slab(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(lo.0, lo.1, z),
            Point3::new(hi.0, lo.1, z),
            Point3::new(hi.0, hi.1, z),
            Point3::new(lo.0, hi.1, z),
        ]
    };
    BRep::polygon_frustum(&rect(lo.2), &rect(hi.2), tol).expect("slab")
}

/// The half ball `z >= 0` and the half ball `x >= 0` (its sphere turned
/// a quarter about y, so its poles lie on x), radius 5.
fn half_balls(tol: Tolerance) -> (BRep, BRep) {
    let ball = || BRep::sphere(5.0, tol).expect("ball");
    let top = boolean(
        BoolOperation::Subtraction,
        ball(),
        slab((-6.0, -6.0, -6.0), (6.0, 6.0, 0.0), tol),
        tol,
        None,
    )
    .expect("top half");
    let turn = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        FRAC_PI_2,
    )
    .expect("turn");
    let side = boolean(
        BoolOperation::Subtraction,
        ball().transform(&turn).expect("turned"),
        slab((-6.0, -6.0, -6.0), (0.0, 6.0, 6.0), tol),
        tol,
        None,
    )
    .expect("side half");
    (top, side)
}

fn measures(body: &BRep, tol: &Tolerance) -> (f64, f64) {
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let volume = popkern_measure::body_volume(body, tol).expect("volume");
    let area: f64 = body
        .faces
        .keys()
        .map(|f| popkern_measure::face_area(body, f, tol).expect("area"))
        .sum();
    (volume, area)
}

fn close(got: f64, want: f64) -> bool {
    (got - want).abs() <= 1e-9 * want
}

#[test]
fn the_union_is_three_quarters_of_the_ball() {
    let tol = Tolerance::default();
    let (top, side) = half_balls(tol);
    let body = boolean(BoolOperation::Union, top, side, tol, None).expect("serves");
    let (volume, area) = measures(&body, &tol);
    // Three quarters of the ball; three quarters of its sphere and two
    // half discs.
    assert!(close(volume, 125.0 * PI), "volume {volume}");
    assert!(close(area, 100.0 * PI), "area {area}");
}

#[test]
fn a_half_balls_hemisphere_measures_half_the_sphere() {
    let tol = Tolerance::default();
    let (top, side) = half_balls(tol);
    for half in [top, side] {
        let (volume, area) = measures(&half, &tol);
        assert!(close(volume, 250.0 * PI / 3.0), "volume {volume}");
        assert!(close(area, 75.0 * PI), "area {area}");
    }
}
