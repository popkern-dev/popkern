//! Viviani's figure eight with the ball's pole at its double point
//! (PK-1139): a ball of radius 10 and a cylinder of radius 5 on x = 5,
//! internally tangent at (10, 0, 0), the ball turned so its pole sits
//! there. Three layers failed in turn: the ball's seam meridian, which
//! ends at the pole, touched the cylinder 7.5e-5 off it (a touch at a
//! curve's own end, located by golden section); the clip pinned the
//! record ends there and re-authored them through 17 samples, 1.5e-3
//! off both surfaces; and with both fixed, the sphere face's chart walk
//! read the lobe halves leaving the pole (all along the meridian in the
//! chart) as duplicate segments and fell back to a 3D outline that kept
//! the lobes in the subtraction, serving an inside-out shell; the
//! union's sphere face, through the pole three times and its antipode
//! once, then meshed 3.0 short and measured 33.7 low (the mesher's pole
//! runs, and the loop normalizer moving a seam side onto its twin).
//! Each operation must give the corpus's closed form, as at every other
//! placement; under the general policy the cut serves glued (its
//! material-sense integral once read the pole runs folded, PK-1145).

use std::f64::consts::{FRAC_PI_2, PI};

use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean_with_options,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn operands(tol: Tolerance) -> (BRep, BRep) {
    let pole_at_double_point = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    )
    .expect("axis");
    let ball = [
        Transform3::rotation_x(-FRAC_PI_2),
        Transform3::rotation_y(PI),
        pole_at_double_point,
    ]
    .iter()
    .try_fold(BRep::sphere(10.0, tol).expect("sphere"), |b, t| {
        b.transform(t)
    })
    .expect("ball");
    let rod = [
        Transform3::rotation_z(PI),
        Transform3::translation(5.0, 0.0, -12.0),
    ]
    .iter()
    .try_fold(BRep::cylinder(5.0, 24.0, tol).expect("cylinder"), |b, t| {
        b.transform(t)
    })
    .expect("rod");
    (ball, rod)
}

fn run(op: BoolOperation, policy: NonManifoldPolicy) -> Result<BRep, BoolError> {
    let tol = Tolerance::default();
    let (ball, rod) = operands(tol);
    boolean_with_options(
        op,
        ball,
        rod,
        tol,
        None,
        &BoolOptions {
            non_manifold: policy,
        },
    )
}

fn assert_volume(body: &BRep, want: f64) {
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let v = popkern_measure::body_volume(body, &Tolerance::default()).expect("volume");
    assert!((v - want).abs() < 1e-6 * want, "volume {v} for {want}");
}

#[test]
fn the_common_and_union_serve_their_closed_forms() {
    let common = run(BoolOperation::Intersection, NonManifoldPolicy::Decline).expect("common");
    assert_volume(&common, 1_205.506_213_504_306_5);
    let union = run(BoolOperation::Union, NonManifoldPolicy::Decline).expect("union");
    assert_volume(&union, 4_868.239_583_435_96);
}

#[test]
fn the_cut_touches_itself_at_the_double_point() {
    match run(BoolOperation::Subtraction, NonManifoldPolicy::Decline) {
        Err(BoolError::NonManifoldResult { contacts }) => {
            let at = Point3::new(10.0, 0.0, 0.0);
            assert!(
                contacts.iter().any(|c| c.location.distance_to(&at) < 1e-6),
                "{contacts:?}"
            );
        }
        other => panic!(
            "expected the contact at the double point, got {:?}",
            other.map(|b| b.faces.len())
        ),
    }
    // The general policy glues the cut where it touches itself: its
    // material-sense integral reads the sphere face's pole runs as the
    // mesher does (PK-1145).
    let general = run(BoolOperation::Subtraction, NonManifoldPolicy::General).expect("general cut");
    assert_volume(&general, 2_983.283_991_282_084_6);
}

/// Viviani's common at kilometre scale with the rod's seam through the
/// double point: its band face's trims are two closed lobe edges whose
/// branches mirror each other in u near the pinch. The quadrature's
/// boundary slivers located chord ends on the mirror branch (0.083 rad
/// away in u, inside a uv band sized by chords 4.8 long in v) and fanned
/// over it: the volume read 2.3e-4 high, its area exact (PK-1146).
#[test]
fn the_kilometre_common_with_the_rods_seam_at_the_double_point_measures_exactly() {
    let base = Tolerance::default();
    let tol = base.with_linear(base.linear() * 1e3);
    let ball = [
        Transform3::rotation_x(-FRAC_PI_2),
        Transform3::rotation_y(PI),
        Transform3::scale_uniform(1e3),
    ]
    .iter()
    .try_fold(BRep::sphere(10.0, base).expect("sphere"), |b, t| {
        b.transform(t)
    })
    .expect("ball");
    let seam_at_double_point =
        Transform3::rotation_around_axis(Point3::new(5.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0), PI)
            .expect("axis");
    let rod = [
        Transform3::rotation_z(PI),
        Transform3::translation(5.0, 0.0, -12.0),
        seam_at_double_point,
        Transform3::scale_uniform(1e3),
    ]
    .iter()
    .try_fold(
        BRep::cylinder(5.0, 24.0, base).expect("cylinder"),
        |b, t| b.transform(t),
    )
    .expect("rod");
    let options = BoolOptions {
        non_manifold: NonManifoldPolicy::Decline,
    };
    let common = boolean_with_options(BoolOperation::Intersection, ball, rod, tol, None, &options)
        .expect("common");
    let want = 1_205.506_213_504_306_5e9;
    let v = popkern_measure::body_volume(&common, &tol).expect("volume");
    assert!((v - want).abs() < 1e-7 * want, "volume {v} for {want}");
}
