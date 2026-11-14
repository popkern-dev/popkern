//! A ball fused to a rod whose wall passes beside the ball's chart
//! pole (PK-1176). `BRep::sphere` has its poles on z; a rod of radius 1
//! along z with its wall `delta` from the south pole meets the ball in
//! a loop passing that close to the pole (beside it for `delta > 0`,
//! round it for `delta < 0`). Within 3e-4 of the pole the union served
//! the rod alone (the ball lost), meshed the ball as slivers, or
//! declined, by three defects: the cylinder × sphere scan handed the
//! loop to the marching path, the chart walk refused the seam where
//! the loop crossed it beside the pole, and the post-merge spur pass
//! retired the seam's stub to the pole with the pole.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Gauss–Legendre nodes and weights on [−1, 1], 16 points.
const GL16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_5),
    (0.617_876_244_402_643_7, 0.149_595_988_816_576_7),
    (0.755_404_408_355_003, 0.124_628_971_255_533_9),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_09),
];

fn integrate(a: f64, b: f64, panels: usize, f: impl Fn(f64) -> f64) -> f64 {
    let h = (b - a) / panels as f64;
    (0..panels)
        .map(|k| {
            let (mid, half) = (a + h * (k as f64 + 0.5), 0.5 * h);
            GL16.iter()
                .map(|&(x, w)| w * half * (f(mid - half * x) + f(mid + half * x)))
                .sum::<f64>()
        })
        .sum()
}

/// The union's volume: the ball (radius 3) and the rod (radius 1, z 8
/// to 13.5, its axis `1 + delta` from the ball's) less their common
/// part, which over the rod's disc runs from the ball's underside to
/// the rod's top: `∫ (0.5 + √(9 − ρ²)) dA`, `ρ` the distance to the
/// ball's axis (under 3 across the disc).
fn union_volume(delta: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let d = 1.0 + delta;
    let common = integrate(0.0, 1.0, 4, |r| {
        r * integrate(0.0, 2.0 * pi, 8, |t| {
            let rho2 = d * d + r * r + 2.0 * d * r * t.cos();
            0.5 + (9.0 - rho2).sqrt()
        })
    });
    4.0 / 3.0 * pi * 27.0 + pi * 5.5 - common
}

fn united(delta: f64, turn: f64) -> Result<BRep, popkern_bool::BoolError> {
    let tol = Tolerance::default();
    let spin = Transform3::rotation_around_axis(
        popkern_core::geometry::Point3::new(5.0, 5.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        turn,
    )
    .expect("axis");
    let ball = BRep::sphere(3.0, tol)
        .and_then(|b| b.transform(&Transform3::translation(5.0, 5.0, 13.0)))
        .and_then(|b| b.transform(&spin))
        .expect("ball");
    let rod = BRep::cylinder(1.0, 5.5, tol)
        .and_then(|b| b.transform(&Transform3::translation(6.0 + delta, 5.0, 8.0)))
        .and_then(|b| b.transform(&spin))
        .expect("rod");
    boolean(BoolOperation::Union, ball, rod, tol, None)
}

fn assert_right(body: &BRep, delta: f64, turn: f64) {
    let tol = Tolerance::default();
    let issues = body.validate_all();
    assert!(issues.is_empty(), "delta {delta} turn {turn}: {issues:?}");
    let v = popkern_measure::body_volume(body, &tol).expect("volume");
    let want = union_volume(delta);
    assert!(
        (v - want).abs() < 1e-6 * want,
        "delta {delta} turn {turn}: volume {v} for {want}"
    );
}

/// The rod's wall from 1e-2 down to 2e-6 beside the pole (as built,
/// the seam through the wall's nearest point, and turned about z, the
/// pole on the axis and the seam turned off it), and from 2e-6 to 1e-3
/// round it as built, and 0.3 and 0.5 round it turned: every union is
/// valid, its volume the closed form's to 1e-6 (the fine measuring
/// mesh's resolution).
#[test]
fn a_rod_beside_a_balls_pole_unites_with_it() {
    let beside = [1e-2, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 2e-6];
    let round = [-2e-6, -3e-6, -1e-5, -1e-4, -1e-3];
    for delta in beside.iter().chain(&round) {
        let body = united(*delta, 0.0).unwrap_or_else(|e| panic!("delta {delta}: {e}"));
        assert_right(&body, *delta, 0.0);
    }
    for turn in [0.4, 2.0] {
        for delta in beside.iter().chain(&[-0.3, -0.5]) {
            let body =
                united(*delta, turn).unwrap_or_else(|e| panic!("delta {delta} turn {turn}: {e}"));
            assert_right(&body, *delta, turn);
        }
    }
}

/// A rod round a turned ball's pole within 1e-3 of it: the ring's
/// sphere pcurve sweeps half a turn of longitude within a stretch as
/// short as that distance. The pcurve recompute's uniform samples could
/// not follow it (an edge whose trims sat 1e-2 off its curve) and the
/// union declined; refined locally where its image strays, the pcurve
/// follows the sweep and every union serves valid, its volume the
/// closed form's to 1e-6 (PK-1178: the measuring mesh splits the ring's
/// sweep from its curve's own points and reads only true poles as
/// poles).
#[test]
fn a_rod_round_a_turned_balls_pole_unites_with_it() {
    for turn in [0.4, 2.0] {
        for delta in [-2e-6, -1e-5, -1e-4, -1e-3] {
            let body =
                united(delta, turn).unwrap_or_else(|e| panic!("delta {delta} turn {turn}: {e}"));
            assert_right(&body, delta, turn);
        }
    }
}

/// The ball alone turned by `seam` about its axis (its chart seam at
/// `seam` from the rod's side of the pole), the rod in place: the turns
/// above move ball and rod together and so leave the seam at 0, −0.4
/// and −2.0 from the rod.
fn united_at_seam(delta: f64, seam: f64) -> Result<BRep, popkern_bool::BoolError> {
    let tol = Tolerance::default();
    let ball = BRep::sphere(3.0, tol)
        .and_then(|b| b.transform(&Transform3::rotation_z(seam)))
        .and_then(|b| b.transform(&Transform3::translation(5.0, 5.0, 13.0)))
        .expect("ball");
    let rod = BRep::cylinder(1.0, 5.5, tol)
        .and_then(|b| b.transform(&Transform3::translation(6.0 + delta, 5.0, 8.0)))
        .expect("rod");
    boolean(BoolOperation::Union, ball, rod, tol, None)
}

/// Sixteen seam turns, none on a symmetry of the scene.
fn seams() -> impl Iterator<Item = f64> {
    (0..16).map(|k| k as f64 * std::f64::consts::FRAC_PI_8 + 0.05)
}

/// Round the pole at every seam turn (PK-1185): the ring's sphere
/// pcurve was left 1e-2 off where uniform doubling could not move it
/// (refined locally only after a doubling helped, now also on a chart
/// with a singular row, its stall rule waiting while a straying
/// interval still turns past a quarter period), and the chart walk
/// misread the seam beside the pole (its border window's slack a chart
/// fraction where the chart is stretched; a run end keyed by a sample
/// the sweep had turned half a period) and kept the polar cap or
/// declined. Every union is valid, its volume the closed form's to 1e-6.
#[test]
fn a_rod_round_a_balls_pole_unites_at_every_seam_turn() {
    for delta in [-2e-6, -1e-4, -1e-3] {
        for seam in seams() {
            let body = united_at_seam(delta, seam)
                .unwrap_or_else(|e| panic!("delta {delta} seam {seam}: {e}"));
            assert_right(&body, delta, seam);
        }
    }
}

/// Beside the pole at every seam turn (PK-1185): the seam crossing the
/// ring both beside the pole and far off misfit its border window (the
/// window's own bound read at the stretched end) and kept the cap, or
/// placed the arc between the crossings a period off on a balanced but
/// wrong branch and served an invalid body. Every union is valid with
/// the closed-form areas of its sphere and cylinder faces to 1e-7
/// relative (the volume the measuring mesh reads beside the pole is
/// PK-1186's).
#[test]
fn a_rod_beside_a_balls_pole_unites_at_every_seam_turn() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    for delta in [1e-3, 1e-4, 1e-5, 3e-6, 2e-6] {
        let d = 1.0 + delta;
        let rho2 = |r: f64, t: f64| d * d + r * r + 2.0 * d * r * t.cos();
        let sphere = 36.0 * pi
            - integrate(0.0, 1.0, 4, |r| {
                r * integrate(0.0, 2.0 * pi, 8, |t| 3.0 / (9.0 - rho2(r, t)).sqrt())
            });
        let wall = integrate(0.0, 2.0 * pi, 8, |t| 5.0 - (9.0 - rho2(1.0, t)).sqrt());
        for seam in seams() {
            let body = united_at_seam(delta, seam)
                .unwrap_or_else(|e| panic!("delta {delta} seam {seam}: {e}"));
            let issues = body.validate_all();
            assert!(issues.is_empty(), "delta {delta} seam {seam}: {issues:?}");
            for (fh, f) in body.faces.iter() {
                let want = match f.surface.as_ref().map(|s| &s.geometry) {
                    Some(SurfaceGeometry::Sphere { .. }) => sphere,
                    Some(SurfaceGeometry::Cylinder { .. }) => wall,
                    _ => pi,
                };
                let area = popkern_measure::face_area(&body, fh, &tol).expect("area");
                assert!(
                    (area - want).abs() < 1e-7 * want,
                    "delta {delta} seam {seam}: area {area} for {want}"
                );
            }
        }
    }
}
