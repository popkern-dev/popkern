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
/// short as that distance, which the pcurve recompute's uniform samples
/// cannot follow (an edge whose trims sit 1e-2 off its curve), so the
/// union declines; it never serves a wrong body (PK-1176's remaining
/// band, filed as its follow-up).
#[test]
fn a_rod_round_a_turned_balls_pole_declines_or_serves_right() {
    for turn in [0.4, 2.0] {
        for delta in [-2e-6, -1e-5, -1e-4, -1e-3] {
            if let Ok(body) = united(delta, turn) {
                assert_right(&body, delta, turn);
            }
        }
    }
}
