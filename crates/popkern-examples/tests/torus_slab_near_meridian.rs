//! A torus less a slab beside the meridian edges of its patches
//! (PK-1173). `BRep::torus` splits its tube at the meridian planes x = 0
//! and y = 0; a slab across the tube there cuts it in two. Within 0.023
//! of a meridian plane the plane × torus solver read each slab face as
//! missing the torus, and the cut served the torus unchanged; with both
//! slab faces on one side the walk filed the strip between the meridian
//! edge and the nearer face as a hole of the far piece, and the cut
//! declined (after first serving half the torus).

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const BIG_R: f64 = 5.0;
const SMALL_R: f64 = 1.5;

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

/// `∫ f` over `[a, b]` on `panels` Gauss–Legendre panels.
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

/// The area of the torus's section by the plane at `x` (|x| under
/// R − r): both ovals, `4 ∫ r² cos²θ ρ / √(ρ² − x²) dθ` with
/// `ρ = R + r sin θ` over the tube's half-turn (smooth in θ).
fn section_area(x: f64) -> f64 {
    let pi = std::f64::consts::PI;
    4.0 * integrate(-0.5 * pi, 0.5 * pi, 8, |t| {
        let rho = BIG_R + SMALL_R * t.sin();
        SMALL_R * SMALL_R * t.cos() * t.cos() * rho / (rho * rho - x * x).sqrt()
    })
}

/// The torus less the slab `x0 ≤ x ≤ x1` (|y| ≤ 10, |z| ≤ 5), both
/// turned `turn` about the axis, and the volume the slab removes.
fn cut(x0: f64, x1: f64, turn: f64) -> (BRep, f64) {
    let tol = Tolerance::default();
    let spin = Transform3::rotation_z(turn);
    let ring = BRep::torus(BIG_R, SMALL_R, tol)
        .and_then(|b| b.transform(&spin))
        .expect("torus");
    let slab = BRep::prism(
        &[
            Point3::new(x0, -10.0, -5.0),
            Point3::new(x1, -10.0, -5.0),
            Point3::new(x1, 10.0, -5.0),
            Point3::new(x0, 10.0, -5.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        tol,
    )
    .and_then(|b| b.transform(&spin))
    .expect("slab");
    let body = boolean(BoolOperation::Subtraction, ring, slab, tol, None).expect("served");
    (body, integrate(x0, x1, 4, section_area))
}

/// Slabs across the meridian edges, beside them on either side, and
/// both faces on one side, as built and turned off and back onto the
/// patches' meridians: each serves a valid body, the
/// slab's volume gone to 1e-9 of the torus.
#[test]
fn a_slab_beside_the_meridian_edges_is_subtracted() {
    let tol = Tolerance::default();
    let torus = 2.0 * std::f64::consts::PI.powi(2) * BIG_R * SMALL_R * SMALL_R;
    let slabs = [
        (-0.005, 0.005),
        (-0.0005, 0.0005),
        (-0.004, 0.006),
        (-0.001, 0.009),
        (0.001, 0.011),
        (-0.011, -0.001),
        (0.2, 0.21),
    ];
    for turn in [0.0, 0.4, std::f64::consts::FRAC_PI_2] {
        for (x0, x1) in slabs {
            let (body, removed) = cut(x0, x1, turn);
            let issues = body.validate_all();
            assert!(issues.is_empty(), "[{x0}, {x1}] turn {turn}: {issues:?}");
            let v = popkern_measure::body_volume(&body, &tol).expect("volume");
            let gone = torus - v;
            assert!(
                (gone - removed).abs() < 1e-9 * torus,
                "[{x0}, {x1}] turn {turn}: removed {gone}, the slab holds {removed}"
            );
        }
    }
}
