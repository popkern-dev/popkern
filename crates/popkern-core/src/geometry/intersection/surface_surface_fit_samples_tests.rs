//! Fit-ladder pins for [`fit_march_samples`] (PK-673).

use super::fit_march_samples;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::marching::MarchSample;

const RADIUS: f64 = 1.229_367_061_317_363;
const PITCH: f64 = 0.5;
/// Axial ripple on the helix: amplitude and waves per turn. Enough
/// detail that the 1/5 and 2/5 least-squares rungs miss the residual
/// target and the ladder climbs, as it does on built thread walls.
const RIPPLE: f64 = 2.0e-4;
const WAVES: f64 = 24.0;

/// The rippled helix at angle `theta`.
fn curve_at(theta: f64) -> Point3 {
    Point3::new(
        RADIUS * theta.cos(),
        RADIUS * theta.sin(),
        PITCH * theta / std::f64::consts::TAU + RIPPLE * (WAVES * theta).sin(),
    )
}

/// Distance from `q` to the rippled helix near angle `hint` (dense
/// bracket, then golden-free ternary refinement; the curve is smooth
/// and the hint is within a sample interval).
fn curve_distance(q: &Point3, hint: f64) -> f64 {
    let mut best = (f64::MAX, hint);
    for k in -200..=200 {
        let t = hint + k as f64 * 1.0e-4;
        let d = curve_at(t).distance_to(q);
        if d < best.0 {
            best = (d, t);
        }
    }
    let (mut lo, mut hi) = (best.1 - 1.0e-4, best.1 + 1.0e-4);
    for _ in 0..80 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve_at(m1).distance_to(q) < curve_at(m2).distance_to(q) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    curve_at(0.5 * (lo + hi)).distance_to(q)
}

/// One turn of exact samples at the radial arm's density (256 per
/// wrap), the first interval cut short by a seam crossing the way the
/// radial level-curve splitter emits them. The fitted curve must stay
/// on the sampled curve BETWEEN samples, not only at them: the
/// least-squares rung at 4/5 of the sample count meets the residual
/// target at every sample yet rings inside the clamped end spans
/// (PK-673: 1.0e-4 off the flank in the first interval of the nut's
/// crest × flank records); the ladder now interpolates past half the
/// sample count.
#[test]
fn exact_samples_fit_without_end_ringing() {
    let step = std::f64::consts::TAU / 256.0;
    let mut thetas = vec![0.0, 0.45 * step];
    let mut t = 0.45 * step;
    while t + step < std::f64::consts::TAU {
        t += step;
        thetas.push(t);
    }
    thetas.push(std::f64::consts::TAU);
    let samples: Vec<MarchSample> = thetas
        .iter()
        .map(|&th| MarchSample {
            u1: th,
            v1: 0.0,
            u2: th,
            v2: 0.0,
            point: curve_at(th),
        })
        .collect();
    let fit = fit_march_samples(&samples, false).expect("fit");
    let mut params = vec![0.0];
    let mut total = 0.0;
    for w in samples.windows(2) {
        total += w[0].point.distance_to(&w[1].point);
        params.push(total);
    }
    let mut worst = 0.0_f64;
    for i in 0..samples.len() - 1 {
        for f in [0.25, 0.5, 0.75] {
            let tm = (params[i] + f * (params[i + 1] - params[i])) / total;
            let q = fit.curve_3d.point_at(tm);
            let hint = thetas[i] + f * (thetas[i + 1] - thetas[i]);
            worst = worst.max(curve_distance(&q, hint));
        }
    }
    assert!(
        worst < 1.0e-6,
        "fit ({} control points) strays {worst:.3e} from the sampled curve between samples",
        fit.curve_3d.control_points.len()
    );
}
