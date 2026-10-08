//! Singular crossings of the cylinder×cylinder ruling scan.
//!
//! Along the smaller cylinder's ruling at angle θ the other cylinder
//! is met at the roots of a quadratic; its discriminant `D(θ)` decides
//! the branch structure. A sign change of `D` is a smooth turnaround
//! (the curve's θ-extreme, where the low and high roots join). A
//! *double* zero — `D` touching zero at a local minimum without
//! changing sign — is different: the two root branches CROSS there.
//! It is the singular point of a tangential contact: equal-radius
//! cylinders with intersecting axes (the Steinmetz tee, whose
//! intersection is two ellipses crossing at two points), or a small
//! cylinder internally tangent to a larger one.
//!
//! Treating a touch as a turnaround joins the branches into one
//! self-touching loop that the clip stage cannot classify (measured,
//! PK-727: an equal-radius tee served its operands unmerged, volume =
//! the sum of both cylinders). The exact structure is OPEN branches
//! between consecutive crossings — low root and high root over each
//! arc of θ — all meeting at the crossing points, which become
//! valence-4 vertices of the intersection graph.
//!
//! A minimum counts as a crossing when the two roots there are closer
//! than the linear tolerance (`|D| < (a·tol)²`, since the root
//! separation is `√D / a`), or `D` is within its own rounding floor
//! at the minimum: a near-miss within tolerance is the same
//! point to every downstream consumer, and splitting there is the
//! only spelling that stays manifold.

use std::f64::consts::TAU;

/// The ruling-scan quadratic at θ: `(D, b, c)`.
pub(super) type Quad<'a> = &'a dyn Fn(f64) -> (f64, f64, f64);

/// Golden-section refinement of a minimum of `D` within `[lo, hi]`.
fn refine_minimum(quad: Quad<'_>, mut lo: f64, mut hi: f64) -> f64 {
    let g = 0.5 * (5.0_f64.sqrt() - 1.0);
    let mut x1 = hi - g * (hi - lo);
    let mut x2 = lo + g * (hi - lo);
    let (mut f1, mut f2) = (quad(x1).0, quad(x2).0);
    for _ in 0..80 {
        if f1 <= f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - g * (hi - lo);
            f1 = quad(x1).0;
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + g * (hi - lo);
            f2 = quad(x2).0;
        }
    }
    0.5 * (lo + hi)
}

/// The crossing angles in `[0, 2π)`, ascending: local minima of the
/// sampled discriminant whose refined value lies within the
/// root-separation band.
pub(super) fn singular_crossings(
    quad: Quad<'_>,
    discs: &[f64],
    a_coef: f64,
    linear: f64,
) -> Vec<f64> {
    let n = discs.len();
    let step = TAU / n as f64;
    let band = (a_coef * linear) * (a_coef * linear);
    let mut out: Vec<f64> = Vec::new();
    for i in 0..n {
        let prev = discs[(i + n - 1) % n];
        let next = discs[(i + 1) % n];
        if !(discs[i] <= prev && discs[i] <= next) {
            continue;
        }
        let theta = i as f64 * step;
        let star = refine_minimum(quad, theta - step, theta + step).rem_euclid(TAU);
        let (d, b, c) = quad(star);
        // `D = b² − 4ac` cancels at a double root: its rounding floor
        // is a few ulps of the terms, which can exceed the tolerance
        // band itself (an r = 40 tee: terms ~6e3, floor ~1e-12·6e3).
        let floor = 64.0 * f64::EPSILON * (b * b + 4.0 * (a_coef * c).abs());
        if d.abs() >= band.max(floor) {
            continue;
        }
        // A flat minimum spans two samples: keep one crossing.
        let dup = out.iter().any(|&s| {
            let d = (s - star).abs();
            d.min(TAU - d) < 2.0 * step
        });
        if !dup {
            out.push(star);
        }
    }
    out.sort_by(f64::total_cmp);
    out
}

/// The open branches between consecutive crossings `crossings` as
/// `(θ, t)` sequences, or `None` when the discriminant is negative
/// anywhere away from them (a crossing beside a turnaround: not this
/// regime). Each branch starts and ends exactly at the crossing
/// points (`t = −b / 2a` there, the double root).
pub(super) fn crossing_branches(
    quad: Quad<'_>,
    discs: &[f64],
    crossings: &[f64],
    a_coef: f64,
    samples_per_turn: usize,
) -> Option<Vec<Vec<(f64, f64)>>> {
    let n = discs.len();
    let step = TAU / n as f64;
    let near = |theta: f64| {
        crossings.iter().any(|&s| {
            let d = (s - theta).abs();
            d.min(TAU - d) <= 2.0 * step
        })
    };
    if (0..n).any(|i| discs[i] <= 0.0 && !near(i as f64 * step)) {
        return None;
    }
    let root = |theta: f64, low: bool| -> f64 {
        let (disc, b, _) = quad(theta);
        let s = disc.max(0.0).sqrt();
        if low {
            (-b - s) / (2.0 * a_coef)
        } else {
            (-b + s) / (2.0 * a_coef)
        }
    };
    let double_root = |theta: f64| -> f64 { -quad(theta).1 / (2.0 * a_coef) };
    let mut branches = Vec::new();
    for (j, &start) in crossings.iter().enumerate() {
        let end = match crossings.get(j + 1) {
            Some(&e) => e,
            None => crossings[0] + TAU,
        };
        let span = end - start;
        let m = ((span / TAU * samples_per_turn as f64).ceil() as usize).max(16);
        for low in [true, false] {
            let mut pts = Vec::with_capacity(m + 1);
            pts.push((start, double_root(start)));
            for i in 1..m {
                let theta = start + span * i as f64 / m as f64;
                pts.push((theta, root(theta, low)));
            }
            pts.push((end, double_root(end)));
            branches.push(pts);
        }
    }
    Some(branches)
}
