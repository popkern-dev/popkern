//! Lifting a trim polygon that runs through a chart pole (PK-936).
//!
//! A loop on a periodic chart that passes a pole (a cone's apex, a
//! sphere's pole: u is arbitrary there) can cover the whole period, and
//! where it enters and leaves the pole its u is the curve's own, not the
//! polygon's. A cone wall less a notch through its seam runs the rim
//! from 2π to 0, down the seam to the apex, and back up the seam's other
//! copy, which the folded convention also writes at u = 0: chained, the
//! polygon misses closure by a period; folded about its circular mean
//! it is a self-crossing scribble (measured: the cone minus a box on its
//! axis, wall loop at u ∈ [−3.53, 2.55], the grid meshed 0.614 of the
//! wall's 6.410). The 7/8 sphere face of a cube-corner union runs a
//! seam spine to the north pole and back, then enters and leaves its
//! octant hole along the south pole row: it meshed nothing.
//!
//! In the chart the loop runs along the pole row between the samples
//! either side of the pole, which costs nothing in 3D; the run's length
//! is the loop's choice. As in the boundary-first mesher
//! (`resolve_pole_crossings`), each run takes the short step or a period
//! more or less, and the loop keeps the first choice that closes with
//! the winding the face flag names.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Azimuth gap, as a fraction of the period, below which the samples
/// either side of a pole sit on one meridian (a seam pair over the
/// pole): the run may be a slit or a whole period either way.
const POLE_VERTEX_GAP: f64 = 1e-6;
/// Most pole runs a loop may enumerate (3 choices each).
const MAX_POLE_RUNS: usize = 4;

/// Rewrite `points`/`points_3d` with explicit pole-row runs so the loop
/// closes in u and its signed uv area is positive exactly when `ccw`.
/// Returns `false`, leaving the inputs untouched, when the loop visits
/// no pole or no choice of runs closes that way.
pub(crate) fn lift_through_pole(
    points: &mut Vec<(f64, f64)>,
    points_3d: &mut Vec<Point3>,
    surface: &SurfaceGeometry,
    period_u: f64,
    ccw: bool,
) -> bool {
    let n = points.len();
    if n < 3 || points_3d.len() != n || period_u <= 0.0 {
        return false;
    }
    let pole: Vec<bool> = points
        .iter()
        .map(|&(u, v)| is_pole(surface, u, v, period_u))
        .collect();
    let regular: Vec<usize> = (0..n).filter(|&i| !pole[i]).collect();
    if regular.len() < 2 || regular.len() == n {
        return false;
    }
    let wrap = |d: f64| d - period_u * (d / period_u).round();
    // Consecutive regular pairs (cyclic) with the short step between
    // them and, when pole samples sit between, the pole's sample.
    let pairs: Vec<(usize, usize, f64, Option<usize>)> = (0..regular.len())
        .map(|k| {
            let (i, j) = (regular[k], regular[(k + 1) % regular.len()]);
            let at_pole = (i + 1) % n;
            let run = (at_pole != j).then_some(at_pole);
            (i, j, wrap(points[j].0 - points[i].0), run)
        })
        .collect();
    let runs: Vec<usize> = (0..pairs.len()).filter(|&k| pairs[k].3.is_some()).collect();
    if runs.len() > MAX_POLE_RUNS {
        return false;
    }
    let options = |k: usize| -> Vec<f64> {
        let d = pairs[k].2;
        if d.abs() <= POLE_VERTEX_GAP * period_u {
            vec![0.0, 1.0, -1.0]
        } else {
            vec![0.0, -d.signum()]
        }
    };
    let combos: usize = runs.iter().map(|&k| options(k).len()).product();
    for code in 0..combos {
        let mut turns = vec![0.0; pairs.len()];
        let mut c = code;
        for &k in &runs {
            let opts = options(k);
            turns[k] = opts[c % opts.len()];
            c /= opts.len();
        }
        let (uv, xyz, total) = assemble(points, points_3d, &pairs, &turns, period_u);
        if (total / period_u).round() == 0.0 && (signed_area(&uv) > 0.0) == ccw {
            *points = uv;
            *points_3d = xyz;
            return true;
        }
    }
    false
}

/// The loop walked from its first regular sample: each pair's step
/// plus its turn, with a pole pair's run written as its two pole-row
/// ends. Returns the polygon and the total u travelled.
fn assemble(
    points: &[(f64, f64)],
    points_3d: &[Point3],
    pairs: &[(usize, usize, f64, Option<usize>)],
    turns: &[f64],
    period_u: f64,
) -> (Vec<(f64, f64)>, Vec<Point3>, f64) {
    let (mut uv, mut xyz) = (Vec::new(), Vec::new());
    let mut u = points[pairs[0].0].0;
    let mut total = 0.0;
    for (&(i, _, d, run), &t) in pairs.iter().zip(turns) {
        uv.push((u, points[i].1));
        xyz.push(points_3d[i]);
        let step = d + t * period_u;
        if let Some(p) = run {
            let v_pole = points[p].1;
            uv.push((u, v_pole));
            uv.push((u + step, v_pole));
            xyz.push(points_3d[p]);
            xyz.push(points_3d[p]);
        }
        u += step;
        total += step;
    }
    (uv, xyz, total)
}

/// `u_range` widened, on a periodic u axis, to cover the outer trim:
/// the pcurve box knows nothing of the branch a pole-row lift moves the
/// loop to (the 7/8 sphere face's grid spanned u 3.80..5.72, its lifted
/// trim −2.44..0.70). A range that already covers it is returned as is.
pub(crate) fn cover_outer_trim(
    u_range: (f64, f64),
    trim_loops: &[crate::trim::TrimLoop],
    surface: &SurfaceGeometry,
) -> (f64, f64) {
    let Some(p) = surface.periodicity().0.filter(|p| *p > 0.0) else {
        return u_range;
    };
    let Some(outer) = trim_loops.iter().find(|l| l.is_outer) else {
        return u_range;
    };
    let lo = outer
        .points
        .iter()
        .map(|q| q.0)
        .fold(f64::INFINITY, f64::min);
    let hi = outer
        .points
        .iter()
        .map(|q| q.0)
        .fold(f64::NEG_INFINITY, f64::max);
    let slack = 1e-9 * p;
    if lo < u_range.0 - slack || hi > u_range.1 + slack {
        (lo.min(u_range.0), hi.max(u_range.1))
    } else {
        u_range
    }
}

/// Signed shoelace area of a closed uv polygon.
pub(crate) fn signed_area(points: &[(f64, f64)]) -> f64 {
    let n = points.len();
    (0..n)
        .map(|i| {
            let (a, b) = (points[i], points[(i + 1) % n]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        / 2.0
}

/// Whether `(u, v)` lies on a collapsed u-row: the surface point does
/// not move as u turns.
pub(crate) fn is_pole(surface: &SurfaceGeometry, u: f64, v: f64, period_u: f64) -> bool {
    let p = surface.point_at(u, v);
    let scale = 1.0 + (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    [1.0 / 3.0, 2.0 / 3.0].iter().all(|f| {
        let q = surface.point_at(u + f * period_u, v);
        (q - p).length() < 1e-9 * scale
    })
}

#[cfg(test)]
#[path = "trim_pole_lift_tests.rs"]
mod tests;
