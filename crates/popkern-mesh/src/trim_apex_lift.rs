//! Lifting a trim polygon that runs through an apex (PK-936).
//!
//! A loop on a cone's periodic chart that passes its apex (a geometric
//! singularity: u is arbitrary there) can cover the whole period, and
//! where it enters and leaves the apex its u is the curve's own, not the
//! polygon's. A cone wall less a notch through its seam runs the rim
//! from 2π to 0, down the seam to the apex, and back up the seam's other
//! copy, which the folded convention also writes at u = 0: chained, the
//! polygon misses closure by a period; folded about its circular mean
//! it is a self-crossing scribble (measured: the cone minus a box on its
//! axis, wall loop at u ∈ [−3.53, 2.55], the grid meshed 0.614 of the
//! wall's 6.410).
//!
//! In the chart the loop runs along the apex row between the samples
//! either side of the apex, which costs nothing in 3D; the run's length
//! is the loop's choice: each run takes the short step or a period more
//! or less, and the loop keeps the first choice that closes with the
//! winding the face flag names. Which samples lie at the apex is the
//! chart service's answer (the caller's `at_apex`). A chart pole never
//! reaches here: a face reaching one meshes on the atlas, whose pieces
//! hold no singular row, on every rung.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;

/// Azimuth gap, as a fraction of the period, below which the samples
/// either side of the apex sit on one meridian (a seam pair over the
/// apex): the run may be a slit or a whole period either way.
const APEX_VERTEX_GAP: f64 = 1e-6;
/// Most apex runs a loop may enumerate (3 choices each).
const MAX_APEX_RUNS: usize = 4;

/// Rewrite `points`/`points_3d` with explicit apex-row runs so the loop
/// closes in u and its signed uv area is positive exactly when `ccw`.
/// `at_apex` says which positions lie at an apex. Returns `false`,
/// leaving the inputs untouched, when the loop visits no apex or no
/// choice of runs closes that way.
pub(crate) fn lift_through_apex(
    points: &mut Vec<(f64, f64)>,
    points_3d: &mut Vec<Point3>,
    at_apex: &dyn Fn(&Point3) -> bool,
    period_u: f64,
    ccw: bool,
) -> bool {
    let n = points.len();
    if n < 3 || points_3d.len() != n || period_u <= 0.0 {
        return false;
    }
    let apex: Vec<bool> = points_3d.iter().map(at_apex).collect();
    let regular: Vec<usize> = (0..n).filter(|&i| !apex[i]).collect();
    if regular.len() < 2 || regular.len() == n {
        return false;
    }
    let wrap = |d: f64| d - period_u * (d / period_u).round();
    // Consecutive regular pairs (cyclic) with the short step between
    // them and, when apex samples sit between, the apex's sample.
    let pairs: Vec<(usize, usize, f64, Option<usize>)> = (0..regular.len())
        .map(|k| {
            let (i, j) = (regular[k], regular[(k + 1) % regular.len()]);
            let at_apex_sample = (i + 1) % n;
            let run = (at_apex_sample != j).then_some(at_apex_sample);
            (i, j, wrap(points[j].0 - points[i].0), run)
        })
        .collect();
    let runs: Vec<usize> = (0..pairs.len()).filter(|&k| pairs[k].3.is_some()).collect();
    if runs.len() > MAX_APEX_RUNS {
        return false;
    }
    let options = |k: usize| -> Vec<f64> {
        let d = pairs[k].2;
        if d.abs() <= APEX_VERTEX_GAP * period_u {
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
/// plus its turn, with an apex pair's run written as its two apex-row
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
            let v_apex = points[p].1;
            uv.push((u, v_apex));
            uv.push((u + step, v_apex));
            xyz.push(points_3d[p]);
            xyz.push(points_3d[p]);
        }
        u += step;
        total += step;
    }
    (uv, xyz, total)
}

/// `u_range` widened, on a periodic u axis, to cover the outer trim:
/// the pcurve box knows nothing of the branch an apex-row lift moves the
/// loop to (the 7/8 sphere face's grid spanned u 3.80..5.72, its lifted
/// trim −2.44..0.70). A range that already covers it is returned as is.
pub(crate) fn cover_outer_trim(
    u_range: (f64, f64),
    trim_loops: &[crate::trim::TrimLoop],
    surface: &SurfaceGeometry,
) -> (f64, f64) {
    let Some(p) = popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(surface)
        .0
        .filter(|p| *p > 0.0)
    else {
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

#[cfg(test)]
#[path = "trim_apex_lift_tests.rs"]
mod tests;
