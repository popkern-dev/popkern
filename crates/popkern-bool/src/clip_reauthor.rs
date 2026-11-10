//! Re-authoring a clipped piece through samples of itself with its
//! ends pinned, at the density the piece needs.
//!
//! The clip's endpoint refinement and its far-junction retrim move a
//! piece's end by interpolating samples of the piece with the end
//! sample replaced: the correction is confined to the last segment,
//! where a bare end-control-point move would shift the whole final
//! span (PK-464). Both once took a fixed 17 samples, which re-authors
//! a dense intersection fit as 16 cubic spans: a quarter of Viviani's
//! curve, its record ending at the sphere's pole, came back 1.5e-3 off
//! both surfaces between samples and the boolean declined (PK-1139).
//! The sample count now doubles from 17 until the refit reproduces the
//! piece between its samples within the budget.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::fitting::fit_curve;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

/// The first sample count tried (the old fixed count).
const FIRST_SAMPLES: usize = 17;
/// The most samples a re-authored piece takes.
const MAX_SAMPLES: usize = 1025;

/// A cubic interpolating the run `run_at` returns over `range`, its
/// samples pinned by `pin` (which replaces end samples), with as few
/// samples as reproduce the run: fitted unpinned, the refit lies
/// within `budget` of the run between every two consecutive samples,
/// distances taken through `to_space` (identity for a 3D curve, the
/// surface for a chart curve). The density is decided unpinned
/// because a pin's correction rings into the segments next to it by
/// a fixed fraction of the move whatever the density: it is the
/// intended bend, not a failure to reproduce. `run_at` maps
/// parameters to points in one call so a chart run can be unwrapped
/// as a whole. Past `MAX_SAMPLES` the densest fit is returned; `None`
/// when no fit can be made.
pub(crate) fn reauthor_pinned_run(
    range: (f64, f64),
    run_at: &dyn Fn(&[f64]) -> Vec<Point3>,
    pin: &dyn Fn(&mut [Point3]),
    to_space: &dyn Fn(&Point3) -> Point3,
    budget: f64,
) -> Option<NurbsCurve> {
    let (t0, t1) = range;
    let mut n = FIRST_SAMPLES;
    loop {
        // The run at twice the density: even entries are fitted, odd
        // ones (the midpoints) check the fit.
        let fine: Vec<f64> = (0..2 * n - 1)
            .map(|k| t0 + (t1 - t0) * k as f64 / (2 * n - 2) as f64)
            .collect();
        let run = run_at(&fine);
        if run.len() != fine.len() {
            return None;
        }
        let params: Vec<f64> = fine.iter().step_by(2).copied().collect();
        let mut pts: Vec<Point3> = run.iter().step_by(2).copied().collect();
        let dense_enough = n >= MAX_SAMPLES || {
            let unpinned = fit_curve(&pts, &params, 3, n).ok()?;
            (0..n - 1).all(|j| {
                let k = 2 * j + 1;
                to_space(&unpinned.point_at(fine[k])).distance_to(&to_space(&run[k])) <= budget
            })
        };
        if dense_enough {
            pin(&mut pts);
            return fit_curve(&pts, &params, 3, n).ok();
        }
        n = 2 * n - 1;
    }
}

#[cfg(test)]
#[path = "clip_reauthor_tests.rs"]
mod tests;
