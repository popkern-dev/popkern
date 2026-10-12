//! Overlap detection for surface-curve intersections.
//!
//! Detects segments where a curve lies on a surface by sampling the curve
//! and projecting each sample point onto the surface via a Newton closest-point
//! iteration.  Contiguous runs of 3 or more samples within linear tolerance are
//! reported as `SurfaceCurveIntersection::Overlap`.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Number of evenly-spaced samples along the curve used for overlap detection.
const NUM_SAMPLES: usize = 30;

/// Maximum Newton iterations for the closest-point projection.
const MAX_ITER: usize = 10;

/// Minimum consecutive on-surface samples needed to declare an overlap.
const MIN_RUN: usize = 3;

// ---------------------------------------------------------------------------
// Closest-point projection
// ---------------------------------------------------------------------------

/// Run a single Gauss-Newton closest-point projection from the given start.
///
/// Minimises `|P - S(u,v)|²` via Newton-Raphson on the gradient equations:
///
///   ∂/∂u = -2·(P - S)·S_u = 0
///   ∂/∂v = -2·(P - S)·S_v = 0
///
/// Uses the Gauss-Newton approximation (second-order terms of S dropped).
/// Returns the squared distance at the converged (or final) iterate.
fn project_from(
    surface: &SurfaceGeometry,
    p_x: f64,
    p_y: f64,
    p_z: f64,
    u_init: f64,
    v_init: f64,
) -> f64 {
    use crate::geometry::Vec3;

    let mut u = u_init;
    let mut v = v_init;

    for _ in 0..MAX_ITER {
        let s = surface.point_at(u, v);
        let diff = Vec3::new(p_x - s.x, p_y - s.y, p_z - s.z);

        let su = surface.partial_u(u, v);
        let sv = surface.partial_v(u, v);

        let g_u = -diff.dot(&su);
        let g_v = -diff.dot(&sv);

        if g_u * g_u + g_v * g_v < 1e-20 {
            break;
        }

        let j00 = su.dot(&su);
        let j01 = su.dot(&sv);
        let j11 = sv.dot(&sv);
        let det = j00 * j11 - j01 * j01;

        if det.abs() < 1e-30 {
            break;
        }

        let inv = 1.0 / det;
        u += (-g_u * j11 + g_v * j01) * inv;
        v += (g_u * j01 - g_v * j00) * inv;
    }

    let s = surface.point_at(u, v);
    let dx = p_x - s.x;
    let dy = p_y - s.y;
    let dz = p_z - s.z;
    dx * dx + dy * dy + dz * dz
}

/// Project a 3-D point `p` onto `surface`, returning the minimum squared
/// distance found.
///
/// To avoid getting stuck in local minima (e.g. the antipodal closest-point
/// on a closed surface such as a cylinder), we try a small coarse grid of
/// starting points spread over the supplied domain, run [`project_from`] from
/// each, and keep the smallest residual.  The cost is modest because the grid
/// is small (3×3 = 9 trials) and the inner Newton loop is at most `MAX_ITER`
/// (10) iterations.
fn closest_point_sq_dist(
    surface: &SurfaceGeometry,
    p_x: f64,
    p_y: f64,
    p_z: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> f64 {
    const GRID: usize = 3; // seeds per axis → 9 starts total
    let (u0, u1) = u_range;
    let (v0, v1) = v_range;

    let mut best = f64::INFINITY;

    for iu in 0..GRID {
        for iv in 0..GRID {
            let fu = if GRID > 1 {
                iu as f64 / (GRID - 1) as f64
            } else {
                0.5
            };
            let fv = if GRID > 1 {
                iv as f64 / (GRID - 1) as f64
            } else {
                0.5
            };
            let u_start = u0 + fu * (u1 - u0);
            let v_start = v0 + fv * (v1 - v0);
            let sq = project_from(surface, p_x, p_y, p_z, u_start, v_start);
            if sq < best {
                best = sq;
            }
        }
    }

    best
}

/// Box around a positive-weight NURBS surface's control points, widened
/// by `pad`: the surface lies in the control points' convex hull, so a
/// point outside the box is farther than `pad` from every surface
/// point and its closest-point solve can be skipped (most samples of
/// a classification ray against a face it crosses lie outside; the
/// solves were a fifth of a lone screw's build, PK-799). `None` (no
/// rejection) for analytic surfaces and non-positive weights.
fn control_hull_box(surface: &SurfaceGeometry, pad: f64) -> Option<(Point3, Point3)> {
    let SurfaceGeometry::Nurbs(ns) = surface else {
        return None;
    };
    if ns.control_points.is_empty() || ns.weights.iter().any(|&w| w <= 0.0) {
        return None;
    }
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in &ns.control_points {
        lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    Some((
        Point3::new(lo.x - pad, lo.y - pad, lo.z - pad),
        Point3::new(hi.x + pad, hi.y + pad, hi.z + pad),
    ))
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Detect overlapping segments where the curve lies on the surface.
///
/// Samples the curve at [`NUM_SAMPLES`] evenly-spaced parameter values within
/// `t_range`.  For each sample the closest point on the surface is found (see
/// [`closest_point_sq_dist`]).  Contiguous runs of [`MIN_RUN`] or more samples
/// whose distance to the surface is below `tolerance.linear()` are returned as
/// [`SurfaceCurveIntersection::Overlap`] results, one per run.
pub(crate) fn detect_surface_curve_overlaps(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<SurfaceCurveIntersection> {
    let (t0, t1) = t_range;

    let lin_sq = tolerance.linear() * tolerance.linear();

    // --- Step 1: evaluate the distance from each sample to the surface -------

    let mut on_surface = [false; NUM_SAMPLES];
    let hull = control_hull_box(surface, tolerance.linear());
    let on = |i: usize| {
        let frac = if NUM_SAMPLES > 1 {
            i as f64 / (NUM_SAMPLES - 1) as f64
        } else {
            0.0
        };
        let p = curve.point_at(t0 + frac * (t1 - t0));
        if hull.is_some_and(|(lo, hi)| {
            p.x < lo.x || p.y < lo.y || p.z < lo.z || p.x > hi.x || p.y > hi.y || p.z > hi.z
        }) {
            return false;
        }
        closest_point_sq_dist(surface, p.x, p.y, p.z, u_range, v_range) < lin_sq
    };
    // Every run of MIN_RUN consecutive samples contains an index that is
    // a multiple of MIN_RUN, so when none of those is on the surface
    // there is no overlap: the remaining closest-point solves are
    // skipped, with the same result. That is the common case (an edge,
    // or a classification ray, against a face it crosses), and the
    // solves were a fifth of a thread cut's intersection stage (PK-762).
    let mut any_stride = false;
    for i in (0..NUM_SAMPLES).step_by(MIN_RUN) {
        on_surface[i] = on(i);
        any_stride |= on_surface[i];
    }
    if !any_stride {
        return Vec::new();
    }
    for (i, slot) in on_surface.iter_mut().enumerate() {
        if i % MIN_RUN != 0 {
            *slot = on(i);
        }
    }

    // --- Step 2: find contiguous runs of MIN_RUN+ on-surface samples ---------

    let mut results = Vec::new();
    let mut i = 0;

    while i < NUM_SAMPLES {
        if on_surface[i] {
            // Start of a potential run.
            let run_start = i;
            while i < NUM_SAMPLES && on_surface[i] {
                i += 1;
            }
            let run_end = i - 1; // inclusive

            let run_len = run_end - run_start + 1;
            if run_len >= MIN_RUN {
                // Map back to curve parameter values.
                let t_first = t0 + (run_start as f64 / (NUM_SAMPLES - 1) as f64) * (t1 - t0);
                let t_last = t0 + (run_end as f64 / (NUM_SAMPLES - 1) as f64) * (t1 - t0);

                results.push(SurfaceCurveIntersection::Overlap {
                    t_range: (t_first, t_last),
                });
            }
        } else {
            i += 1;
        }
    }

    results
}

#[cfg(test)]
#[path = "surface_curve_overlap_tests.rs"]
mod tests;
