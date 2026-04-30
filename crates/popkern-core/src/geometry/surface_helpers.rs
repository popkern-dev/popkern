//! Numerical helpers for [`super::surface`].
//!
//! Split out of `surface.rs` to keep that file under the 500-line
//! structural cap.  Two groups:
//!
//! * Numerical utilities reusable by any surface — Gauss-Newton
//!   closest-point projection and grid-based AABB sampling.
//! * Analytic-to-NURBS conversion helpers used by
//!   `SurfaceGeometry::to_nurbs_in_range` (uniform knot vectors,
//!   sub-rectangle cropping, fallback degree-3 sampling).

use crate::geometry::Point3;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Gauss-Newton closest-point projection for any `SurfaceEval`.
///
/// Multi-start from 5 candidates (center + 4 quadrant midpoints) within the
/// parameter range. Returns the (u, v) that minimizes distance to `point`,
/// or `None` if all starts fail to converge.
pub fn closest_point_gauss_newton(surface: &dyn SurfaceEval, point: &Point3) -> Option<(f64, f64)> {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    // Clamp infinite ranges for the search
    let u_lo = u0.max(-100.0);
    let u_hi = u1.min(100.0);
    let v_lo = v0.max(-100.0);
    let v_hi = v1.min(100.0);

    // 5x5 grid of starts plus domain corners for periodic surfaces.
    // The 5x5 grid covers the interior; corners catch periodic boundary
    // cases where the closest point is near u=0 or v=0 on a torus/sphere.
    let n_grid = 5;
    let mut starts = Vec::with_capacity(n_grid * n_grid + 4);
    for i in 0..n_grid {
        let u = u_lo + (u_hi - u_lo) * (i as f64 + 0.5) / (n_grid as f64);
        for j in 0..n_grid {
            let v = v_lo + (v_hi - v_lo) * (j as f64 + 0.5) / (n_grid as f64);
            starts.push((u, v));
        }
    }
    // Add domain corners (important for periodic surfaces)
    starts.push((u_lo, v_lo));
    starts.push((u_lo, v_hi));
    starts.push((u_hi, v_lo));
    starts.push((u_hi, v_hi));

    let max_iter = 20;
    let step_tol = 1e-12;
    let max_dist = 1e6;

    let mut best: Option<(f64, f64, f64)> = None; // (u, v, dist)

    for (u_start, v_start) in starts {
        let mut u = u_start;
        let mut v = v_start;

        for _ in 0..max_iter {
            let p = surface.point_at(u, v);
            let diff = p - *point;
            let dist = diff.length();
            if dist > max_dist {
                break;
            }

            let su = surface.partial_u(u, v);
            let sv = surface.partial_v(u, v);

            // 2x2 normal equation: J^T J [du, dv]^T = -J^T r
            let a00 = su.dot(&su);
            let a01 = su.dot(&sv);
            let a11 = sv.dot(&sv);
            let b0 = diff.dot(&su);
            let b1 = diff.dot(&sv);

            let det = a00 * a11 - a01 * a01;
            if det.abs() < 1e-30 {
                break;
            }

            let du = -(a11 * b0 - a01 * b1) / det;
            let dv = -(a00 * b1 - a01 * b0) / det;

            u += du;
            v += dv;

            if du * du + dv * dv < step_tol {
                let final_p = surface.point_at(u, v);
                let final_dist = final_p.distance_to(point);
                match &best {
                    Some((_, _, bd)) if final_dist >= *bd => {}
                    _ => best = Some((u, v, final_dist)),
                }
                break;
            }
        }

        // Even if we didn't converge, record the best distance found
        let p = surface.point_at(u, v);
        let dist = p.distance_to(point);
        match &best {
            Some((_, _, bd)) if dist >= *bd => {}
            _ => best = Some((u, v, dist)),
        }
    }

    best.map(|(u, v, _)| (u, v))
}

/// Compute AABB by sampling the surface on a grid.
pub fn bounding_box_by_sampling(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Aabb {
    let n = 8;
    let mut points = Vec::with_capacity((n + 1) * (n + 1));

    for i in 0..=n {
        let u = u_range.0 + (u_range.1 - u_range.0) * (i as f64) / (n as f64);
        for j in 0..=n {
            let v = v_range.0 + (v_range.1 - v_range.0) * (j as f64) / (n as f64);
            points.push(surface.point_at(u, v));
        }
    }

    Aabb::from_points(&points).unwrap_or(Aabb {
        min: points
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0)),
        max: points
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0)),
    })
}

// ---------------------------------------------------------------------------
// Analytic-to-NURBS conversion
// ---------------------------------------------------------------------------

/// Build a uniform clamped knot vector for `n_pts` control points and `degree`.
///
/// The knot vector has `degree` repeated knots at each end and uniform interior
/// knots spanning `[range_start, range_end]`.
pub(super) fn uniform_clamped_knots(
    n_pts: usize,
    degree: usize,
    range_start: f64,
    range_end: f64,
) -> Vec<f64> {
    // n_pts control points, degree `degree` => n_pts + degree + 1 knots
    // interior spans: n_pts - degree
    let n_interior = n_pts - degree; // number of distinct interior spans (>= 1)
    let mut knots = Vec::with_capacity(n_pts + degree + 1);
    // `degree + 1` repeated start knots
    for _ in 0..=degree {
        knots.push(range_start);
    }
    // interior knots
    for i in 1..n_interior {
        let t = range_start + (range_end - range_start) * (i as f64) / (n_interior as f64);
        knots.push(t);
    }
    // `degree + 1` repeated end knots
    for _ in 0..=degree {
        knots.push(range_end);
    }
    knots
}

/// Crop a full-range NURBS surface to the `(u_range, v_range)`
/// sub-rectangle by splitting along each boundary and keeping the
/// middle piece. Used by Gregory `to_nurbs_in_range` sub-range
/// handling (§5.3 item #9) so the rational-exact form survives
/// region restriction instead of degrading to `sample_to_nurbs`.
/// Falls back to the input if any split fails (e.g., a range
/// endpoint equals the existing boundary — no slice needed on
/// that side).
pub(super) fn crop_nurbs_surface_to_range(
    surface: NurbsSurface,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    let mut s = surface;
    // u_lo: split at u_range.0, keep right.
    if u_range.0 > 0.0 {
        match s.split_at_u(u_range.0) {
            Ok((_left, right)) => s = right,
            Err(_) => return s,
        }
    }
    // u_hi: split at u_range.1, keep left. Guard: only split if
    // u_range.1 is strictly inside the current surface's u span.
    if u_range.1 < 1.0 {
        match s.split_at_u(u_range.1) {
            Ok((left, _right)) => s = left,
            Err(_) => return s,
        }
    }
    if v_range.0 > 0.0 {
        match s.split_at_v(v_range.0) {
            Ok((_bottom, top)) => s = top,
            Err(_) => return s,
        }
    }
    if v_range.1 < 1.0 {
        match s.split_at_v(v_range.1) {
            Ok((bottom, _top)) => s = bottom,
            Err(_) => return s,
        }
    }
    s
}

/// Sample the surface on an n_u × n_v grid and use the sampled points as NURBS
/// control points (degree-3 interpolating approximation).
pub(super) fn sample_to_nurbs(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    let n_u: usize = 9;
    let n_v: usize = 9;
    let degree: usize = 3;

    // Build parameter sample positions (uniform over range)
    let u_params: Vec<f64> = (0..n_u)
        .map(|i| u_range.0 + (u_range.1 - u_range.0) * (i as f64) / ((n_u - 1) as f64))
        .collect();
    let v_params: Vec<f64> = (0..n_v)
        .map(|j| v_range.0 + (v_range.1 - v_range.0) * (j as f64) / ((n_v - 1) as f64))
        .collect();

    // Evaluate grid: control_points[i_u * n_v + i_v]
    let control_points: Vec<Point3> = u_params
        .iter()
        .flat_map(|&u| v_params.iter().map(move |&v| surface.point_at(u, v)))
        .collect();

    let weights = vec![1.0_f64; n_u * n_v];

    let knots_u = uniform_clamped_knots(n_u, degree, u_range.0, u_range.1);
    let knots_v = uniform_clamped_knots(n_v, degree, v_range.0, v_range.1);

    NurbsSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights,
    }
}
