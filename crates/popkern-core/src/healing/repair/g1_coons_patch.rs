//! G1-continuous Coons patch construction.
//!
//! Extends the bilinear Coons patch with Hermite-blended cross-boundary
//! derivative corrections to achieve tangent-plane (G1) continuity at
//! all four boundary curves.
//!
//! The standard (G0) Coons patch interpolates position only. The G1
//! variant additionally enforces that the surface tangent plane matches
//! prescribed cross-derivatives at each boundary, eliminating visible
//! crease artifacts when the fill surface meets adjacent faces.

use crate::geometry::Vec3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_surface::NurbsSurface;

use super::coons_patch::{
    CoonsBoundary, CoonsPatchResult, build_clamped_knots, coons_eval, measure_boundary_deviation,
    vec_to_point,
};

/// Cross-boundary derivative information for G1 Coons patch.
///
/// Each field holds sampled derivative vectors along one boundary.
/// These specify the surface tangent plane direction perpendicular to
/// the boundary, enabling G1 (tangent-plane) continuity at the patch
/// edges.
#[derive(Debug, Clone)]
pub struct CrossDerivatives {
    /// Cross-derivative along bottom boundary (dS/dv at v=0).
    pub bottom: Vec<Vec3>,
    /// Cross-derivative along top boundary (dS/dv at v=1).
    pub top: Vec<Vec3>,
    /// Cross-derivative along left boundary (dS/du at u=0).
    pub left: Vec<Vec3>,
    /// Cross-derivative along right boundary (dS/du at u=1).
    pub right: Vec<Vec3>,
}

/// Hermite basis function h10(t) = t^3 - 2t^2 + t.
fn hermite_h10(t: f64) -> f64 {
    t * t * t - 2.0 * t * t + t
}

/// Hermite basis function h11(t) = t^3 - t^2.
fn hermite_h11(t: f64) -> f64 {
    t * t * t - t * t
}

/// Linearly interpolate cross-derivative samples at parameter `t` in [0,1].
fn interpolate_cross_derivs(samples: &[Vec3], t: f64) -> Vec3 {
    if samples.is_empty() {
        return Vec3::zero();
    }
    if samples.len() == 1 {
        return samples[0];
    }
    let n = samples.len() - 1;
    let idx_f = t * n as f64;
    let idx = (idx_f.floor() as usize).min(n - 1);
    let frac = idx_f - idx as f64;
    samples[idx] * (1.0 - frac) + samples[idx + 1] * frac
}

/// Construct a G1-continuous Coons patch from boundary curves and cross-derivatives.
///
/// This is the bicubic Coons patch that interpolates both position and
/// tangent plane at the boundaries. The cross-derivatives at each boundary
/// specify the surface normal direction perpendicular to the boundary curve.
///
/// The construction adds a Hermite-blended cross-derivative correction to
/// the standard bilinear Coons patch, ensuring that the resulting surface
/// matches the prescribed tangent planes along all four boundaries.
///
/// Returns `None` if boundary curves don't close, cross-derivative arrays
/// are empty, or the grid dimensions are too small.
pub fn construct_g1_coons_patch(
    boundary: &CoonsBoundary,
    cross_derivs: &CrossDerivatives,
    n_u: usize,
    n_v: usize,
    tolerance: f64,
) -> Option<CoonsPatchResult> {
    // Validate cross-derivative arrays are non-empty.
    if cross_derivs.bottom.is_empty()
        || cross_derivs.top.is_empty()
        || cross_derivs.left.is_empty()
        || cross_derivs.right.is_empty()
    {
        return None;
    }

    if n_u < 4 || n_v < 4 {
        return None;
    }

    // Extract and verify corner closure (same as G0 path).
    let (t0_b, t1_b) = boundary.bottom.parameter_range();
    let (t0_t, t1_t) = boundary.top.parameter_range();
    let (t0_l, t1_l) = boundary.left.parameter_range();
    let (t0_r, t1_r) = boundary.right.parameter_range();

    let p00 = boundary.bottom.point_at(t0_b);
    let p10 = boundary.bottom.point_at(t1_b);
    let p01 = boundary.top.point_at(t0_t);
    let p11 = boundary.top.point_at(t1_t);

    let l_start = boundary.left.point_at(t0_l);
    let l_end = boundary.left.point_at(t1_l);
    let r_start = boundary.right.point_at(t0_r);
    let r_end = boundary.right.point_at(t1_r);

    if p00.distance_to(&l_start) > tolerance
        || p10.distance_to(&r_start) > tolerance
        || p01.distance_to(&l_end) > tolerance
        || p11.distance_to(&r_end) > tolerance
    {
        return None;
    }

    let corners = [p00, p10, p01, p11];

    // Corner twist vectors (cross-derivatives at corners).
    let twist_00 = interpolate_cross_derivs(&cross_derivs.bottom, 0.0);
    let twist_10 = interpolate_cross_derivs(&cross_derivs.bottom, 1.0);
    let twist_01 = interpolate_cross_derivs(&cross_derivs.top, 0.0);
    let twist_11 = interpolate_cross_derivs(&cross_derivs.top, 1.0);

    // Sample the G1 Coons patch on an n_u x n_v grid.
    let mut control_points = Vec::with_capacity(n_u * n_v);
    let mut weights = Vec::with_capacity(n_u * n_v);

    for i_u in 0..n_u {
        let u = i_u as f64 / (n_u - 1) as f64;
        for i_v in 0..n_v {
            let v = i_v as f64 / (n_v - 1) as f64;

            // Step 1: Bilinear Coons patch (G0).
            let g0 = coons_eval(
                &boundary.bottom,
                &boundary.top,
                &boundary.left,
                &boundary.right,
                &corners,
                u,
                v,
            );

            // Step 2: Hermite-blended cross-derivative correction.
            let cd_bottom = interpolate_cross_derivs(&cross_derivs.bottom, u);
            let cd_top = interpolate_cross_derivs(&cross_derivs.top, u);
            let cd_left = interpolate_cross_derivs(&cross_derivs.left, v);
            let cd_right = interpolate_cross_derivs(&cross_derivs.right, v);

            // Cross-derivative correction in v direction.
            let v_correction = cd_bottom * hermite_h10(v) + cd_top * hermite_h11(v);
            // Cross-derivative correction in u direction.
            let u_correction = cd_left * hermite_h10(u) + cd_right * hermite_h11(u);

            // Twist correction at corners (bilinear blend of corner twists).
            let twist_correction = (twist_00 * ((1.0 - u) * hermite_h10(v))
                + twist_10 * (u * hermite_h10(v)))
                + (twist_01 * ((1.0 - u) * hermite_h11(v)) + twist_11 * (u * hermite_h11(v)));

            let pt = g0 + v_correction + u_correction - twist_correction;
            control_points.push(vec_to_point(pt));
            weights.push(1.0);
        }
    }

    // Build the NURBS surface.
    let degree = 3;
    let knots_u = build_clamped_knots(n_u, degree, 0.0, 1.0);
    let knots_v = build_clamped_knots(n_v, degree, 0.0, 1.0);

    let surface = NurbsSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights,
    };

    let max_boundary_deviation = measure_boundary_deviation(&surface, boundary);

    Some(CoonsPatchResult {
        surface,
        max_boundary_deviation,
        corners,
    })
}

#[cfg(test)]
#[path = "g1_coons_patch_tests.rs"]
mod tests;
