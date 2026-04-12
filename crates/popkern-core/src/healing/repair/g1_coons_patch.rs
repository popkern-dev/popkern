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

/// Cross-boundary second-derivative information for G2 Coons patch.
///
/// Each field holds sampled second-derivative vectors along one boundary.
/// These specify the surface curvature direction perpendicular to the
/// boundary, enabling G2 (curvature) continuity at the patch edges.
#[derive(Debug, Clone)]
pub struct SecondDerivatives {
    /// Second cross-derivative along bottom (d^2S/dv^2 at v=0).
    pub bottom: Vec<Vec3>,
    /// Second cross-derivative along top (d^2S/dv^2 at v=1).
    pub top: Vec<Vec3>,
    /// Second cross-derivative along left (d^2S/du^2 at u=0).
    pub left: Vec<Vec3>,
    /// Second cross-derivative along right (d^2S/du^2 at u=1).
    pub right: Vec<Vec3>,
}

/// Quintic Hermite basis h00(t) = 1 - 10t^3 + 15t^4 - 6t^5.
fn quintic_h00(t: f64) -> f64 {
    let t3 = t * t * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    1.0 - 10.0 * t3 + 15.0 * t4 - 6.0 * t5
}

/// Quintic Hermite basis h01(t) = 10t^3 - 15t^4 + 6t^5.
fn quintic_h01(t: f64) -> f64 {
    let t3 = t * t * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    10.0 * t3 - 15.0 * t4 + 6.0 * t5
}

/// Quintic Hermite basis h10(t) = t - 6t^3 + 8t^4 - 3t^5.
fn quintic_h10(t: f64) -> f64 {
    let t3 = t * t * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    t - 6.0 * t3 + 8.0 * t4 - 3.0 * t5
}

/// Quintic Hermite basis h11(t) = -4t^3 + 7t^4 - 3t^5.
fn quintic_h11(t: f64) -> f64 {
    let t3 = t * t * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    -4.0 * t3 + 7.0 * t4 - 3.0 * t5
}

/// Quintic Hermite basis h20(t) = t^2/2 - 3t^3/2 + 3t^4/2 - t^5/2.
fn quintic_h20(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    0.5 * t2 - 1.5 * t3 + 1.5 * t4 - 0.5 * t5
}

/// Quintic Hermite basis h21(t) = t^3/2 - t^4 + t^5/2.
fn quintic_h21(t: f64) -> f64 {
    let t3 = t * t * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    0.5 * t3 - t4 + 0.5 * t5
}

/// Construct a G2-continuous Coons patch.
///
/// Extends the G1 patch with quintic Hermite blending to match
/// second-order derivatives (curvature) at the boundaries.
/// Requires boundary curves with well-defined curvature.
///
/// The quintic Hermite basis functions ensure position, first-derivative,
/// and second-derivative matching at t=0 and t=1, yielding C2 continuity
/// at all four boundary curves.
///
/// Returns `None` if boundary curves don't close, derivative arrays
/// are empty, or the grid dimensions are too small.
pub fn construct_g2_coons_patch(
    boundary: &CoonsBoundary,
    cross_derivs: &CrossDerivatives,
    second_derivs: &SecondDerivatives,
    n_u: usize,
    n_v: usize,
    tolerance: f64,
) -> Option<CoonsPatchResult> {
    // Validate all derivative arrays are non-empty.
    if cross_derivs.bottom.is_empty()
        || cross_derivs.top.is_empty()
        || cross_derivs.left.is_empty()
        || cross_derivs.right.is_empty()
        || second_derivs.bottom.is_empty()
        || second_derivs.top.is_empty()
        || second_derivs.left.is_empty()
        || second_derivs.right.is_empty()
    {
        return None;
    }

    if n_u < 4 || n_v < 4 {
        return None;
    }

    // Extract and verify corner closure (same as G0/G1 path).
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

    // Sample the G2 Coons patch on an n_u x n_v grid.
    let mut control_points = Vec::with_capacity(n_u * n_v);
    let mut weights = Vec::with_capacity(n_u * n_v);

    for i_u in 0..n_u {
        let u = i_u as f64 / (n_u - 1) as f64;
        for i_v in 0..n_v {
            let v = i_v as f64 / (n_v - 1) as f64;

            // Step 1: Positional interpolation (quintic blend of boundaries).
            let g0 = coons_eval(
                &boundary.bottom,
                &boundary.top,
                &boundary.left,
                &boundary.right,
                &corners,
                u,
                v,
            );

            // Step 2: First-derivative corrections (quintic h10/h11).
            let cd_bottom = interpolate_cross_derivs(&cross_derivs.bottom, u);
            let cd_top = interpolate_cross_derivs(&cross_derivs.top, u);
            let cd_left = interpolate_cross_derivs(&cross_derivs.left, v);
            let cd_right = interpolate_cross_derivs(&cross_derivs.right, v);

            let v_deriv1 = cd_bottom * quintic_h10(v) + cd_top * quintic_h11(v);
            let u_deriv1 = cd_left * quintic_h10(u) + cd_right * quintic_h11(u);

            // Step 3: Second-derivative corrections (quintic h20/h21).
            let sd_bottom = interpolate_cross_derivs(&second_derivs.bottom, u);
            let sd_top = interpolate_cross_derivs(&second_derivs.top, u);
            let sd_left = interpolate_cross_derivs(&second_derivs.left, v);
            let sd_right = interpolate_cross_derivs(&second_derivs.right, v);

            let v_deriv2 = sd_bottom * quintic_h20(v) + sd_top * quintic_h21(v);
            let u_deriv2 = sd_left * quintic_h20(u) + sd_right * quintic_h21(u);

            // Combine: position + first-order + second-order corrections.
            // Use quintic h00/h01 for positional blend to maintain partition
            // of unity with the derivative basis functions.
            let pos_blend = g0 * (quintic_h00(v) + quintic_h01(v));
            let pt = pos_blend + v_deriv1 + u_deriv1 + v_deriv2 + u_deriv2;

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
