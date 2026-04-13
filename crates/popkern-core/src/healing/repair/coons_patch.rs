//! Coons patch construction for hole filling.
//!
//! Constructs a surface interpolating four boundary curves using the
//! bilinear Coons patch formula. The result is a NURBS surface that
//! achieves G0 continuity with all boundary curves.
//!
//! The G1-continuous variant (in [`super::g1_coons_patch`]) additionally
//! matches tangent planes at the boundaries using Hermite blending of
//! cross-boundary derivatives.
//!
//! # Coons Patch Theory
//!
//! Given four boundary curves forming a 4-sided hole:
//! - C_0(u): bottom boundary (v=0)
//! - C_1(u): top boundary (v=1)
//! - D_0(v): left boundary (u=0)
//! - D_1(v): right boundary (u=1)
//!
//! The bilinear Coons patch is:
//! ```text
//! S(u,v) = (1-v)*C_0(u) + v*C_1(u)
//!        + (1-u)*D_0(v) + u*D_1(v)
//!        - [(1-u)*(1-v)*P_00 + u*(1-v)*P_10 + (1-u)*v*P_01 + u*v*P_11]
//! ```
//!
//! Where P_00, P_10, P_01, P_11 are the four corner points. This gives
//! G0 continuity with all four boundary curves by construction.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};
use crate::topology::BRep;

use super::hole_fill::DetectedHole;

/// Input boundary curves for a Coons patch.
///
/// The four curves must form a closed loop: C_0(0) = D_0(0),
/// C_0(1) = D_1(0), C_1(0) = D_0(1), C_1(1) = D_1(1).
#[derive(Debug, Clone)]
pub struct CoonsBoundary {
    /// Bottom boundary curve C_0(u), from u=0 to u=1.
    pub bottom: NurbsCurve,
    /// Top boundary curve C_1(u), from u=0 to u=1.
    pub top: NurbsCurve,
    /// Left boundary curve D_0(v), from v=0 to v=1.
    pub left: NurbsCurve,
    /// Right boundary curve D_1(v), from v=0 to v=1.
    pub right: NurbsCurve,
}

/// Result of Coons patch construction.
#[derive(Debug, Clone)]
pub struct CoonsPatchResult {
    /// The constructed fill surface.
    pub surface: NurbsSurface,
    /// Maximum deviation from the boundary curves.
    pub max_boundary_deviation: f64,
    /// Corner point positions used: [P_00, P_10, P_01, P_11].
    pub corners: [Point3; 4],
}

/// Origin used for point-to-vector conversion.
const ORIGIN: Point3 = Point3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};

/// Convert a Point3 to a Vec3 (position vector from origin).
pub(super) fn point_to_vec(p: Point3) -> Vec3 {
    p - ORIGIN
}

/// Convert a Vec3 back to a Point3 (position from origin).
pub(super) fn vec_to_point(v: Vec3) -> Point3 {
    ORIGIN + v
}

/// Evaluate the bilinear Coons patch formula at (u, v).
///
/// Returns the surface point as a Vec3 (position vector from origin).
pub(super) fn coons_eval(
    bottom: &NurbsCurve,
    top: &NurbsCurve,
    left: &NurbsCurve,
    right: &NurbsCurve,
    corners: &[Point3; 4],
    u: f64,
    v: f64,
) -> Vec3 {
    let (t0_b, t1_b) = bottom.parameter_range();
    let (t0_t, t1_t) = top.parameter_range();
    let (t0_l, t1_l) = left.parameter_range();
    let (t0_r, t1_r) = right.parameter_range();

    let t_bottom = t0_b + u * (t1_b - t0_b);
    let t_top = t0_t + u * (t1_t - t0_t);
    let t_left = t0_l + v * (t1_l - t0_l);
    let t_right = t0_r + v * (t1_r - t0_r);

    let c0 = point_to_vec(bottom.point_at(t_bottom));
    let c1 = point_to_vec(top.point_at(t_top));
    let d0 = point_to_vec(left.point_at(t_left));
    let d1 = point_to_vec(right.point_at(t_right));

    let p00 = point_to_vec(corners[0]);
    let p10 = point_to_vec(corners[1]);
    let p01 = point_to_vec(corners[2]);
    let p11 = point_to_vec(corners[3]);

    // S(u,v) = (1-v)*C_0(u) + v*C_1(u)
    //        + (1-u)*D_0(v) + u*D_1(v)
    //        - [(1-u)*(1-v)*P_00 + u*(1-v)*P_10 + (1-u)*v*P_01 + u*v*P_11]
    let ruled_uv = c0 * (1.0 - v) + c1 * v;
    let ruled_vu = d0 * (1.0 - u) + d1 * u;
    let bilinear = p00 * ((1.0 - u) * (1.0 - v))
        + p10 * (u * (1.0 - v))
        + p01 * ((1.0 - u) * v)
        + p11 * (u * v);

    ruled_uv + ruled_vu - bilinear
}

/// Build a clamped knot vector for `n` control points at the given degree.
///
/// The knot vector spans `[range_start, range_end]` with uniform interior
/// spacing and multiplicity `degree + 1` at both ends.
pub(super) fn build_clamped_knots(
    n: usize,
    degree: usize,
    range_start: f64,
    range_end: f64,
) -> Vec<f64> {
    let num_knots = n + degree + 1;
    let mut knots = Vec::with_capacity(num_knots);

    // Interior knot count
    let interior = if n > degree + 1 { n - degree - 1 } else { 0 };

    // Leading clamped knots
    for _ in 0..=degree {
        knots.push(range_start);
    }

    // Uniform interior knots
    for i in 1..=interior {
        let t = range_start + (range_end - range_start) * (i as f64) / ((interior + 1) as f64);
        knots.push(t);
    }

    // Trailing clamped knots
    for _ in 0..=degree {
        knots.push(range_end);
    }

    knots
}

/// Construct a bilinear Coons patch from four boundary curves.
///
/// The resulting NURBS surface interpolates all four boundary curves
/// at sampled points, achieving G0 continuity by construction.
/// The surface is sampled on an `n_u x n_v` grid and the sampled
/// points are used as control points for a degree-3 NURBS surface.
///
/// Returns `None` if the boundary curves don't form a valid closed loop
/// (corner points don't match within tolerance), or if the grid
/// dimensions are too small (minimum 4 in each direction for degree 3).
pub fn construct_coons_patch(
    boundary: &CoonsBoundary,
    n_u: usize,
    n_v: usize,
    tolerance: f64,
) -> Option<CoonsPatchResult> {
    // Degree-3 surface needs at least degree+1 = 4 control points per direction.
    if n_u < 4 || n_v < 4 {
        return None;
    }

    // Step 1: Extract corner points
    let (t0_b, t1_b) = boundary.bottom.parameter_range();
    let (t0_t, t1_t) = boundary.top.parameter_range();
    let (t0_l, t1_l) = boundary.left.parameter_range();
    let (t0_r, t1_r) = boundary.right.parameter_range();

    let p00 = boundary.bottom.point_at(t0_b);
    let p10 = boundary.bottom.point_at(t1_b);
    let p01 = boundary.top.point_at(t0_t);
    let p11 = boundary.top.point_at(t1_t);

    // Step 2: Verify closure -- corners must match within tolerance
    let l_start = boundary.left.point_at(t0_l);
    let l_end = boundary.left.point_at(t1_l);
    let r_start = boundary.right.point_at(t0_r);
    let r_end = boundary.right.point_at(t1_r);

    // Adaptive corner-match check: rejects the patch if any of the
    // 4 corners is farther than `tolerance` from its expected
    // boundary endpoint. Routed through adaptive_dist_sq_lt
    // (negated since we want > tolerance) so the rejection is
    // deterministic at the boundary.
    let tol_sq = tolerance * tolerance;
    let corner_check = |a: &crate::geometry::Point3, b: &crate::geometry::Point3| -> bool {
        crate::adaptive_precision::adaptive_dist_sq_lt((a.x, a.y, a.z), (b.x, b.y, b.z), tol_sq)
    };
    if !corner_check(&p00, &l_start) {
        return None;
    }
    if !corner_check(&p10, &r_start) {
        return None;
    }
    if !corner_check(&p01, &l_end) {
        return None;
    }
    if !corner_check(&p11, &r_end) {
        return None;
    }

    let corners = [p00, p10, p01, p11];

    // Step 3: Sample the Coons patch on an n_u x n_v grid
    let mut control_points = Vec::with_capacity(n_u * n_v);
    let mut weights = Vec::with_capacity(n_u * n_v);

    for i_u in 0..n_u {
        let u = i_u as f64 / (n_u - 1) as f64;
        for i_v in 0..n_v {
            let v = i_v as f64 / (n_v - 1) as f64;
            let pt_vec = coons_eval(
                &boundary.bottom,
                &boundary.top,
                &boundary.left,
                &boundary.right,
                &corners,
                u,
                v,
            );
            control_points.push(vec_to_point(pt_vec));
            weights.push(1.0);
        }
    }

    // Step 4: Build the NURBS surface
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

    // Step 5: Measure boundary deviation
    let max_boundary_deviation = measure_boundary_deviation(&surface, boundary);

    Some(CoonsPatchResult {
        surface,
        max_boundary_deviation,
        corners,
    })
}

/// Measure the maximum deviation of the NURBS surface from the boundary curves.
///
/// Samples each boundary (u=0, u=1, v=0, v=1) at 20 points and computes
/// the maximum distance to the corresponding original boundary curve.
pub(super) fn measure_boundary_deviation(surface: &NurbsSurface, boundary: &CoonsBoundary) -> f64 {
    use crate::geometry::surface::SurfaceEval;

    let n_samples = 20;
    let mut max_dev = 0.0_f64;

    let (t0_b, t1_b) = boundary.bottom.parameter_range();
    let (t0_t, t1_t) = boundary.top.parameter_range();
    let (t0_l, t1_l) = boundary.left.parameter_range();
    let (t0_r, t1_r) = boundary.right.parameter_range();

    let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();

    for i in 0..=n_samples {
        let frac = i as f64 / n_samples as f64;

        // v=0 boundary (bottom)
        let u_param = u_lo + frac * (u_hi - u_lo);
        let surf_pt = surface.point_at(u_param, v_lo);
        let curve_pt = boundary.bottom.point_at(t0_b + frac * (t1_b - t0_b));
        max_dev = max_dev.max(surf_pt.distance_to(&curve_pt));

        // v=1 boundary (top)
        let surf_pt = surface.point_at(u_param, v_hi);
        let curve_pt = boundary.top.point_at(t0_t + frac * (t1_t - t0_t));
        max_dev = max_dev.max(surf_pt.distance_to(&curve_pt));

        // u=0 boundary (left)
        let v_param = v_lo + frac * (v_hi - v_lo);
        let surf_pt = surface.point_at(u_lo, v_param);
        let curve_pt = boundary.left.point_at(t0_l + frac * (t1_l - t0_l));
        max_dev = max_dev.max(surf_pt.distance_to(&curve_pt));

        // u=1 boundary (right)
        let surf_pt = surface.point_at(u_hi, v_param);
        let curve_pt = boundary.right.point_at(t0_r + frac * (t1_r - t0_r));
        max_dev = max_dev.max(surf_pt.distance_to(&curve_pt));
    }

    max_dev
}

/// Construct a Coons patch for a detected 4-sided hole.
///
/// Given a detected 4-sided hole (from [`super::hole_fill::detect_holes`]),
/// extracts the boundary curves, constructs a Coons patch, and returns the
/// fill surface. The BRep topology is not modified -- adding a face requires
/// Euler operators and is left to the caller.
///
/// Returns `None` if the hole is not 4-sided, any edge lacks curve geometry,
/// or the boundary curves fail the closure check.
pub fn fill_four_sided_hole(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Option<CoonsPatchResult> {
    if hole.num_sides != 4 {
        return None;
    }

    // Extract NURBS curves from the four boundary edges
    let mut curves: Vec<NurbsCurve> = Vec::with_capacity(4);
    for &edge_handle in &hole.boundary_edges {
        let edge_data = brep.edges.get(edge_handle)?;
        let curve_on_edge = edge_data.curve.as_ref()?;
        let nurbs = curve_on_edge
            .geometry
            .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end));
        curves.push(nurbs);
    }

    // Assign edges to boundary roles: bottom, right, top, left
    // (following the loop order of the detected hole)
    let boundary = CoonsBoundary {
        bottom: curves[0].clone(),
        right: curves[1].clone(),
        top: curves[2].clone(),
        left: curves[3].clone(),
    };

    construct_coons_patch(&boundary, 8, 8, tolerance)
}

#[cfg(test)]
#[path = "coons_patch_tests.rs"]
mod tests;
