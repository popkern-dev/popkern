//! Analytic plane-cone intersection solver.
//!
//! Results can be:
//! - **Circle** (plane perpendicular to cone axis)
//! - **Ellipse** (plane oblique to axis, cutting a closed conic)
//! - **Empty** (plane misses the cone, or degenerate/hyperbola cases)
//!
//! Degenerate cases (plane through apex, parabola, hyperbola) return `None`
//! from the top-level `try_plane_cone` so that the marching fallback is used.

use popkern_math::PortableMath;
use std::f64::consts::FRAC_1_SQRT_2;

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try to analytically intersect a plane and a cone.
///
/// Returns `Some(intersections)` when the surface pair is a plane and a
/// cone (in either order), or `None` for other surface combinations or
/// degenerate/non-elliptic cases that should fall through to marching.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_plane_cone(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
        ) => intersect_impl(
            *origin,
            *normal,
            *u_axis,
            *v_axis,
            u1_range,
            v1_range,
            *apex,
            *axis,
            *half_angle,
            u2_range,
            v2_range,
            tolerance,
            false,
        ),
        (
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
        ) => intersect_impl(
            *origin,
            *normal,
            *u_axis,
            *v_axis,
            u2_range,
            v2_range,
            *apex,
            *axis,
            *half_angle,
            u1_range,
            v1_range,
            tolerance,
            true,
        ),
        _ => None,
    }
}

/// Dispatch to perpendicular or oblique based on the angle between the
/// plane normal and the cone axis.  Returns `None` for degenerate cases.
#[allow(clippy::too_many_arguments)]
fn intersect_impl(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    _pur: (f64, f64),
    _pvr: (f64, f64),
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    _cur: (f64, f64),
    _cvr: (f64, f64),
    tol: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    // Check if the plane passes through (or very near) the apex.
    let apex_dist = (apex - po).dot(&pn).abs();
    if apex_dist < tol.linear() {
        // Degenerate: plane through apex produces lines, not a conic.
        return None;
    }

    let cos_alpha = pn.dot(&ca).abs();
    if cos_alpha > 1.0 - tol.angular() {
        Some(case_perpendicular(
            po, pn, pu, pv, apex, ca, half_angle, tol, swapped,
        ))
    } else {
        case_oblique(po, pn, pu, pv, apex, ca, half_angle, tol, swapped)
    }
}

#[path = "plane_cone_oblique.rs"]
mod plane_cone_oblique;
use plane_cone_oblique::case_oblique;

// ── Perpendicular (circle) ──────────────────────────────────────────────

/// Plane normal ≈ parallel to cone axis: intersection is a circle.
#[allow(clippy::too_many_arguments)]
fn case_perpendicular(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    tol: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    let axis_dot = ca.dot(&pn);
    if axis_dot.abs() < tol.angular() {
        return vec![];
    }

    // v = signed distance along axis from apex to plane intersection.
    let v = (po - apex).dot(&pn) / axis_dot;
    if v < tol.linear() {
        // Plane intersects at or before the apex (v <= 0): no intersection
        // on the positive half of the cone.
        return vec![];
    }

    let radius = v * half_angle.ptan();
    let center = apex + ca * v;
    let (p1, p2) = orthonormal_basis(ca);

    let curve_3d = nurbs_ellipse(center, p1, p2, radius, radius);
    let uv_plane = plane_uv_conic(center, p1, p2, radius, radius, po, pu, pv);
    let uv_cone = cone_uv_curve(&curve_3d, apex, ca);

    let (uv1, uv2) = if swapped {
        (uv_cone, uv_plane)
    } else {
        (uv_plane, uv_cone)
    };
    vec![SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1: uv1,
        curve_uv2: uv2,
        is_closed: true,
    })]
}

// ── NURBS builders ──────────────────────────────────────────────────────

/// Degree-2 rational NURBS ellipse (or circle when `sa == sb`).
fn nurbs_ellipse(c: Point3, a1: Vec3, a2: Vec3, sa: f64, sb: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            c + a1 * sa,
            c + a1 * sa + a2 * sb,
            c + a2 * sb,
            c - a1 * sa + a2 * sb,
            c - a1 * sa,
            c - a1 * sa - a2 * sb,
            c - a2 * sb,
            c + a1 * sa - a2 * sb,
            c + a1 * sa,
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Project an ellipse/circle's control points to a plane's UV space.
#[allow(clippy::too_many_arguments)]
fn plane_uv_conic(
    c: Point3,
    a1: Vec3,
    a2: Vec3,
    sa: f64,
    sb: f64,
    po: Point3,
    pu: Vec3,
    pv: Vec3,
) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    let offsets = [
        a1 * sa,
        a1 * sa + a2 * sb,
        a2 * sb,
        a1 * (-sa) + a2 * sb,
        a1 * (-sa),
        a1 * (-sa) + a2 * (-sb),
        a2 * (-sb),
        a1 * sa + a2 * (-sb),
        a1 * sa,
    ];
    let cps: Vec<Point3> = offsets
        .iter()
        .map(|o| plane_uv_pt(c + *o, po, pu, pv))
        .collect();
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: cps,
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Project a 3D point to plane UV, returned as `Point3(u, v, 0)`.
///
/// Uses the correct inverse parameterization:
/// `u = (p-origin)·u_axis / |u_axis|²`, `v = (p-origin)·v_axis / |v_axis|²`.
fn plane_uv_pt(pt: Point3, o: Point3, pu: Vec3, pv: Vec3) -> Point3 {
    let d = pt - o;
    let u_len_sq = pu.dot(&pu);
    let v_len_sq = pv.dot(&pv);
    let u = if u_len_sq > 1e-30 {
        d.dot(&pu) / u_len_sq
    } else {
        d.dot(&pu)
    };
    let v = if v_len_sq > 1e-30 {
        d.dot(&pv) / v_len_sq
    } else {
        d.dot(&pv)
    };
    Point3::new(u, v, 0.0)
}

/// Project a 3D curve's control points to cone UV, unwrapping u.
fn cone_uv_curve(curve: &NurbsCurve, apex: Point3, ca: Vec3) -> NurbsCurve {
    let cone = SurfaceGeometry::Cone {
        apex,
        axis: ca,
        half_angle: 0.0, // half_angle not used by project_point_to_uv for cones
    };
    let mut cps: Vec<Point3> = curve
        .control_points
        .iter()
        .map(|pt| {
            let (u, v) = project_point_to_uv(&cone, pt).unwrap_or((0.0, 0.0));
            Point3::new(u, v, 0.0)
        })
        .collect();
    // Unwrap u so the NURBS curve is continuous across the TAU seam.
    let tau = std::f64::consts::TAU;
    for i in 1..cps.len() {
        let diff = cps[i].x - cps[i - 1].x;
        if diff > tau * 0.5 {
            cps[i].x -= tau;
        } else if diff < -tau * 0.5 {
            cps[i].x += tau;
        }
    }
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: cps,
        weights: curve.weights.clone(),
    }
}

#[cfg(test)]
#[path = "plane_cone_tests.rs"]
mod tests;
