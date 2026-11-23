//! Analytic plane-cone intersection solver.
//!
//! Results can be:
//! - **Circle** (plane perpendicular to cone axis)
//! - **Ellipse** (plane oblique to axis, cutting a closed conic)
//! - **Generator lines** (plane through the apex: `plane_cone_apex`)
//! - **Empty** (plane misses the cone, or degenerate/hyperbola cases)
//!
//! Parabola and hyperbola cases return `None` from the top-level
//! `try_plane_cone` so that the marching fallback is used.

use popkern_math::PortableMath;
use std::f64::consts::FRAC_1_SQRT_2;

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
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
    cvr: (f64, f64),
    tol: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    // A plane through (or very near) the apex cuts generator lines.
    let apex_dist = (apex - po).dot(&pn).abs();
    if apex_dist < tol.linear() {
        return Some(plane_cone_apex::case_through_apex(
            (po, pn, pu, pv),
            apex,
            ca,
            half_angle,
            cvr,
            tol,
            swapped,
        ));
    }

    let cos_alpha = pn.dot(&ca).abs();
    // The circle is the cone's cross-section at the plane's station, so
    // at a tilt alpha it sits off the plane by up to r·sin(alpha): snap
    // only while that stays within the linear tolerance, as plane x
    // cylinder does (PK-845: a tapped hole's mouth cone 1e-4 rad off an
    // inclined face's normal came out a circle 3.3e-4 off the face, and
    // the tap declined on curve-surface validation).
    let sin_alpha = (1.0 - cos_alpha * cos_alpha).max(0.0).sqrt();
    let station = ((po - apex).dot(&pn) / pn.dot(&ca)).abs();
    let circle_radius = station * half_angle.ptan();
    if cos_alpha > 1.0 - tol.angular() && circle_radius * sin_alpha <= tol.linear() {
        case_perpendicular(po, pn, pu, pv, apex, ca, half_angle, cvr, tol, swapped)
    } else {
        case_oblique(po, pn, pu, pv, apex, ca, half_angle, cvr, tol, swapped)
    }
}

#[path = "plane_cone_apex.rs"]
mod plane_cone_apex;
#[path = "plane_cone_oblique.rs"]
mod plane_cone_oblique;
use plane_cone_oblique::case_oblique;

// ── Perpendicular (circle) ──────────────────────────────────────────────

/// Plane normal ≈ parallel to cone axis: intersection is a circle;
/// `None` when its cone pcurve cannot be fitted.
#[allow(clippy::too_many_arguments)]
fn case_perpendicular(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    cvr: (f64, f64),
    tol: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let axis_dot = ca.dot(&pn);
    if axis_dot.abs() < tol.angular() {
        return Some(vec![]);
    }

    // v = signed distance along axis from apex to plane intersection.
    let v = (po - apex).dot(&pn) / axis_dot;
    if v < tol.linear() {
        // Plane intersects at or before the apex (v <= 0): no intersection
        // on the positive half of the cone.
        return Some(vec![]);
    }

    let radius = v * half_angle.ptan();
    let center = apex + ca * v;
    let (p1, p2) = orthonormal_basis(ca);

    let curve_3d = nurbs_ellipse(center, p1, p2, radius, radius);
    let uv_plane = plane_uv_conic(center, p1, p2, radius, radius, po, pu, pv);
    let uv_cone = cone_uv_curve(&curve_3d, apex, ca, half_angle, cvr, tol)?;

    let (uv1, uv2) = if swapped {
        (uv_cone, uv_plane)
    } else {
        (uv_plane, uv_cone)
    };
    Some(vec![SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1: uv1,
        curve_uv2: uv2,
        is_closed: true,
    })])
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

/// The section's pcurve on the cone, co-parameterized with `curve`
/// (`surface(uv(t))` within a thousandth of the linear tolerance of
/// `curve(t)` wherever the curve passes the cone's v window `heights`:
/// the plane x cylinder sections' Hermite fit on the cone's chart). Its control points projected to the chart with the curve's own
/// weights were no such curve: the chart is not affine, so the image
/// strayed off the section between them and its parameter off the
/// curve's (1.3e-2 on a box face's ellipse on a frustum, where a clip
/// crossing read on the pcurve landed 7.6e-3 from the junction its
/// record must reach, and the record was dropped). `None` when the
/// chart cannot be fitted (the marching fallback serves).
fn cone_uv_curve(
    curve: &NurbsCurve,
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    heights: (f64, f64),
    tol: &Tolerance,
) -> Option<NurbsCurve> {
    crate::geometry::intersection::plane_cylinder_pcurve::revolved_conic_pcurve(
        curve,
        apex,
        ca,
        (0.0, half_angle.ptan()),
        Some(heights),
        tol,
    )
}

#[cfg(test)]
#[path = "plane_cone_tests.rs"]
mod tests;
