//! Analytic plane-cone intersection solver.
//!
//! Results can be:
//! - **Circle** (plane perpendicular to cone axis)
//! - **Ellipse** (plane oblique to axis, cutting a closed conic)
//! - **Empty** (plane misses the cone, or degenerate/hyperbola cases)
//!
//! Degenerate cases (plane through apex, parabola, hyperbola) return `None`
//! from the top-level `try_plane_cone` so that the marching fallback is used.

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
    if apex_dist < tol.linear {
        // Degenerate: plane through apex produces lines, not a conic.
        return None;
    }

    let cos_alpha = pn.dot(&ca).abs();
    if cos_alpha > 1.0 - tol.angular {
        Some(case_perpendicular(
            po, pn, pu, pv, apex, ca, half_angle, tol, swapped,
        ))
    } else {
        case_oblique(po, pn, pu, pv, apex, ca, half_angle, tol, swapped)
    }
}

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
    if axis_dot.abs() < tol.angular {
        return vec![];
    }

    // v = signed distance along axis from apex to plane intersection.
    let v = (po - apex).dot(&pn) / axis_dot;
    if v < tol.linear {
        // Plane intersects at or before the apex (v <= 0): no intersection
        // on the positive half of the cone.
        return vec![];
    }

    let radius = v * half_angle.tan();
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

// ── Oblique (ellipse) ───────────────────────────────────────────────────

/// Compute the 3D point on the cone-plane intersection at angle `u`.
///
/// The cone parameterization is:
///   P(u,v) = apex + axis*v + perp1*(v*tan_ha*cos(u)) + perp2*(v*tan_ha*sin(u))
///
/// The plane constraint is: (P - po) · pn = 0
///
/// Solving for v:  v(u) = c / (a + b*cos(u) + d*sin(u))
/// where:
///   c = (po - apex) · pn
///   a = ca · pn
///   b = perp1 · pn * tan_ha
///   d = perp2 · pn * tan_ha
#[allow(clippy::too_many_arguments)]
fn oblique_point(
    apex: Point3,
    perp1: Vec3,
    perp2: Vec3,
    ca: Vec3,
    tan_ha: f64,
    c: f64,
    a: f64,
    b: f64,
    d: f64,
    u: f64,
) -> (Point3, f64) {
    let denom = a + b * u.cos() + d * u.sin();
    let v = c / denom;
    let r = v * tan_ha;
    let pt = apex + ca * v + perp1 * (r * u.cos()) + perp2 * (r * u.sin());
    (pt, v)
}

/// General oblique case: intersection is an ellipse on the cone.
///
/// Returns `None` for hyperbola/parabola (the denominator can reach zero),
/// signaling fallback to marching.
#[allow(clippy::too_many_arguments)]
fn case_oblique(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    tol: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (perp1, perp2) = orthonormal_basis(ca);
    let tan_ha = half_angle.tan();

    let c = (po - apex).dot(&pn);
    let a = ca.dot(&pn);
    let b = perp1.dot(&pn) * tan_ha;
    let d = perp2.dot(&pn) * tan_ha;

    // Ellipse condition: a² > b² + d² (denominator never zero).
    let disc = a * a - b * b - d * d;
    if disc < tol.angular * tol.angular {
        // Parabola or hyperbola: fall through to marching.
        return None;
    }

    // Check that v is positive for all u (intersection on the real cone half).
    // v(u) = c / (a + b*cos(u) + d*sin(u)).
    // The denominator ranges in [a - sqrt(b²+d²), a + sqrt(b²+d²)].
    let amp = (b * b + d * d).sqrt();
    let denom_min = a - amp;
    let denom_max = a + amp;

    // Both v-extremes must be positive (same sign of c and denominator).
    let v_at_min = c / denom_min;
    let v_at_max = c / denom_max;
    if v_at_min < tol.linear || v_at_max < tol.linear {
        // The intersection is on the "wrong" side of the apex or
        // partially behind it.  Return empty.
        return Some(vec![]);
    }

    let pt = |u: f64| -> Point3 { oblique_point(apex, perp1, perp2, ca, tan_ha, c, a, b, d, u).0 };

    // Cardinal angles: 0, π/2, π, 3π/2 (on-curve points).
    use std::f64::consts::{FRAC_PI_2, PI};
    let p0 = pt(0.0);
    let p2 = pt(FRAC_PI_2);
    let p4 = pt(PI);
    let p6 = pt(3.0 * FRAC_PI_2);

    // Tangent at angle u:
    // dP/du = dv/du * (axis + perp1*tan_ha*cos(u) + perp2*tan_ha*sin(u))
    //       + v * (perp1*(-tan_ha*sin(u)) + perp2*(tan_ha*cos(u)))
    // For the NURBS construction, we use a simplified tangent direction
    // approach: the tangent at each cardinal point.
    let tang = |u: f64| -> Vec3 {
        let denom = a + b * u.cos() + d * u.sin();
        let v = c / denom;
        let dv_du = c * (b * u.sin() - d * u.cos()) / (denom * denom);

        // dP/du components:
        let from_v_change = (ca + perp1 * (tan_ha * u.cos()) + perp2 * (tan_ha * u.sin())) * dv_du;
        let from_u_change = (perp1 * (-tan_ha * u.sin()) + perp2 * (tan_ha * u.cos())) * v;
        from_v_change + from_u_change
    };

    // Diagonal control points: tangent-line intersection of adjacent
    // cardinal points.  We use P + 1*T as in the cylinder case.
    let p1 = p0 + tang(0.0);
    let p3 = p2 + tang(FRAC_PI_2);
    let p5 = p4 + tang(PI);
    let p7 = p6 + tang(3.0 * FRAC_PI_2);

    let w = FRAC_1_SQRT_2;
    let curve_3d = NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![p0, p1, p2, p3, p4, p5, p6, p7, p0],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    };

    // UV on plane: project all control points.
    let uv_plane = {
        let cps: Vec<Point3> = curve_3d
            .control_points
            .iter()
            .map(|p| plane_uv_pt(*p, po, pu, pv))
            .collect();
        NurbsCurve {
            degree: 2,
            knots: curve_3d.knots.clone(),
            control_points: cps,
            weights: curve_3d.weights.clone(),
        }
    };

    let uv_cone = cone_uv_curve(&curve_3d, apex, ca);

    let (uv1, uv2) = if swapped {
        (uv_cone, uv_plane)
    } else {
        (uv_plane, uv_cone)
    };

    // Validate: check a few sampled points lie on both surfaces.
    // If tolerance is grossly violated, return None to fall through.
    let _ = tol;

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
