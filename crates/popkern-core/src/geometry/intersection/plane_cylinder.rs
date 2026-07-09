//! Analytic plane-cylinder intersection solver.
//!
//! Results can be:
//! - **Circle** (plane perpendicular to cylinder axis)
//! - **Ellipse** (plane oblique to cylinder axis)
//! - **Two lines** (plane parallel to cylinder axis and within radius)
//! - **Empty** (plane parallel to axis but beyond radius, or tangent)

use std::f64::consts::FRAC_1_SQRT_2;

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try to analytically intersect a plane and a cylinder.
///
/// Returns `Some(intersections)` when the surface pair is a plane and a
/// cylinder (in either order), or `None` for other surface combinations.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_plane_cylinder(
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
            SurfaceGeometry::Cylinder {
                origin: co,
                axis,
                radius,
            },
        ) => Some(intersect_impl(
            *origin, *normal, *u_axis, *v_axis, u1_range, v1_range, *co, *axis, *radius, u2_range,
            v2_range, tolerance, false,
        )),
        (
            SurfaceGeometry::Cylinder {
                origin: co,
                axis,
                radius,
            },
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
        ) => Some(intersect_impl(
            *origin, *normal, *u_axis, *v_axis, u2_range, v2_range, *co, *axis, *radius, u1_range,
            v1_range, tolerance, true,
        )),
        _ => None,
    }
}

/// Dispatch to perpendicular / parallel / oblique based on the angle.
#[allow(clippy::too_many_arguments)]
fn intersect_impl(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    _pur: (f64, f64),
    _pvr: (f64, f64),
    co: Point3,
    ca: Vec3,
    cr: f64,
    _cur: (f64, f64),
    _cvr: (f64, f64),
    tol: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    let cos_alpha = pn.dot(&ca).abs();
    if cos_alpha > 1.0 - tol.angular() {
        case_perpendicular(po, pn, pu, pv, co, ca, cr, tol, swapped)
    } else if cos_alpha < tol.angular() {
        case_parallel(po, pn, pu, pv, co, ca, cr, _cvr, tol, swapped)
    } else {
        case_oblique(po, pn, pu, pv, co, ca, cr, tol, swapped)
    }
}

// ── Perpendicular (circle) ──────────────────────────────────────────────

/// Plane normal ≈ parallel to cylinder axis: intersection is a circle.
#[allow(clippy::too_many_arguments)]
fn case_perpendicular(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    co: Point3,
    ca: Vec3,
    cr: f64,
    tol: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    let axis_dot = ca.dot(&pn);
    if axis_dot.abs() < tol.angular() {
        return vec![];
    }
    let t = (po - co).dot(&pn) / axis_dot;
    let center = co + ca * t;
    let (p1, p2) = orthonormal_basis(ca);

    let curve_3d = nurbs_ellipse(center, p1, p2, cr, cr);
    let uv_plane = plane_uv_conic(center, p1, p2, cr, cr, po, pu, pv);
    let uv_cyl = nurbs_hline(t);

    let (uv1, uv2) = if swapped {
        (uv_cyl, uv_plane)
    } else {
        (uv_plane, uv_cyl)
    };
    vec![SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1: uv1,
        curve_uv2: uv2,
        is_closed: true,
    })]
}

// ── Parallel (two lines / empty) ────────────────────────────────────────

/// Plane normal ≈ perpendicular to cylinder axis.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_arguments)]
fn case_parallel(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    co: Point3,
    ca: Vec3,
    cr: f64,
    cvr: (f64, f64),
    tol: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    let d = (co - po).dot(&pn);
    if d.abs() > cr + tol.linear() {
        return vec![];
    }
    if d.abs() > cr - tol.linear() {
        return vec![];
    } // tangent

    let (p1, p2) = orthonormal_basis(ca);
    let n1 = pn.dot(&p1);
    let n2 = pn.dot(&p2);
    let amp = (n1 * n1 + n2 * n2).sqrt();
    if amp < tol.angular() {
        return vec![];
    }

    let ratio = (-d / cr) / amp;
    if ratio.abs() > 1.0 {
        return vec![];
    }
    // Solve d + cr·(n1·cos u + n2·sin u) = 0. With n1 = amp·cos(base)
    // and n2 = amp·sin(base) this is amp·cos(u − base) = −d/cr, so
    // u = base ± acos(ratio). (The previous π ± δ − base form placed
    // the lines up to a diameter off the plane — the trim-honesty
    // gate rejected them and the real intersection went missing.)
    let base = f64::atan2(n2, n1);
    let delta = ratio.clamp(-1.0, 1.0).acos();

    let angles = [norm_angle(base + delta), norm_angle(base - delta)];

    // Span the caller's axial window with headroom; the ±10 constant
    // truncated the lines to 20 mm and the clip lost everything
    // beyond it.
    let pad = (cvr.1 - cvr.0).abs().max(1.0);
    let (v_lo, v_hi) = (cvr.0 - pad, cvr.1 + pad);

    angles
        .iter()
        .map(|&u_ang| {
            let lp = co + p1 * (cr * u_ang.cos()) + p2 * (cr * u_ang.sin());
            let a = lp + ca * v_lo;
            let b = lp + ca * v_hi;
            let curve_3d = nurbs_line(a, b);
            let uv_cyl = nurbs_line(Point3::new(u_ang, v_lo, 0.0), Point3::new(u_ang, v_hi, 0.0));
            let uv_plane = nurbs_line(plane_uv_pt(a, po, pu, pv), plane_uv_pt(b, po, pu, pv));
            let (uv1, uv2) = if swapped {
                (uv_cyl, uv_plane)
            } else {
                (uv_plane, uv_cyl)
            };
            SurfaceSurfaceIntersection::Curve(CurveOnSurface {
                curve_3d,
                curve_uv1: uv1,
                curve_uv2: uv2,
                is_closed: false,
            })
        })
        .collect()
}

// ── Oblique (ellipse) ───────────────────────────────────────────────────

/// Compute the 3D point on the cylinder-plane intersection at angle `u`.
///
/// `d0 = (po - co) . pn`, `n1 = perp1 . pn`, `n2 = perp2 . pn`,
/// `nd = ca . pn`.  Height: `h(u) = (d0 - cr*cos(u)*n1 - cr*sin(u)*n2) / nd`.
#[allow(clippy::too_many_arguments)]
fn oblique_point(
    co: Point3,
    perp1: Vec3,
    perp2: Vec3,
    ca: Vec3,
    cr: f64,
    d0: f64,
    n1: f64,
    n2: f64,
    nd: f64,
    u: f64,
) -> Point3 {
    let h = (d0 - cr * u.cos() * n1 - cr * u.sin() * n2) / nd;
    co + perp1 * (cr * u.cos()) + perp2 * (cr * u.sin()) + ca * h
}

/// General oblique case: intersection is an ellipse on the cylinder.
///
/// Builds a degree-2 rational NURBS curve with 9 control points that lies
/// exactly on both the cylinder and the plane.  Cardinal control points
/// (weight 1) sit on the curve; diagonal control points (weight 1/sqrt(2))
/// are tangent-line intersections, exactly as in the standard rational
/// Bezier circle construction.
#[allow(clippy::too_many_arguments)]
fn case_oblique(
    po: Point3,
    pn: Vec3,
    pu: Vec3,
    pv: Vec3,
    co: Point3,
    ca: Vec3,
    cr: f64,
    tol: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    let nd = ca.dot(&pn);
    if nd.abs() < tol.angular() {
        return vec![];
    }

    let (perp1, perp2) = orthonormal_basis(ca);
    let d0 = (po - co).dot(&pn);
    let n1 = perp1.dot(&pn);
    let n2 = perp2.dot(&pn);

    let pt = |u: f64| oblique_point(co, perp1, perp2, ca, cr, d0, n1, n2, nd, u);

    // Cardinal angles: 0, pi/2, pi, 3pi/2, 2pi (on-curve points).
    use std::f64::consts::{FRAC_PI_2, PI};
    let p0 = pt(0.0);
    let p2 = pt(FRAC_PI_2);
    let p4 = pt(PI);
    let p6 = pt(3.0 * FRAC_PI_2);

    // Tangent at angle u: d/du of the curve.
    // radial part: perp1*(-cr*sin u) + perp2*(cr*cos u)
    // height part: ca * (cr*sin(u)*n1 - cr*cos(u)*n2) / nd
    let tang = |u: f64| -> Vec3 {
        let dr = perp1 * (-cr * u.sin()) + perp2 * (cr * u.cos());
        let dh = ca * ((cr * u.sin() * n1 - cr * u.cos() * n2) / nd);
        dr + dh
    };

    // Diagonal control points: tangent-line intersection of adjacent
    // cardinal points.  For the standard circle NURBS, these are at
    // P_i + 1*T_i = P_{i+1} - 1*T_{i+1} (the tangent parameter is 1
    // because the radial part works out to r*sqrt(2) at 45 degrees).
    // We simply compute P_cardinal + 1 * tangent_at_cardinal, which
    // gives the correct off-curve control point for any height variation
    // that is linear in cos/sin.
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

    let uv_cyl = cylinder_uv_curve(&curve_3d, co, ca, cr);

    let (uv1, uv2) = if swapped {
        (uv_cyl, uv_plane)
    } else {
        (uv_plane, uv_cyl)
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

/// Horizontal line in cylinder UV at constant v.
fn nurbs_hline(v: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, v, 0.0),
            Point3::new(std::f64::consts::TAU, v, 0.0),
        ],
        weights: vec![1.0, 1.0],
    }
}

/// Degree-1 NURBS line.
fn nurbs_line(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Project a 3D point to plane UV, returned as `Point3(u, v, 0)`.
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

/// Normalize angle to [0, TAU).
fn norm_angle(a: f64) -> f64 {
    let tau = std::f64::consts::TAU;
    let mut r = a % tau;
    if r < 0.0 {
        r += tau;
    }
    r
}

/// Project a 3D curve's control points to cylinder UV, unwrapping u.
fn cylinder_uv_curve(curve: &NurbsCurve, co: Point3, ca: Vec3, cr: f64) -> NurbsCurve {
    let cyl = SurfaceGeometry::Cylinder {
        origin: co,
        axis: ca,
        radius: cr,
    };
    let mut cps: Vec<Point3> = curve
        .control_points
        .iter()
        .map(|pt| {
            let (u, v) = project_point_to_uv(&cyl, pt).unwrap_or((0.0, 0.0));
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
#[path = "plane_cylinder_tests.rs"]
mod tests;
