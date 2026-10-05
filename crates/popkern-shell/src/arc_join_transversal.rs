//! Transversal connector authoring for mixed corners (PK-326 E1):
//! the connector edge of the two-sector peel lies on BOTH cut-wedge
//! surfaces. Exact conic arms (plane × cylinder ellipse, equal-radius
//! cylinder × cylinder bisector-plane ellipse) with a marched
//! alternating-projection NURBS fit fallback under a measured
//! on-surface residual gate.

use crate::curved_recompute::project_to_surface;
use popkern_core::geometry::analytic_curves::circle_basis;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_math::PortableMath;

/// The transversal intersection arc between two join/offset surfaces,
/// from `a` to `b` (PK-326 mixed corners): the connector edge of the
/// two-sector peel lies on BOTH cut-wedge surfaces. Exact conic arms:
/// plane × cylinder is an ellipse in the plane (a circle when the axis
/// is normal to it), and two equal-radius cylinders with intersecting
/// axes intersect in two planar ellipses (the bisector planes — pick
/// the branch through `a` and `b`). General pairs march the
/// intersection by alternating projection and fit a NURBS with a
/// measured residual gate.
pub(crate) fn transversal_arc(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    a: Point3,
    b: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let gate = (tolerance.linear() * 10.0).max(a.distance_to(&b) * 1.0e-9);
    let exact = plane_cylinder_arc(s1, s2, a, b)
        .or_else(|| plane_cylinder_arc(s2, s1, a, b))
        .or_else(|| cylinder_cylinder_arc(s1, s2, a, b, tolerance));
    if let Some(curve) = exact {
        // Authoring self-check: endpoints reproduce and the midpoint
        // sits on both surfaces.
        let c0 = curve.geometry.point_at(curve.t_start);
        let c1 = curve.geometry.point_at(curve.t_end);
        let ends_ok = (c0.distance_to(&a).max(c1.distance_to(&b)) <= gate)
            || (c0.distance_to(&b).max(c1.distance_to(&a)) <= gate);
        let mid = curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end));
        let on_both = crate::curved_recompute::surface_distance(s1, &mid)
            .is_some_and(|d| d <= gate)
            && crate::curved_recompute::surface_distance(s2, &mid).is_some_and(|d| d <= gate);
        if ends_ok && on_both {
            return Some(curve);
        }
    }
    marched_intersection_fit(s1, s2, a, b, gate)
}

/// Plane × cylinder intersection through `a`/`b`: authored directly in
/// the [`CurveGeometry::Ellipse`] basis form — the cylinder's angular
/// parameter is the ellipse parameter, so endpoint angles are exact.
fn plane_cylinder_arc(
    plane: &SurfaceGeometry,
    cylinder: &SurfaceGeometry,
    a: Point3,
    b: Point3,
) -> Option<CurveOnEdge> {
    let SurfaceGeometry::Plane { origin, normal, .. } = plane else {
        return None;
    };
    let SurfaceGeometry::Cylinder {
        origin: c_origin,
        axis,
        radius,
    } = cylinder
    else {
        return None;
    };
    let n = normal.normalized()?;
    let ax = axis.normalized()?;
    let denom = n.dot(&ax);
    if denom.abs() < 1.0e-9 {
        return None; // axis parallel to the plane: line pair, not a conic
    }
    let k = n.dot(&Vec3::new(origin.x, origin.y, origin.z));
    let (e1, e2) = circle_basis(ax);
    let o_vec = Vec3::new(c_origin.x, c_origin.y, c_origin.z);
    let center = *c_origin + ax * ((k - n.dot(&o_vec)) / denom);
    let u_basis = (e1 - ax * (n.dot(&e1) / denom)) * *radius;
    let v_basis = (e2 - ax * (n.dot(&e2) / denom)) * *radius;

    // Endpoint angles from the CYLINDER parametrization (shared with
    // the ellipse parameter by construction).
    let angle_of = |p: Point3| -> Option<f64> {
        let w = p - *c_origin;
        let radial = w - ax * w.dot(&ax);
        let r_len = radial.length();
        if r_len < 1.0e-12 {
            return None;
        }
        Some(radial.dot(&e2).patan2(radial.dot(&e1)))
    };
    let theta_a = angle_of(a)?;
    let theta_b = angle_of(b)?;
    minor_sweep_curve(
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: 0.0,
            end_angle: 0.0,
        },
        theta_a,
        theta_b,
    )
}

/// Equal-radius cylinders with intersecting, non-parallel axes meet in
/// two planar ellipses lying in the bisector planes; reduce to the
/// plane × cylinder arm with the bisector plane containing `a` and `b`.
fn cylinder_cylinder_arc(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    a: Point3,
    b: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let SurfaceGeometry::Cylinder {
        origin: o1,
        axis: ax1,
        radius: r1,
    } = s1
    else {
        return None;
    };
    let SurfaceGeometry::Cylinder {
        origin: o2,
        axis: ax2,
        radius: r2,
    } = s2
    else {
        return None;
    };
    let tol = tolerance.linear();
    if (r1 - r2).abs() > tol {
        return None;
    }
    let a1 = ax1.normalized()?;
    let a2 = ax2.normalized()?;
    let cross = a1.cross(&a2);
    if cross.length() < 1.0e-9 {
        return None; // parallel axes: no transversal conic
    }
    // Closest points between the two axis lines must coincide.
    let w = *o2 - *o1;
    let n_len_sq = cross.dot(&cross);
    let t1 = w.cross(&a2).dot(&cross) / n_len_sq;
    let t2 = w.cross(&a1).dot(&cross) / n_len_sq;
    let p1 = *o1 + a1 * t1;
    let p2 = *o2 + a2 * t2;
    if p1.distance_to(&p2) > tol {
        return None; // skew axes: non-planar quartic, use the fit
    }
    let apex = Point3::new(
        0.5 * (p1.x + p2.x),
        0.5 * (p1.y + p2.y),
        0.5 * (p1.z + p2.z),
    );
    // The two bisector planes have normals a1−a2 and a1+a2; pick the
    // one containing both endpoints.
    for n in [a1 - a2, a1 + a2] {
        let Some(n) = n.normalized() else { continue };
        let da = (a - apex).dot(&n).abs();
        let db = (b - apex).dot(&n).abs();
        if da <= tol * 10.0 && db <= tol * 10.0 {
            let plane = SurfaceGeometry::Plane {
                origin: apex,
                normal: n,
                u_axis: circle_basis(n).0,
                v_axis: circle_basis(n).1,
            };
            return plane_cylinder_arc(&plane, s1, a, b);
        }
    }
    None
}

/// Wrap the sweep to the minor arc and stamp the angle range.
fn minor_sweep_curve(geometry: CurveGeometry, theta_a: f64, theta_b: f64) -> Option<CurveOnEdge> {
    use std::f64::consts::{PI, TAU};
    let mut delta = (theta_b - theta_a) % TAU;
    if delta > PI {
        delta -= TAU;
    } else if delta < -PI {
        delta += TAU;
    }
    if delta.abs() >= PI - 1.0e-9 || delta.abs() < 1.0e-12 {
        return None; // antipodal or degenerate sweep
    }
    let (start, end) = if delta >= 0.0 {
        (theta_a, theta_a + delta)
    } else {
        (theta_b, theta_b - delta)
    };
    let geometry = match geometry {
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: start,
            end_angle: end,
        },
        other => other,
    };
    Some(CurveOnEdge {
        geometry,
        t_start: start,
        t_end: end,
    })
}

/// General fallback: march the intersection by alternating projection
/// from chord seeds and fit a NURBS through the converged samples,
/// gating the on-surface residual.
fn marched_intersection_fit(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    a: Point3,
    b: Point3,
    gate: f64,
) -> Option<CurveOnEdge> {
    const N: usize = 17;
    let mut samples: Vec<Point3> = Vec::with_capacity(N);
    for i in 0..N {
        let f = i as f64 / (N - 1) as f64;
        let mut p = Point3::new(
            a.x + f * (b.x - a.x),
            a.y + f * (b.y - a.y),
            a.z + f * (b.z - a.z),
        );
        for _ in 0..40 {
            let q1 = project_to_surface(s1, &p)?;
            let q2 = project_to_surface(s2, &q1)?;
            let moved = p.distance_to(&q2);
            p = q2;
            if moved < 1.0e-14 {
                break;
            }
        }
        if crate::curved_recompute::surface_distance(s1, &p)? > gate
            || crate::curved_recompute::surface_distance(s2, &p)? > gate
        {
            return None;
        }
        samples.push(p);
    }
    // Pin the exact endpoints.
    samples[0] = a;
    samples[N - 1] = b;
    let nurbs = popkern_blend::center_curve::fit_nurbs_through_points(&samples).ok()?;
    let (t0, t1) = nurbs.parameter_range();
    for i in 0..16 {
        let t = t0 + (i as f64 + 0.5) / 16.0 * (t1 - t0);
        let p = nurbs.point_at(t);
        if crate::curved_recompute::surface_distance(s1, &p)? > gate * 2.0
            || crate::curved_recompute::surface_distance(s2, &p)? > gate * 2.0
        {
            return None;
        }
    }
    Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: t0,
        t_end: t1,
    })
}
