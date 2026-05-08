//! Sphere × cone non-axial dispatch (4th-order case).
//!
//! When the cone axis is not aligned with the sphere's polar axis or
//! the sphere centre is off the cone axis, the existing axial
//! dispatcher (`analytic_trig_sphere_cone`) returns `None` and the
//! intersection is a genuine 4th-degree algebraic curve. This
//! dispatcher samples the curve directly via the **closed-form
//! quadratic in `v`** that the cone-sphere implicit pair admits at
//! each cone-`u` angle:
//!
//!   |apex + v·a + v·t·(cos u·perp1 + sin u·perp2) − sphere_centre|² = r_sphere²
//!
//! Substituting `A = sphere_centre − apex` and `t = tan(half_angle)`:
//!
//!   v²·sec²θ − 2v·(A·a + t·(A·perp1·cos u + A·perp2·sin u))
//!     + (|A|² − r_sphere²) = 0
//!
//! Two roots in v per u (when the discriminant is non-negative);
//! anchor-driven branch selection picks one. Output is NURBS sample-
//! and-fit curve_3d plus NURBS pcurves on each surface, with G1
//! endpoint snapping.
//!
//! Output:
//! - `curve_3d`: NURBS via degree-5 interpolation over 65 cone-u
//!   samples. First/last CPs snapped to anchors.
//! - sphere pcurve: NURBS sample-and-fit in (longitude, colatitude).
//! - cone pcurve: NURBS sample-and-fit in (u_cone, v_cone).
//! - `ParameterMap::Identity` via shared cone-u angle.

use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_sphere_cone_general(
    sphere_center: Point3,
    sphere_radius: f64,
    sphere_is_b: bool,
    cone_apex: Point3,
    cone_axis: Vec3,
    half_angle: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let big_a = sphere_center - cone_apex;
    let a_dot_axis = big_a.dot(&cone_axis);
    let a_dot_p1 = big_a.dot(&perp1);
    let a_dot_p2 = big_a.dot(&perp2);
    let a_len_sq = big_a.dot(&big_a);
    let t = half_angle.tan();
    let sec_sq = 1.0 + t * t;
    let r_sq = sphere_radius * sphere_radius;
    let c0 = a_len_sq - r_sq; // constant term in v²·sec²θ − 2v·B(u) + c0 = 0

    // Helper: at a given cone-u angle, compute the linear-in-v
    // coefficient B(u) and the discriminant. Solve the quadratic.
    let v_at = |u: f64, branch_sign: f64| -> Option<f64> {
        let b = a_dot_axis + t * (a_dot_p1 * u.cos() + a_dot_p2 * u.sin());
        let disc = b * b - sec_sq * c0;
        if disc < 0.0 {
            return None;
        }
        let v = (b + branch_sign * disc.sqrt()) / sec_sq;
        if v <= tol.linear() {
            return None; // negative-v half = other nappe
        }
        Some(v)
    };

    // Solve for the cone u and v at an anchor.
    let solve_anchor = |p: Point3| -> Option<(f64, f64)> {
        let dp = p - cone_apex;
        let v_axis = dp.dot(&cone_axis);
        let perp_vec = dp - cone_axis * v_axis;
        let u = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
        Some((u, v_axis))
    };
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let (u_start, _v_start) = solve_anchor(start)?;
    let (u_end_raw, _v_end) = solve_anchor(end)?;
    let u_end = if u_end_raw < u_start {
        u_end_raw + std::f64::consts::TAU
    } else {
        u_end_raw
    };

    // Pick branch by checking which sign of `disc.sqrt()` reproduces
    // the start anchor's v.
    let b_start = a_dot_axis + t * (a_dot_p1 * u_start.cos() + a_dot_p2 * u_start.sin());
    let disc_start = b_start * b_start - sec_sq * c0;
    if disc_start < -tol.linear() {
        return None;
    }
    let sqrt_disc = disc_start.max(0.0).sqrt();
    let v_pos = (b_start + sqrt_disc) / sec_sq;
    let v_neg = (b_start - sqrt_disc) / sec_sq;
    let branch_sign = if (v_pos - _v_start).abs() <= (v_neg - _v_start).abs() {
        1.0
    } else {
        -1.0
    };

    // Sanity-check anchors against branch.
    let validate_anchor = |p: Point3, side: SurfaceSide| -> Result<(), EdgeTripleError> {
        let off_sphere = (p - sphere_center).length() - sphere_radius;
        if off_sphere.abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: off_sphere,
            });
        }
        let dp = p - cone_apex;
        let v_a = dp.dot(&cone_axis);
        let perp_vec = dp - cone_axis * v_a;
        let r_perp = perp_vec.length();
        let off_cone = (r_perp - v_a * t).abs();
        if off_cone > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: off_cone,
            });
        }
        Ok(())
    };
    if let Err(e) = validate_anchor(start, SurfaceSide::A) {
        return Some(Err(e));
    }
    if let Err(e) = validate_anchor(end, SurfaceSide::B) {
        return Some(Err(e));
    }

    // Sample the curve.
    let n_samples = 65;
    let mut pts_3d: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_sphere: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_cone: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = u_end - u_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let u = u_start + frac * span;
        let v = v_at(u, branch_sign)?;
        let r = v * t;
        let p3 = cone_apex + cone_axis * v + perp1 * (r * u.cos()) + perp2 * (r * u.sin());
        pts_3d.push(p3);
        // Sphere uv: (lon, colat) of (p3 − sphere_centre).
        let dp_s = p3 - sphere_center;
        let lon = dp_s.y.atan2(dp_s.x);
        let colat = (dp_s.z / sphere_radius).clamp(-1.0, 1.0).acos();
        pts_sphere.push(Point3::new(lon, colat, 0.0));
        // Cone uv: (u, v).
        pts_cone.push(Point3::new(u, v, 0.0));
        params.push(frac);
    }
    unwrap_lon(&mut pts_sphere, n_samples);

    let mut nurbs_3d = match crate::geometry::fitting::fit_curve(&pts_3d, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    if let Some(first) = nurbs_3d.control_points.first_mut() {
        *first = start;
    }
    if let Some(last) = nurbs_3d.control_points.last_mut() {
        *last = end;
    }
    let nurbs_sphere = match crate::geometry::fitting::fit_curve(&pts_sphere, &params, 5, n_samples)
    {
        Ok(c) => c,
        Err(_) => return None,
    };
    let nurbs_cone = match crate::geometry::fitting::fit_curve(&pts_cone, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };

    let curve_3d = CurveGeometry::Nurbs(reknot_to(nurbs_3d, u_start, u_end));
    let pcurve_sphere = Pcurve::from_nurbs(reknot_to(nurbs_sphere, u_start, u_end));
    let pcurve_cone = Pcurve::from_nurbs(reknot_to(nurbs_cone, u_start, u_end));

    let sphere_surf = SurfaceGeometry::Sphere {
        center: sphere_center,
        radius: sphere_radius,
    };
    let cone_surf = SurfaceGeometry::Cone {
        apex: cone_apex,
        axis: cone_axis,
        half_angle,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if sphere_is_b {
        (pcurve_cone, pcurve_sphere, &cone_surf, &sphere_surf)
    } else {
        (pcurve_sphere, pcurve_cone, &sphere_surf, &cone_surf)
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: u_start,
        t_end: u_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereCone,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn unwrap_lon(pts: &mut [Point3], n: usize) {
    for i in 1..n {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for cp in pts.iter_mut().skip(i) {
                cp.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for cp in pts.iter_mut().skip(i) {
                cp.x += std::f64::consts::TAU;
            }
        }
    }
}

fn reknot_to(
    mut curve: crate::geometry::nurbs_curve::NurbsCurve,
    t0: f64,
    t1: f64,
) -> crate::geometry::nurbs_curve::NurbsCurve {
    let Some(&k_lo) = curve.knots.first() else {
        return curve;
    };
    let Some(&k_hi) = curve.knots.last() else {
        return curve;
    };
    let span_old = k_hi - k_lo;
    if span_old.abs() < 1.0e-30 {
        return curve;
    }
    let span_new = t1 - t0;
    for k in &mut curve.knots {
        *k = t0 + (*k - k_lo) / span_old * span_new;
    }
    curve
}
