//! Cylinder × cone non-coaxial dispatch (4th-order case).
//!
//! Two quadrics with axes neither aligned-and-coincident: the
//! intersection is a 4th-degree algebraic curve. The existing
//! coaxial dispatcher (`analytic_trig_cyl_cone`) returns `None`
//! when either the axes are not parallel OR the cone apex is off the
//! cylinder axis; this dispatcher handles the general case.
//!
//! Closed-form sampling: for each cone-`u` angle, the cylinder
//! constraint `|perp_to_cyl_axis(P − O_cyl)|² = r_cyl²` produces a
//! quadratic in the cone-`v` parameter:
//!
//!   A(u)·v² + B(u)·v + C = 0
//!
//! where
//!
//!   A(u) = sec²θ − (a·a_cyl + t·(p1·a_cyl·cos u + p2·a_cyl·sin u))²
//!   B(u) = 2[(E·a − E·a_cyl·a·a_cyl)
//!          + t·((E·p1 − E·a_cyl·p1·a_cyl)·cos u
//!             + (E·p2 − E·a_cyl·p2·a_cyl)·sin u)]
//!   C    = |E|² − (E·a_cyl)² − r_cyl²
//!
//! with `E = apex − O_cyl`, `t = tan(θ)`, `(p1, p2)` an orthonormal
//! basis perpendicular to the cone axis `a`, and `a_cyl` the
//! cylinder axis (unit). Two roots per u; anchor-driven branch
//! selection picks one.
//!
//! Output:
//! - `curve_3d`: NURBS via degree-5 interpolation over 65 cone-u
//!   samples. First/last CPs snapped to anchors for G1.
//! - cyl pcurve: NURBS sample-and-fit in (u_cyl, v_cyl).
//! - cone pcurve: NURBS sample-and-fit in (u_cone, v_cone).
//! - `ParameterMap::Identity` over cone-u.

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
pub(super) fn try_cyl_cone_general(
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
    cyl_is_b: bool,
    cone_apex: Point3,
    cone_axis: Vec3,
    half_angle: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let big_e = cone_apex - cyl_origin;
    let e_a = big_e.dot(&cone_axis);
    let e_axial = big_e.dot(&cyl_axis);
    let e_p1 = big_e.dot(&perp1);
    let e_p2 = big_e.dot(&perp2);
    let e_len_sq = big_e.dot(&big_e);
    let alpha_axial = cone_axis.dot(&cyl_axis);
    let beta_p1 = perp1.dot(&cyl_axis);
    let beta_p2 = perp2.dot(&cyl_axis);
    let t = half_angle.tan();
    let sec_sq = 1.0 + t * t;
    let r_sq = cyl_radius * cyl_radius;
    let c = e_len_sq - e_axial * e_axial - r_sq;

    // Constants used in B(u): linear-in-(cos u, sin u) coefficients.
    let b_const = 2.0 * (e_a - e_axial * alpha_axial);
    let b_cos = 2.0 * t * (e_p1 - e_axial * beta_p1);
    let b_sin = 2.0 * t * (e_p2 - e_axial * beta_p2);

    let v_at = |u: f64, branch_sign: f64| -> Option<f64> {
        let cu = u.cos();
        let su = u.sin();
        let inner = alpha_axial + t * (beta_p1 * cu + beta_p2 * su);
        let a_u = sec_sq - inner * inner;
        if a_u.abs() < tol.solver() {
            return None; // degenerate at this u (linear-in-v case)
        }
        let b_u = b_const + b_cos * cu + b_sin * su;
        let disc = b_u * b_u - 4.0 * a_u * c;
        if disc < 0.0 {
            return None;
        }
        let v = (-b_u + branch_sign * disc.sqrt()) / (2.0 * a_u);
        if v <= tol.linear() {
            return None;
        }
        Some(v)
    };

    // Solve for cone-(u, v) at an anchor and validate it's on both
    // surfaces.
    let solve_anchor = |p: Point3| -> Option<(f64, f64)> {
        let dp = p - cone_apex;
        let v_axis = dp.dot(&cone_axis);
        let perp_vec = dp - cone_axis * v_axis;
        let u = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
        Some((u, v_axis))
    };
    let validate_anchor = |p: Point3, side: SurfaceSide| -> Result<(), EdgeTripleError> {
        let dp_cyl = p - cyl_origin;
        let v_cyl = dp_cyl.dot(&cyl_axis);
        let perp_cyl = dp_cyl - cyl_axis * v_cyl;
        let off_cyl = (perp_cyl.length() - cyl_radius).abs();
        if off_cyl > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: off_cyl,
            });
        }
        let dp_cone = p - cone_apex;
        let v_cone = dp_cone.dot(&cone_axis);
        let perp_cone = dp_cone - cone_axis * v_cone;
        let off_cone = (perp_cone.length() - v_cone * t).abs();
        if off_cone > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: off_cone,
            });
        }
        Ok(())
    };

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    if let Err(e) = validate_anchor(start, SurfaceSide::A) {
        return Some(Err(e));
    }
    if let Err(e) = validate_anchor(end, SurfaceSide::B) {
        return Some(Err(e));
    }
    let (u_start, v_start_truth) = solve_anchor(start)?;
    let (u_end_raw, _v_end_truth) = solve_anchor(end)?;
    let u_end = if u_end_raw < u_start {
        u_end_raw + std::f64::consts::TAU
    } else {
        u_end_raw
    };

    // Branch sign by which root reproduces the start anchor's v.
    let cu = u_start.cos();
    let su = u_start.sin();
    let inner = alpha_axial + t * (beta_p1 * cu + beta_p2 * su);
    let a_u = sec_sq - inner * inner;
    if a_u.abs() < tol.solver() {
        return None;
    }
    let b_u = b_const + b_cos * cu + b_sin * su;
    let disc = b_u * b_u - 4.0 * a_u * c;
    if disc < -tol.linear() {
        return None;
    }
    let sqrt_disc = disc.max(0.0).sqrt();
    let v_pos = (-b_u + sqrt_disc) / (2.0 * a_u);
    let v_neg = (-b_u - sqrt_disc) / (2.0 * a_u);
    let branch_sign = if (v_pos - v_start_truth).abs() <= (v_neg - v_start_truth).abs() {
        1.0
    } else {
        -1.0
    };

    // Sample the curve.
    let n_samples = 65;
    let mut pts_3d: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_cyl: Vec<Point3> = Vec::with_capacity(n_samples);
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
        // Cyl uv.
        let dp_cyl = p3 - cyl_origin;
        let v_cyl_coord = dp_cyl.dot(&cyl_axis);
        let (cyl_p1, cyl_p2) = crate::geometry::analytic_surfaces::orthonormal_basis(cyl_axis);
        let perp_vec = dp_cyl - cyl_axis * v_cyl_coord;
        let u_cyl = perp_vec.dot(&cyl_p2).atan2(perp_vec.dot(&cyl_p1));
        pts_cyl.push(Point3::new(u_cyl, v_cyl_coord, 0.0));
        // Cone uv.
        pts_cone.push(Point3::new(u, v, 0.0));
        params.push(frac);
    }
    unwrap_u(&mut pts_cyl, n_samples);

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
    let nurbs_cyl = match crate::geometry::fitting::fit_curve(&pts_cyl, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let nurbs_cone = match crate::geometry::fitting::fit_curve(&pts_cone, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };

    let curve_3d = CurveGeometry::Nurbs(reknot_to(nurbs_3d, u_start, u_end));
    let pcurve_cyl = Pcurve::from_nurbs(reknot_to(nurbs_cyl, u_start, u_end));
    let pcurve_cone = Pcurve::from_nurbs(reknot_to(nurbs_cone, u_start, u_end));

    let cyl_surf = SurfaceGeometry::Cylinder {
        origin: cyl_origin,
        axis: cyl_axis,
        radius: cyl_radius,
    };
    let cone_surf = SurfaceGeometry::Cone {
        apex: cone_apex,
        axis: cone_axis,
        half_angle,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if cyl_is_b {
        (pcurve_cone, pcurve_cyl, &cone_surf, &cyl_surf)
    } else {
        (pcurve_cyl, pcurve_cone, &cyl_surf, &cone_surf)
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
            kind: AnalyticPair::CylinderCone,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn unwrap_u(pts: &mut [Point3], n: usize) {
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
    if curve.knots.is_empty() {
        return curve;
    }
    let k_lo = *curve.knots.first().unwrap();
    let k_hi = *curve.knots.last().unwrap();
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
