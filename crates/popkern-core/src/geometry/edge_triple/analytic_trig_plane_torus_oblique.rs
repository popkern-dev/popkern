//! Plane × torus generic-oblique dispatch (4th-order case).
//!
//! Handles plane × torus intersections that are neither
//! axis-perpendicular (handled by `analytic_trig_plane_torus`) nor
//! Villarceau-bitangent (handled by
//! `analytic_trig_plane_torus_villarceau`). The intersection is
//! genuinely 4th-degree.
//!
//! Closed-form sampling: substitute the torus parameterisation
//!
//!   P(u, v) = T_centre + (R + r·cos v)·(cos u·p1 + sin u·p2) + r·sin v·a
//!
//! into the plane equation `(P − po)·n = 0`. With
//!
//!   d   = (po − T_centre)·n
//!   cn  = a·n
//!   p1n = p1·n   p2n = p2·n
//!   X(u) = p1n·cos u + p2n·sin u
//!   M(u) = √(X(u)² + cn²)
//!   φ(u) = atan2(cn, X(u))
//!
//! the constraint reduces to
//!
//!   cos(v − φ(u)) = (d − R·X(u)) / (r·M(u))
//!
//! For each cone-style u where the right-hand side is in `[−1, 1]`,
//! two real `v` values exist:
//!
//!   v = φ(u) ± acos((d − R·X(u)) / (r·M(u)))
//!
//! No Newton needed. Anchor-driven branch selection picks the sign.
//! u-range is restricted to where the rhs fits `[−1, 1]`; outside
//! the curve has no intersection on this branch.
//!
//! Output:
//! - `curve_3d`: NURBS via degree-5 interpolation over 65 u-samples,
//!   anchor-snapped endpoints.
//! - plane pcurve: NURBS sample-and-fit in plane uv.
//! - torus pcurve: NURBS sample-and-fit in (u_torus, v_torus).
//! - `ParameterMap::Identity` over u.

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
pub(super) fn try_plane_torus_oblique(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_is_b: bool,
    torus_center: Point3,
    torus_axis: Vec3,
    major_radius: f64,
    minor_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let n_len_sq = plane_normal.dot(&plane_normal);
    if n_len_sq < tol.solver() {
        return None;
    }
    let n_unit = plane_normal * (1.0 / n_len_sq.sqrt());
    let cn = torus_axis.dot(&n_unit);
    // Reject perpendicular case (handled by axis-perpendicular dispatcher).
    if cn.abs() >= 1.0 - tol.angular() {
        return None;
    }

    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(torus_axis);
    let p1n = perp1.dot(&n_unit);
    let p2n = perp2.dot(&n_unit);
    let d = (plane_origin - torus_center).dot(&n_unit);
    let r = minor_radius;
    let big_r = major_radius;

    let v_at = |u: f64, branch_sign: f64| -> Option<f64> {
        let cu = u.cos();
        let su = u.sin();
        let x = p1n * cu + p2n * su;
        let m = (x * x + cn * cn).sqrt();
        if m < tol.solver() {
            return None; // degenerate
        }
        let arg = (d - big_r * x) / (r * m);
        if !(-1.0..=1.0).contains(&arg) {
            return None;
        }
        let phi = cn.atan2(x);
        Some(phi + branch_sign * arg.acos())
    };

    // Recover (u, v) in torus-surface coords from a 3-D anchor.
    let solve_anchor = |p: Point3| -> Option<(f64, f64)> {
        let dp = p - torus_center;
        let v_axis = dp.dot(&torus_axis);
        let perp_vec = dp - torus_axis * v_axis;
        let u = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
        let tube_r = perp_vec.length();
        // sin v = v_axis / r
        let sin_v = (v_axis / r).clamp(-1.0, 1.0);
        // cos v = (tube_r − R) / r — sign from which side of the major circle
        let cos_v = ((tube_r - big_r) / r).clamp(-1.0, 1.0);
        let v_t = sin_v.atan2(cos_v);
        let v = if v_t < 0.0 {
            v_t + std::f64::consts::TAU
        } else {
            v_t
        };
        Some((u, v))
    };

    let validate_anchor = |p: Point3, side: SurfaceSide| -> Result<(), EdgeTripleError> {
        let off_plane = (p - plane_origin).dot(&n_unit).abs();
        if off_plane > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: off_plane,
            });
        }
        let dp = p - torus_center;
        let v_axis = dp.dot(&torus_axis);
        let perp_vec = dp - torus_axis * v_axis;
        let tube_r = perp_vec.length();
        let off_torus = ((tube_r - big_r).powi(2) + v_axis.powi(2) - r * r).abs();
        if off_torus.sqrt() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: off_torus.sqrt(),
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

    // Branch sign: which v root reproduces the start anchor's v.
    let v_pos = v_at(u_start, 1.0)?;
    let v_neg = v_at(u_start, -1.0)?;
    let canonical = |v: f64| -> f64 {
        let mut x = v;
        while x < 0.0 {
            x += std::f64::consts::TAU;
        }
        while x >= std::f64::consts::TAU {
            x -= std::f64::consts::TAU;
        }
        x
    };
    let dist = |a: f64, b: f64| -> f64 {
        let d = (canonical(a) - canonical(b)).abs();
        d.min(std::f64::consts::TAU - d)
    };
    let branch_sign = if dist(v_pos, v_start_truth) <= dist(v_neg, v_start_truth) {
        1.0
    } else {
        -1.0
    };

    // Sample, build NURBS curve_3d and pcurves.
    let n_samples = 65;
    let mut pts_3d: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_plane: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut pts_torus: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = u_end - u_start;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let u = u_start + frac * span;
        let v = v_at(u, branch_sign)?;
        let cv = v.cos();
        let sv = v.sin();
        let radial = perp1 * u.cos() + perp2 * u.sin();
        let p3 = torus_center + radial * (big_r + r * cv) + torus_axis * (r * sv);
        pts_3d.push(p3);
        // Plane uv.
        let dp_pl = p3 - plane_origin;
        let u_pl = dp_pl.dot(&plane_u) / plane_u_len_sq;
        let v_pl = dp_pl.dot(&plane_v) / plane_v_len_sq;
        pts_plane.push(Point3::new(u_pl, v_pl, 0.0));
        // Torus uv.
        let v_canon = if v < 0.0 {
            v + std::f64::consts::TAU
        } else {
            v
        };
        pts_torus.push(Point3::new(u, v_canon, 0.0));
        params.push(frac);
    }
    unwrap_periodic(&mut pts_torus, n_samples);

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
    let nurbs_plane = match crate::geometry::fitting::fit_curve(&pts_plane, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let nurbs_torus = match crate::geometry::fitting::fit_curve(&pts_torus, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };

    let curve_3d = CurveGeometry::Nurbs(reknot_to(nurbs_3d, u_start, u_end));
    let pcurve_plane = Pcurve::from_nurbs(reknot_to(nurbs_plane, u_start, u_end));
    let pcurve_torus = Pcurve::from_nurbs(reknot_to(nurbs_torus, u_start, u_end));

    let plane_surf = SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        u_axis: plane_u,
        v_axis: plane_v,
    };
    let torus_surf = SurfaceGeometry::Torus {
        center: torus_center,
        axis: torus_axis,
        major_radius,
        minor_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if plane_is_b {
        (pcurve_torus, pcurve_plane, &torus_surf, &plane_surf)
    } else {
        (pcurve_plane, pcurve_torus, &plane_surf, &torus_surf)
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
            kind: AnalyticPair::PlaneTorus,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn unwrap_periodic(pts: &mut [Point3], n: usize) {
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
        let dv = pts[i].y - pts[i - 1].y;
        if dv > std::f64::consts::PI {
            for cp in pts.iter_mut().skip(i) {
                cp.y -= std::f64::consts::TAU;
            }
        } else if dv < -std::f64::consts::PI {
            for cp in pts.iter_mut().skip(i) {
                cp.y += std::f64::consts::TAU;
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
