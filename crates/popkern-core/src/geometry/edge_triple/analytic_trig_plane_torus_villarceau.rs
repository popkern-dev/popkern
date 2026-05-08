//! Plane × torus Villarceau circles dispatch.
//!
//! When the cutting plane is **bitangent** to the torus — passing
//! through the torus centre and tilted by `α = arcsin(r/R)` from the
//! torus's symmetry plane — the quartic plane-torus intersection
//! factors into **two circles** of radius `R` (the major radius).
//! These are the Villarceau circles, named for Yvon Villarceau (1848).
//!
//! Detection criteria:
//! 1. Plane passes through torus centre (`|signed_dist| ≤ tol.linear()`).
//! 2. `|n · axis| ≈ r / R` (within `tol.angular()`).
//!
//! Geometry of the two circles (derived in plane-uv coordinates
//! `(a, b)` with basis `e_y = a × tilt_dir` and `e_in_plane =
//! cos α · perp1_in_symmetry_plane − sin α · axis`):
//!
//! - Both circles have radius `R`.
//! - Circle 1 centre in plane-uv: `(a, b) = (r, 0)` ⇒ 3-D centre =
//!   `T_centre + r · e_y`.
//! - Circle 2 centre in plane-uv: `(a, b) = (−r, 0)` ⇒ 3-D centre =
//!   `T_centre − r · e_y`.
//! - Both circles lie in the cutting plane and have plane-normal as
//!   their circle normal.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc`.
//! - plane pcurve: exact `Pcurve::Elliptical` (circle projected onto
//!   plane uv).
//! - torus pcurve: degree-5 NURBS sample-and-fit. The (u_torus,
//!   v_torus) trace is `atan2` of rational sinusoids.
//! - `ParameterMap::Identity` via shared θ.

use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::analytic_trig_helpers::decompose_2x2;
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
pub(super) fn try_plane_torus_villarceau(
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
    // Detection: plane through torus centre.
    let n_len_sq = plane_normal.dot(&plane_normal);
    if n_len_sq < tol.solver() {
        return None;
    }
    let n_unit = plane_normal * (1.0 / n_len_sq.sqrt());
    let signed_dist = (torus_center - plane_origin).dot(&n_unit);
    if signed_dist.abs() > tol.linear() {
        return None;
    }
    // Detection: the Villarceau plane is tilted by α from the
    // torus's symmetry plane, with `sin α = r / R`. The plane
    // normal therefore makes angle α with the torus axis, giving
    // `|n · axis| = cos α = √(R² − r²) / R`. Equivalently:
    // `sin α = √(1 − (n·axis)²) = r / R`.
    let cn = n_unit.dot(&torus_axis).abs();
    let sin_a_actual = (1.0 - cn * cn).max(0.0).sqrt();
    let sin_a_target = minor_radius / major_radius;
    if (sin_a_actual - sin_a_target).abs() > tol.angular() {
        return None;
    }

    // Tilt direction: projection of plane normal onto symmetry plane,
    // normalised. e_y = axis × tilt_dir is the in-plane axis along
    // which the two Villarceau-circle centres are offset by ±r.
    let tilt_proj = n_unit - torus_axis * n_unit.dot(&torus_axis);
    let tilt_len = tilt_proj.length();
    if tilt_len < tol.solver() {
        return None; // perpendicular plane — handled elsewhere
    }
    let tilt_dir = tilt_proj * (1.0 / tilt_len);
    let e_y = torus_axis
        .cross(&tilt_dir)
        .normalized()
        .unwrap_or(Vec3::zero());
    if e_y.length() < tol.solver() {
        return None;
    }

    // Two Villarceau-circle centres.
    let c1 = torus_center + e_y * minor_radius;
    let c2 = torus_center - e_y * minor_radius;

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    // Pick the branch the start anchor lies on.
    let d1 = (start - c1).length();
    let d2 = (start - c2).length();
    let circle_center = if (d1 - major_radius).abs() <= (d2 - major_radius).abs() {
        c1
    } else {
        c2
    };

    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(n_unit);

    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - circle_center;
        let cu = dp.dot(&arc_basis_u);
        let cv = dp.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - major_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r - major_radius,
            });
        }
        Ok(cv.atan2(cu))
    };
    let theta_start = match theta_of(start, SurfaceSide::A) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = match theta_of(end, SurfaceSide::B) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: n_unit,
        radius: major_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    // Plane pcurve: Elliptical via 2×2 SVD. Same machinery as plane×cyl
    // perpendicular.
    let center_in_plane = circle_center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    let a_uu = major_radius * arc_basis_u.dot(&plane_u) / plane_u_len_sq;
    let a_vu = major_radius * arc_basis_u.dot(&plane_v) / plane_v_len_sq;
    let a_uv = major_radius * arc_basis_v.dot(&plane_u) / plane_u_len_sq;
    let a_vv = major_radius * arc_basis_v.dot(&plane_v) / plane_v_len_sq;
    let (u_radius, v_radius, rotation, phase) = decompose_2x2(a_uu, a_uv, a_vu, a_vv);
    let pcurve_plane = Pcurve::elliptical(
        (center_u, center_v),
        u_radius,
        v_radius,
        rotation,
        theta_start - phase,
        theta_end - phase,
        theta_start,
        theta_end,
    );

    // Torus pcurve: NURBS sample-and-fit in (u_torus, v_torus). The
    // trace is non-trivial — atan2 of rational sinusoids in θ.
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(torus_axis);
    let n_samples = 65;
    let mut torus_pts: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let span = theta_end - theta_start;
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let theta = theta_start + frac * span;
        let p3 = circle_center
            + arc_basis_u * (major_radius * theta.cos())
            + arc_basis_v * (major_radius * theta.sin());
        let dp = p3 - torus_center;
        let v_axis = dp.dot(&torus_axis);
        let perp_vec = dp - torus_axis * v_axis;
        // u_torus = atan2 of perp coords.
        let u_t = perp_vec.dot(&perp2).atan2(perp_vec.dot(&perp1));
        // v_torus from sin v = v_axis / r.
        let sin_v = (v_axis / minor_radius).clamp(-1.0, 1.0);
        // cos v has sign determined by tube radius vs. major: if
        // tube_radius > R then cos v > 0 (outer tube); else < 0.
        let tube_radius = perp_vec.length();
        let cos_v = if tube_radius >= major_radius {
            ((tube_radius - major_radius) / minor_radius).clamp(-1.0, 1.0)
        } else {
            -((major_radius - tube_radius) / minor_radius).clamp(-1.0, 1.0)
        };
        let v_t = sin_v.atan2(cos_v);
        let v_t = if v_t < 0.0 {
            v_t + std::f64::consts::TAU
        } else {
            v_t
        };
        torus_pts.push(Point3::new(u_t, v_t, 0.0));
        params.push(frac);
    }
    // Unwrap u_torus and v_torus jumps for clean NURBS fit.
    unwrap_periodic(&mut torus_pts, n_samples);

    let torus_nurbs = match crate::geometry::fitting::fit_curve(&torus_pts, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => return None,
    };
    let pcurve_torus = Pcurve::from_nurbs(reknot_to(torus_nurbs, theta_start, theta_end));

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
        t_start: theta_start,
        t_end: theta_end,
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
