//! Cone × torus axial dispatch — octic factors into circles.
//!
//! When the cone apex lies on the torus axis and both axes are
//! parallel (= coincident), the algebraic intersection (a degree-8
//! curve from the degree-2 cone × degree-4 torus) **factors into 0,
//! 1, or 2 horizontal circles** at fixed axial heights, by the same
//! axisymmetric-collapse argument as sphere × torus axial.
//!
//! Math: with cone half-angle `θ`, axial offset `h_axis = (T_c −
//! apex) · a`, the cone constraint `(R + r·cos v) = (h_axis +
//! r·sin v) · tan θ` rearranges to
//!
//!   cos(v + θ) = (sin θ · h_axis − R · cos θ) / r
//!
//! Real solutions iff `|rhs| ≤ 1`. Two `v` branches when in range:
//!
//!   v_± = −θ ± acos(rhs)
//!
//! Each `v_±` defines a circle of radius `R + r·cos v_±` at axial
//! height `r·sin v_±` from the torus centre.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc`.
//! - torus pcurve: `Pcurve::Linear` (constant v, linear u).
//! - cone pcurve: `Pcurve::Linear` (constant v_cone = h_axis +
//!   r·sin v_branch, linear u — both surfaces share the same
//!   `orthonormal_basis(axis)`, so frame phases cancel out).
//! - `ParameterMap::Identity` over θ.

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
pub(super) fn try_cone_torus_axial(
    cone_apex: Point3,
    cone_axis: Vec3,
    half_angle: f64,
    cone_is_b: bool,
    torus_center: Point3,
    torus_axis: Vec3,
    major_radius: f64,
    minor_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detection: axes parallel (allow opposite signs).
    let axis_align = cone_axis.dot(&torus_axis);
    if axis_align.abs() < 1.0 - tol.angular() {
        return None;
    }
    // Cone apex on torus axis.
    let dp = cone_apex - torus_center;
    let perp = dp - torus_axis * dp.dot(&torus_axis);
    if perp.length() > tol.linear() {
        return None;
    }

    // Use torus axis as the canonical direction. If cone axis points
    // opposite, flip the cone half-angle interpretation: cone v >= 0
    // is along cone_axis, which may be -torus_axis.
    let sign = axis_align.signum();
    let h_axis = (torus_center - cone_apex).dot(&torus_axis); // signed
    let big_r = major_radius;
    let r = minor_radius;
    let theta = half_angle;
    let st = theta.sin();
    let ct = theta.cos();
    let rhs = (st * h_axis - big_r * ct) / r;
    if !(-1.0..=1.0).contains(&rhs) {
        return None;
    }
    let arc_acos = rhs.acos();
    // v_± = -θ ± acos. Note: derivation assumed cone axis = torus
    // axis (sign = +1). If sign = -1, the cone surface in the
    // canonical frame uses v_cone = -(h_axis + r·sin v) — pick the
    // branch via the anchor.
    let v_plus = -theta + arc_acos;
    let v_minus = -theta - arc_acos;

    let solve_v = |p: Point3| -> f64 {
        let dp = p - torus_center;
        let v_axis = dp.dot(&torus_axis);
        let perp_vec = dp - torus_axis * v_axis;
        let tube_r = perp_vec.length();
        let sin_v = (v_axis / r).clamp(-1.0, 1.0);
        let cos_v = ((tube_r - big_r) / r).clamp(-1.0, 1.0);
        sin_v.atan2(cos_v)
    };
    let dist_v = |a: f64, b: f64| -> f64 {
        let d = (a - b).abs();
        d.min(std::f64::consts::TAU - d)
    };
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let v_start_truth = solve_v(start);
    let v_branch = if dist_v(v_plus, v_start_truth) <= dist_v(v_minus, v_start_truth) {
        v_plus
    } else {
        v_minus
    };
    let circle_radius = big_r + r * v_branch.cos();
    if circle_radius <= tol.linear() {
        return None;
    }
    let axial_offset_from_torus = r * v_branch.sin();
    let circle_center = torus_center + torus_axis * axial_offset_from_torus;
    // v_cone = signed axial distance from apex along cone_axis. With
    // sign = +1 (cone axis = torus axis): v_cone = h_axis + r·sin v.
    // With sign = -1: v_cone = -(h_axis + r·sin v).
    let v_cone = sign * (h_axis + r * v_branch.sin());
    if v_cone <= tol.linear() {
        return None; // intersection lies on the negative-v half (other nappe)
    }

    let (perp1_t, perp2_t) = crate::geometry::analytic_surfaces::orthonormal_basis(torus_axis);
    let (perp1_c, perp2_c) = crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(torus_axis);
    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let d = p - circle_center;
        let cu = d.dot(&arc_basis_u);
        let cv = d.dot(&arc_basis_v);
        let r_real = (cu * cu + cv * cv).sqrt();
        if (r_real - circle_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r_real - circle_radius,
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
        normal: torus_axis,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let v_branch_canon = if v_branch < 0.0 {
        v_branch + std::f64::consts::TAU
    } else {
        v_branch
    };
    let torus_u_offset = arc_basis_u.dot(&perp2_t).atan2(arc_basis_u.dot(&perp1_t));
    let pcurve_torus = Pcurve::linear(
        (theta_start + torus_u_offset, v_branch_canon),
        (theta_end + torus_u_offset, v_branch_canon),
        theta_start,
        theta_end,
    );
    let cone_u_offset = arc_basis_u.dot(&perp2_c).atan2(arc_basis_u.dot(&perp1_c));
    let pcurve_cone = Pcurve::linear(
        (theta_start + cone_u_offset, v_cone),
        (theta_end + cone_u_offset, v_cone),
        theta_start,
        theta_end,
    );

    let cone_surf = SurfaceGeometry::Cone {
        apex: cone_apex,
        axis: cone_axis,
        half_angle,
    };
    let torus_surf = SurfaceGeometry::Torus {
        center: torus_center,
        axis: torus_axis,
        major_radius,
        minor_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if cone_is_b {
        (pcurve_torus, pcurve_cone, &torus_surf, &cone_surf)
    } else {
        (pcurve_cone, pcurve_torus, &cone_surf, &torus_surf)
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
            kind: AnalyticPair::ConeTorus,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
