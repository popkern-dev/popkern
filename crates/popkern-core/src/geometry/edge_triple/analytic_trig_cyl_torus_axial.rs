//! Cylinder × torus axial dispatch — octic factors into circles.
//!
//! When the cylinder axis coincides with the torus axis, the
//! intersection (degree-2 cyl × degree-4 torus = degree-8 generally)
//! collapses to 0, 1, or 2 horizontal circles by the
//! axisymmetric-collapse argument.
//!
//! Math: a point on the cyl at axial height `v_cyl` (relative to
//! cyl origin) lies on the torus iff
//!
//!   r_cyl = R + r·cos v_torus       (radial constraint)
//!   h_offset + v_cyl = r·sin v_torus (axial constraint)
//!
//! where `h_offset = (cyl_origin − T_centre) · axis`. Combining via
//! `cos² + sin² = 1`:
//!
//!   v_cyl_± = -h_offset ± √(r² − (r_cyl − R)²)
//!
//! Real solutions iff `|r_cyl − R| ≤ r` (the cylinder's radius lies
//! within the tube's swept-radius range `[R−r, R+r]`).
//!
//! Output: exact `CurveGeometry::Arc` + `Pcurve::Linear` on both
//! surfaces (constant v on each, linear u).

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
use popkern_math::PortableMath;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cyl_torus_axial(
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
    cyl_is_b: bool,
    torus_center: Point3,
    torus_axis: Vec3,
    major_radius: f64,
    minor_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detection: cyl axis ∥ torus axis.
    let axis_align = cyl_axis.dot(&torus_axis);
    if axis_align.abs() < 1.0 - tol.angular() {
        return None;
    }
    // Cyl origin on torus axis.
    let dp = cyl_origin - torus_center;
    let perp = dp - torus_axis * dp.dot(&torus_axis);
    if perp.length() > tol.linear() {
        return None;
    }

    let big_r = major_radius;
    let r = minor_radius;
    let radial_diff = cyl_radius - big_r;
    let radical = r * r - radial_diff * radial_diff;
    if radical < -tol.linear() {
        return None; // cyl radius outside [R−r, R+r] — no intersection
    }
    let sqrt_radical = radical.max(0.0).sqrt();
    let h_offset = (cyl_origin - torus_center).dot(&torus_axis);
    // v_cyl_± = ±sqrt_radical − h_offset (in cyl-axis frame).
    let v_cyl_plus = -h_offset + sqrt_radical;
    let v_cyl_minus = -h_offset - sqrt_radical;

    // Pick branch by start-anchor's cyl-v.
    let cyl_v_of = |p: Point3| -> f64 { (p - cyl_origin).dot(&cyl_axis) };
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let v_start = cyl_v_of(start);
    let v_branch = if (v_cyl_plus - v_start).abs() <= (v_cyl_minus - v_start).abs() {
        v_cyl_plus
    } else {
        v_cyl_minus
    };

    let circle_radius = cyl_radius;
    let circle_center = cyl_origin + cyl_axis * v_branch;
    // v_torus from sin and cos (sin v_torus = (h_offset + v_branch)/r,
    // cos v_torus = (r_cyl − R)/r).
    let sin_v_torus = ((h_offset + v_branch) / r).clamp(-1.0, 1.0);
    let cos_v_torus = (radial_diff / r).clamp(-1.0, 1.0);
    let v_torus = sin_v_torus.patan2(cos_v_torus);
    let v_torus_canon = if v_torus < 0.0 {
        v_torus + std::f64::consts::TAU
    } else {
        v_torus
    };

    let (perp1_c, perp2_c) = crate::geometry::analytic_surfaces::orthonormal_basis(cyl_axis);
    let (perp1_t, perp2_t) = crate::geometry::analytic_surfaces::orthonormal_basis(torus_axis);
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(cyl_axis);
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
        Ok(cv.patan2(cu))
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
        normal: cyl_axis,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let cyl_u_offset = arc_basis_u.dot(&perp2_c).patan2(arc_basis_u.dot(&perp1_c));
    let pcurve_cyl = Pcurve::linear(
        (theta_start + cyl_u_offset, v_branch),
        (theta_end + cyl_u_offset, v_branch),
        theta_start,
        theta_end,
    );
    let torus_u_offset = arc_basis_u.dot(&perp2_t).patan2(arc_basis_u.dot(&perp1_t));
    let pcurve_torus = Pcurve::linear(
        (theta_start + torus_u_offset, v_torus_canon),
        (theta_end + torus_u_offset, v_torus_canon),
        theta_start,
        theta_end,
    );

    let cyl_surf = SurfaceGeometry::Cylinder {
        origin: cyl_origin,
        axis: cyl_axis,
        radius: cyl_radius,
    };
    let torus_surf = SurfaceGeometry::Torus {
        center: torus_center,
        axis: torus_axis,
        major_radius,
        minor_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if cyl_is_b {
        (pcurve_torus, pcurve_cyl, &torus_surf, &cyl_surf)
    } else {
        (pcurve_cyl, pcurve_torus, &cyl_surf, &torus_surf)
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
            kind: AnalyticPair::CylinderTorus,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
