//! Cylinder × cone coaxial dispatch.
//!
//! When the cylinder and cone share an axis (parallel directions and
//! the cone's apex lies on the cylinder axis), the surfaces meet in
//! exactly one circle: at the height where the cone's radius
//! `v_cone · tan(half_angle)` equals the cylinder's constant radius.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc`.
//! - cyl pcurve: `Pcurve::Linear` (linear u, constant v at the
//!   intersection height projected onto cyl axis).
//! - cone pcurve: `Pcurve::Linear` (linear u, constant v at the
//!   axial height where the cone's radius equals the cyl's).
//! - `ParameterMap::Identity` via shared angular u.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_cyl_cone_coaxial(
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
    // Require axes parallel.
    let parallel = cyl_axis.dot(&cone_axis).abs();
    if parallel < 1.0 - tol.angular() {
        return None;
    }
    // Require cone apex on the cyl axis.
    let dp = cone_apex - cyl_origin;
    let along = dp.dot(&cyl_axis);
    let perp = dp - cyl_axis * along;
    if perp.length() > tol.linear() {
        return None;
    }

    // Intersection circle radius = cyl_radius. Cone v at this radius
    // satisfies v · tan(half_angle) = cyl_radius.
    let t = half_angle.tan();
    if t.abs() < tol.solver() {
        return None; // degenerate cone (zero half-angle)
    }
    let v_cone = cyl_radius / t;
    // 3-D circle centre: apex + v_cone · cone_axis. Check direction
    // of cone_axis matches cyl_axis (positive sense).
    let axis_sign = cone_axis.dot(&cyl_axis).signum();
    let cone_axis_eff = cone_axis * axis_sign;
    let v_cone_eff = v_cone * axis_sign;
    let circle_center = cone_apex + cone_axis_eff * v_cone;
    // v on cylinder: project circle_center − cyl_origin onto cyl_axis.
    let v_cyl = (circle_center - cyl_origin).dot(&cyl_axis);

    let (perp_a, perp_b) = crate::geometry::analytic_surfaces::orthonormal_basis(cyl_axis);
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(cyl_axis);
    let cyl_u_offset = arc_basis_u.dot(&perp_b).atan2(arc_basis_u.dot(&perp_a));
    // Cone uses orthonormal_basis(cone_axis); since cone_axis_eff is
    // parallel to cyl_axis (same or opposite sign), the cone's
    // perpendicular basis MAY differ. Compute its frame-phase offset
    // separately.
    let (cone_perp_a, cone_perp_b) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let cone_u_offset = arc_basis_u
        .dot(&cone_perp_b)
        .atan2(arc_basis_u.dot(&cone_perp_a));

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - circle_center;
        let cu = dp.dot(&arc_basis_u);
        let cv = dp.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - cyl_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r - cyl_radius,
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
        normal: cyl_axis,
        radius: cyl_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    // Cylinder pcurve: u = θ + frame-phase offset, v = constant.
    let pcurve_cyl = Pcurve::linear(
        (theta_start + cyl_u_offset, v_cyl),
        (theta_end + cyl_u_offset, v_cyl),
        theta_start,
        theta_end,
    );
    // Cone pcurve: u = θ + cone-frame offset, v = v_cone (constant).
    // If axis_sign is negative (cone facing opposite direction),
    // cone v = -v_cone_eff so the cone evaluation lands at the same
    // 3-D height. But popkern's cone is defined for v ≥ 0; if
    // axis_sign is negative we have a tip-down cone whose v-coord
    // for the intersection is negative — out of the canonical
    // single-nappe model. Reject for safety.
    if v_cone_eff <= tol.linear() {
        return None;
    }
    let pcurve_cone = Pcurve::linear(
        (theta_start + cone_u_offset, v_cone),
        (theta_end + cone_u_offset, v_cone),
        theta_start,
        theta_end,
    );

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
        t_start: theta_start,
        t_end: theta_end,
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
