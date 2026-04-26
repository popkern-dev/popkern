//! Plane × cone analytic dispatch.

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
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Plane × cone when the plane normal is aligned with the cone's
/// axis (axis-perpendicular slice). The intersection is a circle
/// at height `v` along the cone's axis, with radius
/// `v · tan(half_angle)`.
///
/// Returns `None` for oblique slices (which produce ellipses /
/// parabolas / hyperbolas — handled in a later stage with extended
/// dispatch logic).
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cone_perpendicular(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_is_b: bool,
    cone_apex: Point3,
    cone_axis: Vec3,
    half_angle: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let nd = plane_normal.dot(&cone_axis).abs();
    if nd < 1.0 - tol.angular() {
        return None;
    }

    // Height along cone axis where the plane intersects.
    let height = (plane_origin - cone_apex).dot(&plane_normal) / plane_normal.dot(&cone_axis);
    if height <= tol.linear() {
        return None; // at or below apex — degenerate
    }
    let circle_radius = height * half_angle.tan();
    let circle_center = cone_apex + cone_axis * height;

    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(cone_axis);
    let (cone_basis_u, cone_basis_v) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let cone_u_offset = arc_basis_u
        .dot(&cone_basis_v)
        .atan2(arc_basis_u.dot(&cone_basis_u));

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    let theta_of = |p: Point3| -> Option<f64> {
        let d = p - circle_center;
        let cu = d.dot(&arc_basis_u);
        let cv = d.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - circle_radius).abs() > tol.linear() {
            return None;
        }
        Some(cv.atan2(cu))
    };
    let theta_start = match theta_of(start) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::A,
                distance: (start - circle_center).length() - circle_radius,
            }));
        }
    };
    let theta_end = match theta_of(end) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::B,
                distance: (end - circle_center).length() - circle_radius,
            }));
        }
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    // Cone pcurve: u = θ + offset (angular position around axis),
    // v = height (constant for axis-perpendicular slice).
    let pcurve_cone = Pcurve::linear(
        (theta_start + cone_u_offset, height),
        (theta_end + cone_u_offset, height),
        theta_start,
        theta_end,
    );

    // Plane pcurve: project the 3-D circle onto plane UV.
    let center_in_plane = circle_center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    let a_uu = circle_radius * arc_basis_u.dot(&plane_u) / plane_u_len_sq;
    let a_vu = circle_radius * arc_basis_u.dot(&plane_v) / plane_v_len_sq;
    let a_uv = circle_radius * arc_basis_v.dot(&plane_u) / plane_u_len_sq;
    let a_vv = circle_radius * arc_basis_v.dot(&plane_v) / plane_v_len_sq;
    let (u_radius, v_radius, rotation) = decompose_2x2(a_uu, a_uv, a_vu, a_vv);
    let pcurve_plane = Pcurve::elliptical(
        (center_u, center_v),
        u_radius,
        v_radius,
        rotation,
        theta_start,
        theta_end,
        theta_start,
        theta_end,
    );

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: cone_axis,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let plane_surf = SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        u_axis: plane_u,
        v_axis: plane_v,
    };
    let cone_surf = SurfaceGeometry::Cone {
        apex: cone_apex,
        axis: cone_axis,
        half_angle,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if plane_is_b {
        (pcurve_cone, pcurve_plane, &cone_surf, &plane_surf)
    } else {
        (pcurve_plane, pcurve_cone, &plane_surf, &cone_surf)
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
            kind: AnalyticPair::PlaneCone,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
