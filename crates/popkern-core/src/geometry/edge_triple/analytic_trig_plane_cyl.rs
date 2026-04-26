//! Plane × cylinder analytic dispatch (perpendicular slice).

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

/// Perpendicular plane-cylinder slice: plane normal ≈ cylinder
/// axis. Returns `None` for oblique slices so the caller falls back
/// to the existing NURBS path (oblique gives an ellipse — handled
/// in stage 2b.4 with a more general `Elliptical` analytic).
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cylinder_perpendicular(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_is_b: bool,
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Detect perpendicular: |plane_normal · cyl_axis| ≈ 1.
    let nd = plane_normal.dot(&cyl_axis).abs();
    if nd < 1.0 - tol.angular() {
        return None;
    }

    // Centre of the intersection circle: project cyl_origin onto
    // the plane along cyl_axis.
    let dist_along_axis =
        (plane_origin - cyl_origin).dot(&plane_normal) / plane_normal.dot(&cyl_axis);
    let circle_center = cyl_origin + cyl_axis * dist_along_axis;

    // The Arc uses `circle_basis` for 3-D evaluation; the Cylinder
    // surface uses `orthonormal_basis` for its UV. Compute angles
    // in the Arc's frame, then offset the cyl pcurve's u-values by
    // the angular phase between the two frames so
    // `surf_cyl(pcurve(t))` matches `arc(t)` exactly.
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(cyl_axis);
    let (cyl_basis_u, cyl_basis_v) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cyl_axis);
    let cyl_u_offset = arc_basis_u
        .dot(&cyl_basis_v)
        .atan2(arc_basis_u.dot(&cyl_basis_u));

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None, // closed-loop / seam handled in stage 2b.5
    };

    let theta_of = |p: Point3| -> Option<f64> {
        let d = p - circle_center;
        let cu = d.dot(&arc_basis_u);
        let cv = d.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - cyl_radius).abs() > tol.linear() {
            return None;
        }
        Some(cv.atan2(cu))
    };
    let theta_start = match theta_of(start) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::A,
                distance: (start - circle_center).length() - cyl_radius,
            }));
        }
    };
    let theta_end = match theta_of(end) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::B,
                distance: (end - circle_center).length() - cyl_radius,
            }));
        }
    };

    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    let v_on_cyl = (circle_center - cyl_origin).dot(&cyl_axis);
    let pcurve_cyl = Pcurve::linear(
        (theta_start + cyl_u_offset, v_on_cyl),
        (theta_end + cyl_u_offset, v_on_cyl),
        theta_start,
        theta_end,
    );

    let center_in_plane = circle_center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    let a_uu = cyl_radius * arc_basis_u.dot(&plane_u) / plane_u_len_sq;
    let a_vu = cyl_radius * arc_basis_u.dot(&plane_v) / plane_v_len_sq;
    let a_uv = cyl_radius * arc_basis_v.dot(&plane_u) / plane_u_len_sq;
    let a_vv = cyl_radius * arc_basis_v.dot(&plane_v) / plane_v_len_sq;

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

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: cyl_axis,
        radius: cyl_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let plane_surf = SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        u_axis: plane_u,
        v_axis: plane_v,
    };
    let cyl_surf = SurfaceGeometry::Cylinder {
        origin: cyl_origin,
        axis: cyl_axis,
        radius: cyl_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if plane_is_b {
        (pcurve_cyl, pcurve_plane, &cyl_surf, &plane_surf)
    } else {
        (pcurve_plane, pcurve_cyl, &plane_surf, &cyl_surf)
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
            kind: AnalyticPair::PlaneCylinder,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
