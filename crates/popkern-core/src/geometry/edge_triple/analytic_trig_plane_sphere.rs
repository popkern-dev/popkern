//! Plane × sphere analytic dispatch (polar-axis-aligned slice).

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

/// Plane × sphere when the plane normal is parallel to the sphere's
/// UV polar axis `(0, 0, 1)` (popkern's sphere parameterisation).
/// Returns `None` for non-axis-aligned slices; those need a great-
/// circle pcurve in sphere UV which a Linear variant cannot express
/// — handled in stage 2b.5 via a sampled Nurbs pcurve fallback.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_sphere_polar(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_is_b: bool,
    sphere_center: Point3,
    sphere_radius: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // popkern's sphere UV uses (0,0,1) as polar axis. Detect alignment.
    let polar = Vec3::new(0.0, 0.0, 1.0);
    if plane_normal.dot(&polar).abs() < 1.0 - tol.angular() {
        return None;
    }
    // Distance from sphere centre to the plane (signed along polar).
    let h = (plane_origin - sphere_center).dot(&plane_normal) / plane_normal.dot(&polar);
    if h.abs() >= sphere_radius - tol.linear() {
        return None; // tangent or non-intersecting
    }
    let circle_radius = (sphere_radius * sphere_radius - h * h).sqrt();
    let circle_center = sphere_center + polar * h;

    // Arc uses circle_basis(polar). Sphere u-coordinate corresponds
    // to atan2(p.y - cy, p.x - cx) which is the angle in the global
    // XY plane (NOT in any chosen basis). Compute offset between
    // arc_basis and global X/Y so the Identity map holds.
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(polar);
    let global_x = Vec3::new(1.0, 0.0, 0.0);
    let global_y = Vec3::new(0.0, 1.0, 0.0);
    let sphere_u_offset = arc_basis_u.dot(&global_y).atan2(arc_basis_u.dot(&global_x));

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

    // Sphere pcurve: u = θ + offset (longitude), v = colatitude
    // = acos(h/r_sphere). Linear in (u, v) since v is constant.
    let v_sphere = (h / sphere_radius).clamp(-1.0, 1.0).acos();
    let pcurve_sphere = Pcurve::linear(
        (theta_start + sphere_u_offset, v_sphere),
        (theta_end + sphere_u_offset, v_sphere),
        theta_start,
        theta_end,
    );

    let center_in_plane = circle_center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    let a_uu = circle_radius * arc_basis_u.dot(&plane_u) / plane_u_len_sq;
    let a_vu = circle_radius * arc_basis_u.dot(&plane_v) / plane_v_len_sq;
    let a_uv = circle_radius * arc_basis_v.dot(&plane_u) / plane_u_len_sq;
    let a_vv = circle_radius * arc_basis_v.dot(&plane_v) / plane_v_len_sq;
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
        normal: polar,
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
    let sphere_surf = SurfaceGeometry::Sphere {
        center: sphere_center,
        radius: sphere_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if plane_is_b {
        (pcurve_sphere, pcurve_plane, &sphere_surf, &plane_surf)
    } else {
        (pcurve_plane, pcurve_sphere, &plane_surf, &sphere_surf)
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
            kind: AnalyticPair::PlaneSphere,
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
