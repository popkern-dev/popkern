//! Trig-pair analytic dispatch (plane × cylinder, plane × sphere,
//! plane × cone). Produces analytic primitives for `curve_3d` and
//! both pcurves with `ParameterMap::Identity` — matches commercial-kernel
//! and OpenCASCADE practice for the canonical analytic-pair
//! cases. See ADR §13.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::analytic_trig_cone::try_plane_cone_perpendicular;
use crate::geometry::edge_triple::analytic_trig_cyl_sphere::try_cyl_sphere_axial;
use crate::geometry::edge_triple::analytic_trig_helpers::decompose_2x2;
use crate::geometry::edge_triple::analytic_trig_sphere_sphere::try_sphere_sphere_polar;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try the analytic trig-pair dispatch.
///
/// Returns `Some(Ok(triple))` on success, `Some(Err(...))` on
/// real failure (anchors off circle, etc.), or `None` if the pair
/// isn't a trig pair we handle here (caller should fall through to
/// the NURBS-via-`intersect_surfaces` path).
pub(super) fn try_dispatch_trig(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    use SurfaceGeometry::*;
    match (surf_a, surf_b) {
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
        ) => {
            try_plane_cylinder_perpendicular(*po, *pn, *pu, *pv, false, *co, *ca, *cr, anchors, tol)
        }
        (
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => {
            try_plane_cylinder_perpendicular(*po, *pn, *pu, *pv, true, *co, *ca, *cr, anchors, tol)
        }
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Sphere { center, radius },
        ) => try_plane_sphere_polar(*po, *pn, *pu, *pv, false, *center, *radius, anchors, tol),
        (
            Sphere { center, radius },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => try_plane_sphere_polar(*po, *pn, *pu, *pv, true, *center, *radius, anchors, tol),
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Cone {
                apex,
                axis,
                half_angle,
            },
        ) => try_plane_cone_perpendicular(
            *po,
            *pn,
            *pu,
            *pv,
            false,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        ),
        (
            Cone {
                apex,
                axis,
                half_angle,
            },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => try_plane_cone_perpendicular(
            *po,
            *pn,
            *pu,
            *pv,
            true,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        ),
        (
            Sphere {
                center: ca,
                radius: ra,
            },
            Sphere {
                center: cb,
                radius: rb,
            },
        ) => try_sphere_sphere_polar(*ca, *ra, *cb, *rb, anchors, tol),
        (
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
            Sphere { center, radius },
        ) => try_cyl_sphere_axial(*co, *ca, *cr, false, *center, *radius, anchors, hint, tol),
        (
            Sphere { center, radius },
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
        ) => try_cyl_sphere_axial(*co, *ca, *cr, true, *center, *radius, anchors, hint, tol),
        _ => None,
    }
}

/// Perpendicular plane-cylinder slice: plane normal ≈ cylinder
/// axis. Returns `None` for oblique slices so the caller falls back
/// to the existing NURBS path (oblique gives an ellipse — handled
/// in stage 2b.4 with a more general `Elliptical` analytic).
#[allow(clippy::too_many_arguments)]
fn try_plane_cylinder_perpendicular(
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
    // arc_basis_u expressed in cyl basis = (arc_u · cyl_u, arc_u · cyl_v).
    // The angle of arc_basis_u in cyl frame = atan2(.cyl_v_comp, .cyl_u_comp).
    // surf_cyl(u_offset, _) = cyl_basis_u·cos(u_offset) + cyl_basis_v·sin(u_offset),
    // so to land on arc_basis_u we need u_offset = that angle.
    let cyl_u_offset = arc_basis_u
        .dot(&cyl_basis_v)
        .atan2(arc_basis_u.dot(&cyl_basis_u));

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None, // closed-loop / seam handled in stage 2b.5
    };

    // Compute angles θ_start, θ_end IN THE ARC'S FRAME so that
    // `Arc::point_at(θ) ≈ anchor`.
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

    // Choose CCW arc by default; flip end angle by ±2π if needed
    // to get a single sweep from start to end without wrapping
    // through the seam unnecessarily. Keep direction so
    // `theta_start < theta_end`.
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    // Cyl pcurve: linear in (u, v) space, where u = arc_θ +
    // cyl_u_offset (frame phase) and v = projection along cyl_axis.
    // The triple's parameter is the arc's θ; the pcurve's stored
    // u-coordinate adds the offset so that surf_cyl evaluates to
    // the same 3-D point as the Arc.
    let v_on_cyl = (circle_center - cyl_origin).dot(&cyl_axis);
    let pcurve_cyl = Pcurve::linear(
        (theta_start + cyl_u_offset, v_on_cyl),
        (theta_end + cyl_u_offset, v_on_cyl),
        theta_start,
        theta_end,
    );

    // Plane pcurve: project the 3-D circle onto plane UV. The
    // resulting curve is an ellipse parameterised by θ.
    let center_in_plane = circle_center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    // (cos θ basis · pu, cos θ basis · pv) and likewise for sin θ.
    let a_uu = cyl_radius * arc_basis_u.dot(&plane_u) / plane_u_len_sq;
    let a_vu = cyl_radius * arc_basis_u.dot(&plane_v) / plane_v_len_sq;
    let a_uv = cyl_radius * arc_basis_v.dot(&plane_u) / plane_u_len_sq;
    let a_vv = cyl_radius * arc_basis_v.dot(&plane_v) / plane_v_len_sq;

    // The mapping (cos θ, sin θ) → (u, v) is the linear map A.
    // Decompose into A = R(rot) · diag(u_radius, v_radius) · R(0)
    // — a simple ellipse with rotation. SVD-style decomposition.
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

/// Plane × sphere when the plane normal is parallel to the sphere's
/// UV polar axis `(0, 0, 1)` (popkern's sphere parameterisation).
/// Returns `None` for non-axis-aligned slices; those need a great-
/// circle pcurve in sphere UV which a Linear variant cannot express
/// — handled in stage 2b.5 via a sampled Nurbs pcurve fallback.
#[allow(clippy::too_many_arguments)]
fn try_plane_sphere_polar(
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
    // Sphere u-axis = global X for points on the equator-like circle.
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

    // Plane pcurve: same elliptical projection used in plane×cyl.
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

// `decompose_2x2` lives in `analytic_trig_helpers.rs` so this file
// stays under the 500-line cap.
