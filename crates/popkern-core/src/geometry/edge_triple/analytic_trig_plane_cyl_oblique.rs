//! Plane × cylinder analytic dispatch (oblique slice).
//!
//! When a plane cuts a cylinder at any angle other than perpendicular
//! or parallel-to-axis, the intersection is a 3-D ellipse (proven by
//! the Dandelin spheres construction). This dispatcher returns:
//!
//! - `curve_3d` as `CurveGeometry::Ellipse` with non-orthogonal
//!   coplanar `(u_basis, v_basis)` derived from the cylinder's
//!   `orthonormal_basis(axis)` projected through the cutting plane.
//! - `pcurve` on the cylinder as `Pcurve::sinusoidal` — `u(s) = s`
//!   (the cylinder's natural angular parameter), `v(s)` is the
//!   axial position which traces a single sine wave per orbit.
//! - `pcurve` on the plane as `Pcurve::elliptical` — both `u` and
//!   `v` are linear-in-(cos s, sin s), giving an exact ellipse fit.
//!
//! `ParameterMap::Identity` holds: the `Ellipse` 3-D curve uses the
//! same basis vectors as the cylinder surface, so `s` is the shared
//! parameter on both sides.

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

/// Oblique plane-cylinder slice. Returns `None` when the plane is
/// either perpendicular to the cylinder axis (handled by the
/// perpendicular dispatcher) or parallel to it (gives 0/1/2 lines —
/// distinct geometry, separate dispatcher).
#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_cylinder_oblique(
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
    // Reject the perpendicular case (handled elsewhere) and the
    // parallel case (no closed-curve intersection).
    let cn = plane_normal.dot(&cyl_axis);
    let cn_abs = cn.abs();
    if cn_abs >= 1.0 - tol.angular() {
        return None; // perpendicular
    }
    if cn_abs <= tol.angular() {
        return None; // parallel — distinct geometry
    }

    // Cylinder's orthonormal basis (perpendicular to axis).
    let (eu, ev) = crate::geometry::analytic_surfaces::orthonormal_basis(cyl_axis);

    // Centre of the intersection ellipse: the cylinder axis pierces
    // the plane at this point.
    let d = (plane_origin - cyl_origin).dot(&plane_normal);
    let center = cyl_origin + cyl_axis * (d / cn);

    // Ellipse semi-axis vectors. Both lie in the cutting plane
    // (their components along `plane_normal` are zero by
    // construction). They are equal-length pairs of `r` only in the
    // perpendicular limit; otherwise they pick up an axis-parallel
    // tilt that produces the 3-D ellipse. Generally non-orthogonal
    // when the plane normal is not aligned with `eu` or `ev`.
    let eu_dot_n = eu.dot(&plane_normal);
    let ev_dot_n = ev.dot(&plane_normal);
    let u_basis = (eu - cyl_axis * (eu_dot_n / cn)) * cyl_radius;
    let v_basis = (ev - cyl_axis * (ev_dot_n / cn)) * cyl_radius;

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None, // closed-loop / seam handled later
    };

    // Solve for parameter `s` at each anchor. The cylinder surface
    // parameterisation gives `(anchor - center) · eu = r cos s` and
    // `(anchor - center) · ev = r sin s`, so `s = atan2(d·ev, d·eu)`.
    let s_of = |p: Point3| -> Option<f64> {
        let dp = p - center;
        let cu = dp.dot(&eu);
        let cv = dp.dot(&ev);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - cyl_radius).abs() > tol.linear() {
            return None;
        }
        Some(cv.atan2(cu))
    };
    let s_start = match s_of(start) {
        Some(s) => s,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::A,
                distance: ((start - center).dot(&eu).powi(2) + (start - center).dot(&ev).powi(2))
                    .sqrt()
                    - cyl_radius,
            }));
        }
    };
    let s_end = match s_of(end) {
        Some(s) => s,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::B,
                distance: ((end - center).dot(&eu).powi(2) + (end - center).dot(&ev).powi(2))
                    .sqrt()
                    - cyl_radius,
            }));
        }
    };
    let s_end = if s_end < s_start {
        s_end + std::f64::consts::TAU
    } else {
        s_end
    };

    // Cylinder pcurve: u(s) = s, v(s) is the axial position which
    // varies sinusoidally because the cutting plane is tilted off
    // the cylinder axis.
    //
    // From `surf_cyl(u, v) = co + r cos u eu + r sin u ev + v a`,
    // setting `u = s` and requiring the surface point to lie on the
    // cutting plane gives:
    //   v(s) = d/cn − (r·(eu·n)/cn) cos s − (r·(ev·n)/cn) sin s
    let pcurve_cyl = Pcurve::sinusoidal(
        0.0,
        d / cn,
        -cyl_radius * eu_dot_n / cn,
        -cyl_radius * ev_dot_n / cn,
        s_start,
        s_end,
    );

    // Plane pcurve: project `p(s) − po = (center − po) + cos s · U +
    // sin s · V` onto the plane's `(u_axis, v_axis)`. Both u and v
    // are of the form `A + B cos s + C sin s`, an exact ellipse in
    // plane UV — the existing `Elliptical` variant handles this via
    // the same 2×2 SVD decomposition used in the perpendicular case.
    let center_in_plane = center - plane_origin;
    let plane_u_len_sq = plane_u.dot(&plane_u);
    let plane_v_len_sq = plane_v.dot(&plane_v);
    let center_u = center_in_plane.dot(&plane_u) / plane_u_len_sq;
    let center_v = center_in_plane.dot(&plane_v) / plane_v_len_sq;
    let a_uu = u_basis.dot(&plane_u) / plane_u_len_sq;
    let a_vu = u_basis.dot(&plane_v) / plane_v_len_sq;
    let a_uv = v_basis.dot(&plane_u) / plane_u_len_sq;
    let a_vv = v_basis.dot(&plane_v) / plane_v_len_sq;
    let (u_radius, v_radius, rotation) = decompose_2x2(a_uu, a_uv, a_vu, a_vv);
    let pcurve_plane = Pcurve::elliptical(
        (center_u, center_v),
        u_radius,
        v_radius,
        rotation,
        s_start,
        s_end,
        s_start,
        s_end,
    );

    let curve_3d = CurveGeometry::Ellipse {
        center,
        u_basis,
        v_basis,
        start_angle: s_start,
        end_angle: s_end,
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
        t_start: s_start,
        t_end: s_end,
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
