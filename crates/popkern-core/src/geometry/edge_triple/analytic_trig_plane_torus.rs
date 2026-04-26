//! Plane × torus axis-perpendicular dispatch.
//!
//! When the cutting plane is perpendicular to the torus axis, the
//! intersection is 0, 1, or 2 concentric circles depending on the
//! signed distance `h` between torus centre and plane:
//!
//! - `|h| > minor_radius`: empty (plane misses the tube).
//! - `|h| = minor_radius`: one tangent circle.
//! - `|h| < minor_radius`: two circles at radii `R ± r·cos v_branch`,
//!   where `v_branch ∈ {arcsin(h/r), π − arcsin(h/r)}`.
//!
//! This dispatcher selects whichever branch contains the anchors and
//! emits an exact `CurveGeometry::Arc` plus `Pcurve::Linear` on the
//! torus side (constant v, linear u) and `Pcurve::Elliptical` on the
//! plane side. `ParameterMap::Identity` holds via the shared u angle.

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

#[allow(clippy::too_many_arguments)]
pub(super) fn try_plane_torus_perpendicular(
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
    // Require perpendicular: plane normal parallel to torus axis.
    let nd = plane_normal.dot(&torus_axis).abs();
    if nd < 1.0 - tol.angular() {
        return None;
    }

    // Signed distance from torus centre to plane, along torus axis.
    let h = (plane_origin - torus_center).dot(&plane_normal) / plane_normal.dot(&torus_axis);
    if h.abs() >= minor_radius - tol.linear() {
        return None; // tangent or non-intersecting
    }

    // Two branches: v₁ = arcsin(h/r), v₂ = π − v₁. Same axial height,
    // different tube radius (R ± |r·cos v|).
    let sin_v = (h / minor_radius).clamp(-1.0, 1.0);
    let v_branch_inner = sin_v.asin(); // |cos v| · r  pointing toward the axis side
    let v_branch_outer = std::f64::consts::PI - v_branch_inner;

    // Pick the branch whose tube radius matches the start anchor.
    let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(torus_axis);
    let circle_center = torus_center + torus_axis * h;

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    let radial_at = |p: Point3| -> f64 {
        let dp = p - torus_center;
        let proj = dp - torus_axis * dp.dot(&torus_axis);
        proj.length()
    };
    let r_at_start = radial_at(start);
    let r_outer = major_radius + minor_radius * v_branch_outer.cos();
    let r_inner = major_radius + minor_radius * v_branch_inner.cos();
    // Choose the closer branch.
    let (v_branch, circle_radius) = if (r_at_start - r_outer).abs() <= (r_at_start - r_inner).abs()
    {
        (v_branch_outer, r_outer)
    } else {
        (v_branch_inner, r_inner)
    };
    if circle_radius <= tol.linear() {
        return None; // degenerate inner branch
    }

    // Frame-phase offset: Arc uses `circle_basis(axis)`, torus uses
    // `orthonormal_basis(axis)`. Bake the offset into the Linear
    // torus pcurve's u-coords so `surf_torus(pcurve_u(t), v_branch) ≡
    // arc(t)`.
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(torus_axis);
    let torus_u_offset = arc_basis_u.dot(&perp2).atan2(arc_basis_u.dot(&perp1));

    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - circle_center;
        let cu = dp.dot(&arc_basis_u);
        let cv = dp.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - circle_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r - circle_radius,
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

    // Torus pcurve: linear u (with frame-phase offset), constant v.
    let pcurve_torus = Pcurve::linear(
        (theta_start + torus_u_offset, v_branch),
        (theta_end + torus_u_offset, v_branch),
        theta_start,
        theta_end,
    );

    // Plane pcurve: project the circle onto the cutting plane via
    // 2×2 SVD. Same machinery as plane×cyl perpendicular.
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
