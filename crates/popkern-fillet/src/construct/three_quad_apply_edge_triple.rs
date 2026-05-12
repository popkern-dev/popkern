//! EdgeTriple integration for cap-outer-boundary edges (ADR
//! `2026-04-25-edge-triple-factory.md`, stage 6).
//!
//! `try_build_edge_triple_for_seam` is called from
//! [`super::three_quad_apply::reproject_loop_pcurves_onto_gregory`]
//! (gated on `POPKERN_EDGE_TRIPLE=1`) to construct a coherent
//! curve_3d + both pcurves for cap-vs-analytic-face seam edges.
//! Closes the chord-vs-arc disagreement diagnosed in Task #144.
//!
//! Returns `None` when:
//! - The edge has no 3-D curve (internal cap-to-cap edges from the
//!   poke; reproject_loop falls through to the legacy
//!   `OldSurfPcurve` path).
//! - The other-side surface is non-analytic (NURBS / Gregory).
//!   Stage 7 will add the marching catch-all for those.
//! - The edge midpoint doesn't project to a Gregory UV box
//!   boundary within tolerance (i.e. the edge isn't on a
//!   canonical boundary; this can happen on some degenerate
//!   poke configurations).
//! - Anchor verification fails (cylinder fillet's apex-end vertex
//!   isn't where the Gregory's boundary endpoint is — Bug A
//!   territory; needs the underlying chord-line bug fixed first).

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::edge_triple::{
    BoundarySide, EdgeAnchors, EdgeHint, EdgeTripleOptions, SurfaceSide, compute_edge_geometry,
};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Tolerance for "approximately on a UV-box boundary axis" when
/// classifying which boundary an edge lies on. Gregory UV box is
/// `[0, 1]²` so 0.05 = 5% of axis range. Below this threshold,
/// the projection is treated as on-boundary.
const BOUNDARY_TOL_UV: f64 = 0.05;

/// Try to construct a coherent edge-triple for the half-edge `he`
/// on the cap face, where the other side of the edge is an
/// analytic surface. On success, returns `(curve_3d, t_start,
/// t_end, cap_pcurve, other_pcurve)` ready to write into the BRep.
pub(super) fn try_build_edge_triple_for_seam(
    brep: &BRep,
    he: HalfEdgeHandle,
    cap_face: FaceHandle,
    cap_surf: &SurfaceGeometry,
) -> Option<(CurveGeometry, f64, f64, NurbsCurve, NurbsCurve)> {
    use popkern_core::geometry::curve::CurveEval;

    // Cap surface must be a Gregory boundary or an analytic sphere cap.
    if !matches!(
        cap_surf,
        SurfaceGeometry::Gregory(_) | SurfaceGeometry::Sphere { .. }
    ) {
        return None;
    }

    // Identify the other-side face's surface. It must be analytic
    // (Plane/Cylinder/Sphere/Torus) for stage-5 dispatch to work.
    let twin = brep.half_edges.get(he)?.twin;
    let other_face = brep.loops.get(brep.half_edges.get(twin)?.loop_)?.face;
    if other_face == cap_face {
        return None; // self-loop / seam — not a cap-vs-other case
    }
    let other_surf = brep
        .faces
        .get(other_face)?
        .surface
        .as_ref()?
        .geometry
        .clone();
    if !matches!(
        other_surf,
        SurfaceGeometry::Plane { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. }
    ) {
        return None;
    }

    // Identify which Gregory UV-box boundary this edge sits on by
    // projecting the edge's 3-D midpoint into Gregory UV. Analytic
    // sphere caps do not need a boundary hint; the analytic
    // surface-pair dispatcher constructs the cap pcurve directly.
    let edge = brep.half_edges.get(he)?.edge;
    let coe = brep.edges.get(edge)?.curve.as_ref()?;
    let hint = if matches!(cap_surf, SurfaceGeometry::Gregory(_)) {
        let mid_t = (coe.t_start + coe.t_end) * 0.5;
        let p_mid = coe.geometry.point_at(mid_t);
        let (u_mid, v_mid) = project_point_to_uv(cap_surf, &p_mid)?;
        Some(EdgeHint::SurfaceBoundary {
            side: SurfaceSide::A,
            boundary: identify_uv_boundary(u_mid, v_mid)?,
        })
    } else {
        None
    };
    let _ = he; // reserved for future tracing

    // Build anchors from the BRep vertices at the half-edge endpoints.
    let v_origin = brep.vertices.get(brep.half_edges.get(he)?.origin)?.position;
    let v_dest = brep
        .vertices
        .get(brep.half_edges.get(twin)?.origin)?
        .position;
    let anchors = EdgeAnchors::Open {
        start: v_origin,
        end: v_dest,
    };
    let opts = EdgeTripleOptions::default();
    let triple = compute_edge_geometry(
        cap_surf,
        &other_surf,
        anchors,
        hint,
        &Tolerance::default(),
        &opts,
    )
    .ok()?;

    // Extract NurbsCurves from the Pcurve wrapper. Reproject path
    // expects NurbsCurve specifically; Pcurve currently stores
    // NurbsCurve internally.
    let cap_pc_nurbs = pcurve_to_nurbs(&triple.pcurve_a)?;
    let other_pc_nurbs = pcurve_to_nurbs(&triple.pcurve_b)?;

    Some((
        triple.curve_3d,
        triple.t_start,
        triple.t_end,
        cap_pc_nurbs,
        other_pc_nurbs,
    ))
}

/// Classify a UV point in Gregory `[0, 1]²` to one of the four
/// boundary sides if it's within `BOUNDARY_TOL_UV` of an axis.
/// Returns `None` for interior points.
fn identify_uv_boundary(u: f64, v: f64) -> Option<BoundarySide> {
    let d_umin = u.abs();
    let d_umax = (1.0 - u).abs();
    let d_vmin = v.abs();
    let d_vmax = (1.0 - v).abs();
    let min = d_umin.min(d_umax).min(d_vmin).min(d_vmax);
    if min > BOUNDARY_TOL_UV {
        return None;
    }
    if min == d_umin {
        Some(BoundarySide::UMin)
    } else if min == d_umax {
        Some(BoundarySide::UMax)
    } else if min == d_vmin {
        Some(BoundarySide::VMin)
    } else {
        Some(BoundarySide::VMax)
    }
}

/// Extract the underlying NurbsCurve from a `Pcurve` wrapper.
fn pcurve_to_nurbs(pc: &popkern_core::geometry::pcurve::Pcurve) -> Option<NurbsCurve> {
    Some(pc.to_nurbs())
}
