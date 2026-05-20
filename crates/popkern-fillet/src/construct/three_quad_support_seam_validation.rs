//! Validation helpers for relocation-backed support seams.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

pub(super) fn validate_moved_vertex_curve_coherence(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    tolerance: f64,
) -> bool {
    let Some(twin) = brep.half_edges.get(half_edge).map(|he| he.twin) else {
        return false;
    };
    let Some(start_vertex) = brep.half_edges.get(half_edge).map(|he| he.origin) else {
        return false;
    };
    let Some(end_vertex) = brep.half_edges.get(twin).map(|he| he.origin) else {
        return false;
    };
    [start_vertex, end_vertex]
        .into_iter()
        .all(|vertex| incident_curves_touch_vertex(brep, vertex, tolerance))
}

fn incident_curves_touch_vertex(brep: &BRep, vertex: VertexHandle, tolerance: f64) -> bool {
    let vertex_position = brep.vertices[vertex].position;
    brep.edges.iter().all(|(_, edge)| {
        let primary = edge.half_edge;
        let Some(twin) = brep.half_edges.get(primary).map(|he| he.twin) else {
            return true;
        };
        let Some(a) = brep.half_edges.get(primary).map(|he| he.origin) else {
            return true;
        };
        let Some(b) = brep.half_edges.get(twin).map(|he| he.origin) else {
            return true;
        };
        if a != vertex && b != vertex {
            return true;
        }
        let Some(curve) = edge.curve.as_ref() else {
            return true;
        };
        let p0 = curve.geometry.point_at(curve.t_start);
        let p1 = curve.geometry.point_at(curve.t_end);
        let gap = vertex_position
            .distance_to(&p0)
            .min(vertex_position.distance_to(&p1));
        if gap > tolerance && crate::direct_apex_trace::cap_g0_enabled() {
            eprintln!(
                "[direct-apex-trace cap-g0] moved vertex coherence failed vertex={vertex:?} edge={:?} gap={gap:.6e}",
                edge.half_edge
            );
        }
        gap <= tolerance
    })
}

pub(super) fn validate_support_seam(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    twin: HalfEdgeHandle,
    cap_surface: &SurfaceGeometry,
    tolerance: f64,
    sample_count: usize,
) -> Result<f64, f64> {
    let edge = brep.half_edges[half_edge].edge;
    let curve = brep.edges[edge].curve.as_ref().ok_or(f64::INFINITY)?;
    let cap_pcurve = brep.half_edges[half_edge]
        .pcurve
        .as_ref()
        .ok_or(f64::INFINITY)?;
    let support_face = brep.loops[brep.half_edges[twin].loop_].face;
    let support_surface = &brep.faces[support_face]
        .surface
        .as_ref()
        .ok_or(f64::INFINITY)?
        .geometry;
    let support_pcurve = brep.half_edges[twin].pcurve.as_ref().ok_or(f64::INFINITY)?;
    let mut worst_curve_cap = 0.0_f64;
    let mut worst_curve_support = 0.0_f64;
    let mut worst_cap_support = 0.0_f64;
    for i in 0..sample_count {
        let s = i as f64 / (sample_count - 1) as f64;
        let t = curve.t_start + s * (curve.t_end - curve.t_start);
        let p_curve = curve.geometry.point_at(t);
        let cap_uv = cap_pcurve.point_at(t);
        let support_uv = support_pcurve.point_at(t);
        let p_cap = cap_surface.point_at(cap_uv.x, cap_uv.y);
        let p_support = support_surface.point_at(support_uv.x, support_uv.y);
        worst_curve_cap = worst_curve_cap.max(p_curve.distance_to(&p_cap));
        worst_curve_support = worst_curve_support.max(p_curve.distance_to(&p_support));
        worst_cap_support = worst_cap_support.max(p_cap.distance_to(&p_support));
    }
    let worst = worst_curve_cap
        .max(worst_curve_support)
        .max(worst_cap_support);
    if worst > tolerance && crate::direct_apex_trace::cap_g0_enabled() {
        eprintln!(
            "[direct-apex-trace cap-g0] support seam validation he={half_edge:?}: \
             curve-cap={worst_curve_cap:.6e} curve-support={worst_curve_support:.6e} \
             cap-support={worst_cap_support:.6e}"
        );
    }
    if worst <= tolerance {
        Ok(worst)
    } else {
        Err(worst)
    }
}
