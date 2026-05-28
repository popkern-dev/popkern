//! Support-loop surface checks for exact-boundary commit.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

pub(super) fn assign_loop_valid_support_surface(
    brep: &mut BRep,
    support_face: FaceHandle,
    support_surface: SurfaceGeometry,
) -> Option<SurfaceGeometry> {
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let residual = support_loop_surface_residual(brep, support_face, &support_surface)?;
    if residual > tolerance {
        return None;
    }
    let face = brep.faces.get_mut(support_face)?;
    face.surface = Some(SurfaceOnFace {
        geometry: support_surface.clone(),
        orientation: true,
    });
    Some(support_surface)
}

pub(super) fn support_loop_surface_residual(
    brep: &BRep,
    support_face: FaceHandle,
    support_surface: &SurfaceGeometry,
) -> Option<f64> {
    let loop_handle = brep.faces.get(support_face)?.outer_loop;
    collect_loop_half_edges(brep, loop_handle)?
        .into_iter()
        .filter_map(|half_edge| {
            let edge = brep.edges.get(brep.half_edges.get(half_edge)?.edge)?;
            let curve = edge.curve.as_ref()?;
            Some(
                edge_surface_deviation(
                    &curve.geometry,
                    (curve.t_start, curve.t_end),
                    support_surface,
                    9,
                )
                .max_deviation,
            )
        })
        .reduce(f64::max)
}

fn collect_loop_half_edges(brep: &BRep, loop_handle: LoopHandle) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}
