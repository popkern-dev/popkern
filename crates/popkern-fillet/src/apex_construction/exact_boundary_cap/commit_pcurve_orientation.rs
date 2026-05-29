//! Live half-edge pcurve orientation helpers for exact-boundary commit.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use super::commit_helpers::reversed_pcurve;

pub(super) fn orient_pcurve_to_half_edge_vertices(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    pcurve: Pcurve,
    surface: &SurfaceGeometry,
) -> Pcurve {
    let Some(direct_gap) = pcurve_endpoint_gap(brep, half_edge, &pcurve, surface) else {
        return pcurve;
    };
    let reversed = reversed_pcurve(&pcurve);
    let Some(reversed_gap) = pcurve_endpoint_gap(brep, half_edge, &reversed, surface) else {
        return pcurve;
    };
    if reversed_gap < direct_gap {
        reversed
    } else {
        pcurve
    }
}

fn pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
) -> Option<f64> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let destination = brep
        .vertices
        .get(half_edge_destination_vertex(brep, half_edge)?)?
        .position;
    let (t_start, t_end) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t_start);
    let end_uv = pcurve.point_at(t_end);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    Some(
        start
            .distance_to(&origin)
            .max(end.distance_to(&destination)),
    )
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin| twin.origin)
}
