//! Live outer-seam refresh for materialized exact-boundary sectors.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::super::sector::ExactBoundaryCapSectorInput;
use super::{ExactBoundaryTopologyDeclineReason, ExactBoundaryTopologyDiagnostic};

pub(super) fn refresh_sector_outer_seams(
    brep: &BRep,
    sector: &mut ExactBoundaryCapSectorInput,
) -> Result<(), ExactBoundaryTopologyDiagnostic> {
    let fillet = live_outer_seam(
        brep,
        sector,
        sector.fillet_end_edge,
        sector.fillet_support_face,
    )?;
    let face = live_outer_seam(brep, sector, sector.face_edge, sector.bridge_support_face)?;
    sector.fillet_support_face = fillet.1;
    sector.bridge_support_face = face.1;
    sector.cap_half_edges = [fillet.0, face.0];
    sector.support_half_edges = [fillet.2, face.2];
    sector.outer_seams[0].cap_half_edge = fillet.0;
    sector.outer_seams[0].support_face = fillet.1;
    sector.outer_seams[0].support_half_edge = fillet.2;
    sector.outer_seams[0].cap_pcurve_present = brep
        .half_edges
        .get(fillet.0)
        .is_some_and(|half_edge| half_edge.pcurve.is_some());
    sector.outer_seams[0].support_pcurve_present = brep
        .half_edges
        .get(fillet.2)
        .is_some_and(|half_edge| half_edge.pcurve.is_some());
    sector.outer_seams[1].cap_half_edge = face.0;
    sector.outer_seams[1].support_face = face.1;
    sector.outer_seams[1].support_half_edge = face.2;
    sector.outer_seams[1].cap_pcurve_present = brep
        .half_edges
        .get(face.0)
        .is_some_and(|half_edge| half_edge.pcurve.is_some());
    sector.outer_seams[1].support_pcurve_present = brep
        .half_edges
        .get(face.2)
        .is_some_and(|half_edge| half_edge.pcurve.is_some());
    Ok(())
}

fn live_outer_seam(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
    edge: EdgeHandle,
    expected_support_face: FaceHandle,
) -> Result<(HalfEdgeHandle, FaceHandle, HalfEdgeHandle), ExactBoundaryTopologyDiagnostic> {
    let edge_data = brep.edges.get(edge).ok_or_else(|| {
        ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
            Some(sector.sector_name.clone()),
            None,
        )
    })?;
    let primary = edge_data.half_edge;
    let twin = brep.half_edges.get(primary).map(|half_edge| half_edge.twin);
    for half_edge in [Some(primary), twin] {
        let Some(half_edge) = half_edge else {
            continue;
        };
        if face_for_half_edge(brep, half_edge) != Some(sector.cap_face) {
            continue;
        }
        let support_half_edge = brep
            .half_edges
            .get(half_edge)
            .map(|data| data.twin)
            .ok_or_else(|| {
                ExactBoundaryTopologyDiagnostic::new(
                    ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
                    Some(sector.sector_name.clone()),
                    None,
                )
            })?;
        let support_face = face_for_half_edge(brep, support_half_edge).ok_or_else(|| {
            ExactBoundaryTopologyDiagnostic::new(
                ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
                Some(sector.sector_name.clone()),
                None,
            )
        })?;
        if support_face != sector.cap_face {
            if support_face != expected_support_face
                && let Some(semantic_support_half_edge) =
                    same_oriented_side_on_face(brep, half_edge, expected_support_face)
            {
                return Ok((half_edge, expected_support_face, semantic_support_half_edge));
            }
            let expected_support_is_cylinder = brep
                .faces
                .get(expected_support_face)
                .and_then(|face| face.surface.as_ref())
                .is_some_and(|surface| {
                    matches!(surface.geometry, SurfaceGeometry::Cylinder { .. })
                });
            if expected_support_is_cylinder
                && brep
                    .faces
                    .get(support_face)
                    .and_then(|face| face.surface.as_ref())
                    .is_none()
                && let Some((semantic_support_half_edge, semantic_support_face)) =
                    same_oriented_cylinder_support_side(brep, sector.cap_face, half_edge)
            {
                return Ok((half_edge, semantic_support_face, semantic_support_half_edge));
            }
            return Ok((half_edge, support_face, support_half_edge));
        }
    }
    Err(ExactBoundaryTopologyDiagnostic::new(
        ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
        Some(sector.sector_name.clone()),
        None,
    ))
}

fn same_oriented_side_on_face(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
    expected_support_face: FaceHandle,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(expected_support_face)?;
    let sides = collect_loop_half_edges(brep, face_data.outer_loop)?;
    sides
        .into_iter()
        .find(|side| same_oriented_side(brep, cap_half_edge, *side))
}

fn same_oriented_cylinder_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if same_oriented_side(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn same_oriented_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != support_start {
        return false;
    }
    let Some(cap_end) = super::half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = super::half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    cap_end == support_end
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let loop_data = brep.loops.get(loop_handle)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(cycle);
        }
    }
    None
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    brep.half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face)
}
