//! Live topology-neighborhood extraction for corner/apex diagnostics.

use super::types::*;
use crate::euler::DiagnosticHandle;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

pub(super) fn extract_topology_neighborhood(
    brep: &BRep,
    vertex: VertexHandle,
    incident_edges: &[EdgeHandle],
    selected_edges: &[LocalBlendEdgeSelection],
) -> CornerTopologyNeighborhood {
    let selected = selected_edges
        .iter()
        .filter_map(|selection| match selection.edge {
            DiagnosticHandle::Edge(edge) => Some(edge),
            _ => None,
        })
        .collect::<std::collections::BTreeSet<_>>();
    let mut neighborhood = CornerTopologyNeighborhood::default();
    for edge in incident_edges {
        neighborhood.incident_edges.push(CornerIncidentEdgeFact {
            edge: DiagnosticHandle::Edge(*edge),
            endpoint_role: endpoint_role(brep, *edge, vertex),
            selected_for_blend: selected.contains(edge),
        });
        let pair = adjacent_face_pair(brep, *edge);
        neighborhood.adjacent_face_pairs.push(pair);
        if let (Some(DiagnosticHandle::Face(first)), Some(DiagnosticHandle::Face(second))) =
            (pair.first_face, pair.second_face)
            && face_contains_vertex(brep, first, vertex)
            && face_contains_vertex(brep, second, vertex)
        {
            neighborhood
                .surrounding_face_adjacencies
                .push(CornerFaceAdjacencyFact {
                    first_face: DiagnosticHandle::Face(first),
                    second_face: DiagnosticHandle::Face(second),
                    shared_edge: DiagnosticHandle::Edge(*edge),
                });
        }
    }
    neighborhood.incident_edges.sort();
    neighborhood.adjacent_face_pairs.sort();
    neighborhood.surrounding_face_adjacencies.sort();
    neighborhood
}

pub(super) fn stale_topology_neighborhood(
    source: CornerTopologyFactSource,
) -> CornerTopologyNeighborhood {
    CornerTopologyNeighborhood {
        stale_sources: vec![source],
        ..CornerTopologyNeighborhood::default()
    }
}

fn endpoint_role(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> EdgeEndpointRole {
    let Some(edge_record) = brep.edges.get(edge) else {
        return EdgeEndpointRole::Unknown;
    };
    let Some(half_edge) = brep.half_edges.get(edge_record.half_edge) else {
        return EdgeEndpointRole::Unknown;
    };
    let Some(twin) = brep.half_edges.get(half_edge.twin) else {
        return EdgeEndpointRole::Unknown;
    };
    match (half_edge.origin == vertex, twin.origin == vertex) {
        (true, true) => EdgeEndpointRole::Degenerate,
        (true, false) => EdgeEndpointRole::Start,
        (false, true) => EdgeEndpointRole::End,
        (false, false) => EdgeEndpointRole::Unknown,
    }
}

fn adjacent_face_pair(brep: &BRep, edge: EdgeHandle) -> CornerAdjacentFacePair {
    let Some(edge_record) = brep.edges.get(edge) else {
        return missing_pair(edge);
    };
    let Some(first_half_edge) = brep.half_edges.get(edge_record.half_edge) else {
        return missing_pair(edge);
    };
    let first_face = face_for_half_edge(brep, edge_record.half_edge);
    let second_face = face_for_half_edge(brep, first_half_edge.twin);
    CornerAdjacentFacePair {
        edge: DiagnosticHandle::Edge(edge),
        first_face: first_face.map(DiagnosticHandle::Face),
        second_face: second_face.map(DiagnosticHandle::Face),
    }
}

fn missing_pair(edge: EdgeHandle) -> CornerAdjacentFacePair {
    CornerAdjacentFacePair {
        edge: DiagnosticHandle::Edge(edge),
        first_face: None,
        second_face: None,
    }
}

fn face_for_half_edge(
    brep: &BRep,
    half_edge: crate::topology::HalfEdgeHandle,
) -> Option<FaceHandle> {
    let loop_ = brep.half_edges.get(half_edge)?.loop_;
    Some(brep.loops.get(loop_)?.face)
}

fn face_contains_vertex(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    loop_contains_vertex(brep, brep.faces[face].outer_loop, vertex)
        || brep.faces[face]
            .inner_loops
            .iter()
            .any(|loop_| loop_contains_vertex(brep, *loop_, vertex))
}

fn loop_contains_vertex(
    brep: &BRep,
    loop_: crate::topology::LoopHandle,
    vertex: VertexHandle,
) -> bool {
    let start = brep.loops[loop_].half_edge;
    let mut cursor = start;
    loop {
        let half_edge = &brep.half_edges[cursor];
        if half_edge.origin == vertex {
            return true;
        }
        cursor = half_edge.next;
        if cursor == start {
            return false;
        }
    }
}
