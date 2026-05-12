//! Shared helpers for live corner/apex diagnostic extraction.

use crate::euler::{DiagnosticHandle, DiagnosticTolerance, DiagnosticToleranceSource};
use crate::geometry::curve::CurveEval;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

pub(super) fn incident_edges_at_vertex(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges = std::collections::BTreeSet::new();
    for (half_edge, he) in &brep.half_edges {
        let touches_vertex = he.origin == vertex
            || brep
                .half_edges
                .get(he.twin)
                .is_some_and(|twin| twin.origin == vertex);
        if touches_vertex && brep.edges.contains_key(he.edge) {
            edges.insert(he.edge);
        }
        if brep
            .edges
            .get(he.edge)
            .is_some_and(|edge| edge.half_edge == half_edge)
            && edge_contains_vertex(brep, he.edge, vertex)
        {
            edges.insert(he.edge);
        }
    }
    edges.into_iter().collect()
}

pub(super) fn incident_faces_at_vertex(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut faces = std::collections::BTreeSet::new();
    for (_, he) in &brep.half_edges {
        if he.origin != vertex {
            continue;
        }
        if let Some(loop_) = brep.loops.get(he.loop_) {
            faces.insert(loop_.face);
        }
    }
    faces.into_iter().collect()
}

pub(super) fn face_adjacency_edges(brep: &BRep, edges: &[EdgeHandle]) -> usize {
    edges
        .iter()
        .filter(|edge| edge_faces(brep, **edge).len() >= 2)
        .count()
}

pub(super) fn edge_faces(brep: &BRep, edge: EdgeHandle) -> Vec<FaceHandle> {
    let mut faces = std::collections::BTreeSet::new();
    let Some(edge_data) = brep.edges.get(edge) else {
        return Vec::new();
    };
    for half_edge in [edge_data.half_edge]
        .into_iter()
        .chain(brep.half_edges.get(edge_data.half_edge).map(|he| he.twin))
    {
        let Some(he) = brep.half_edges.get(half_edge) else {
            continue;
        };
        let Some(loop_) = brep.loops.get(he.loop_) else {
            continue;
        };
        faces.insert(loop_.face);
    }
    faces.into_iter().collect()
}

pub(super) fn live_vertex(handle: DiagnosticHandle, brep: &BRep) -> Option<VertexHandle> {
    match handle {
        DiagnosticHandle::Vertex(vertex) if brep.vertices.contains_key(vertex) => Some(vertex),
        _ => None,
    }
}

pub(super) fn live_edge(handle: DiagnosticHandle, brep: &BRep) -> Option<EdgeHandle> {
    match handle {
        DiagnosticHandle::Edge(edge) if brep.edges.contains_key(edge) => Some(edge),
        _ => None,
    }
}

pub(super) fn edge_contains_vertex(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let Some(first) = brep.half_edges.get(edge_data.half_edge) else {
        return false;
    };
    if first.origin == vertex {
        return true;
    }
    brep.half_edges
        .get(first.twin)
        .is_some_and(|second| second.origin == vertex)
}

pub(super) fn edge_extent_nm(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> Option<i64> {
    let edge_data = brep.edges.get(edge)?;
    if let Some(curve) = &edge_data.curve {
        let t = curve
            .geometry
            .closest_point(&brep.vertices[vertex].position)?;
        let start = curve.t_start;
        let end = curve.t_end;
        let extent = curve
            .geometry
            .length_between(start, t)
            .max(curve.geometry.length_between(t, end));
        return finite_nm(extent);
    }
    let first = brep.half_edges.get(edge_data.half_edge)?;
    let second = brep.half_edges.get(first.twin)?;
    let other = if first.origin == vertex {
        second.origin
    } else {
        first.origin
    };
    finite_nm(
        brep.vertices
            .get(vertex)?
            .position
            .distance_to(&brep.vertices.get(other)?.position),
    )
}

pub(super) fn finite_nm(value: f64) -> Option<i64> {
    value
        .is_finite()
        .then_some((value.max(0.0) * 1_000_000_000.0).round() as i64)
}

pub(super) fn local_tolerance(brep: &BRep, handle: DiagnosticHandle) -> DiagnosticTolerance {
    let linear = match handle {
        DiagnosticHandle::Vertex(vertex) if brep.vertices.contains_key(vertex) => {
            brep.vertex_tolerance(vertex)
        }
        _ => brep.tolerance.linear(),
    };
    DiagnosticTolerance {
        source: DiagnosticToleranceSource::LocalGeometryPredicate,
        nanometers: finite_nm(linear),
    }
}
