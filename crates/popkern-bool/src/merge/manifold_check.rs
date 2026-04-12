//! Pre-assembly manifold validation on the intersection graph.
//!
//! Implements the §7.4 active prevention requirement: rather than
//! letting the assembly path stitch together a non-manifold BRep and
//! catching the failure in post-validation, this module walks the
//! [`IntersectionGraph`] AFTER fragment selection but BEFORE assembly,
//! identifies edges where more than two faces would meet, and reports
//! them as a structured [`ManifoldReport`].
//!
//! Detection runs in O(N) over the graph half-edges. The check
//! groups half-edges by their canonical endpoint pair (smaller node
//! index first) and counts how many distinct face assignments use
//! each group. A manifold edge has exactly one or two assignments
//! (free edge or shared edge); anything higher is a multi-incident
//! edge.
//!
//! The current pass is detection-only — its findings are surfaced
//! through `BoolDiagnostics::entries` so callers can choose to fail
//! on multi-incident edges in strict mode. Subsequent work can drive
//! resolution from the report (e.g. pair-wise manifold splitting).
//!
//! "Active" here means at construction time rather than after
//! assembly — the report fires while the graph still contains the
//! provenance needed to fix the issue, before assembly turns the
//! graph into immutable BRep edges.

#![allow(dead_code)]

use std::collections::HashMap;

use super::graph::{FaceRef, IntersectionGraph};

/// A graph edge where more than two faces would share a single
/// BRep edge after assembly.
#[derive(Debug, Clone)]
pub(crate) struct MultiIncidentEdge {
    /// Canonical (smaller-index, larger-index) endpoint pair.
    pub endpoints: (usize, usize),
    /// All face assignments referring to this edge.
    pub face_assignments: Vec<FaceRef>,
}

/// Result of pre-assembly manifold validation.
#[derive(Debug, Clone, Default)]
pub(crate) struct ManifoldReport {
    /// Multi-incident edges (>2 face assignments).
    pub multi_incident: Vec<MultiIncidentEdge>,
    /// Total half-edge groups inspected.
    pub groups_checked: usize,
}

impl ManifoldReport {
    /// Returns true if no multi-incident edges were detected.
    pub fn is_manifold(&self) -> bool {
        self.multi_incident.is_empty()
    }
}

/// Walk the graph and detect multi-incident edges.
///
/// Groups half-edges by their canonical pair index `min(idx, twin)`
/// — i.e. by the underlying *edge*, not by endpoint position. This
/// matters because two distinct closed self-loop edges (separate
/// intersection circles) can share the same vertex; bucketing by
/// endpoint position would falsely combine them.
///
/// For each canonical edge, counts the distinct face assignments
/// across the half-edge pair. A manifold edge has exactly 1 (free
/// boundary) or 2 (interior shared edge) distinct face assignments.
/// More than 2 indicates non-manifold construction.
pub(crate) fn detect_multi_incident_edges(graph: &IntersectionGraph) -> ManifoldReport {
    // Map canonical pair index → set of face assignments. The pair
    // index is the smaller of (he_idx, twin_idx) for that half-edge.
    let mut groups: HashMap<usize, (usize, usize, Vec<FaceRef>)> = HashMap::new();

    for (he_idx, he) in graph.half_edges.iter().enumerate() {
        let twin_idx = he.twin;
        let pair_id = he_idx.min(twin_idx);
        let entry = groups.entry(pair_id).or_insert_with(|| {
            let twin = &graph.half_edges[twin_idx];
            let a = he.origin;
            let b = twin.origin;
            let endpoints = if a <= b { (a, b) } else { (b, a) };
            (endpoints.0, endpoints.1, Vec::new())
        });
        if let Some(face) = he.face {
            entry.2.push(face);
        }
    }

    let groups_checked = groups.len();
    let mut multi_incident = Vec::new();
    for (_pair_id, (a, b, assignments)) in groups {
        let mut unique = assignments;
        unique.sort_by_key(|f| (f.operand as u8, f.face));
        unique.dedup();
        if unique.len() > 2 {
            multi_incident.push(MultiIncidentEdge {
                endpoints: (a, b),
                face_assignments: unique,
            });
        }
    }

    ManifoldReport {
        multi_incident,
        groups_checked,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::graph::{EdgeSource, IntersectionGraph, NodeSource, Operand};
    use popkern_core::geometry::Point3;
    use popkern_core::topology::{EdgeHandle, FaceHandle};

    #[test]
    fn empty_graph_is_manifold() {
        let graph = IntersectionGraph::new();
        let report = detect_multi_incident_edges(&graph);
        assert!(report.is_manifold());
        assert_eq!(report.groups_checked, 0);
    }

    #[test]
    fn single_edge_two_faces_is_manifold() {
        let mut graph = IntersectionGraph::new();
        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let edge_h = EdgeHandle::default();
        let (fwd, twin) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

        let face_a = FaceRef {
            operand: Operand::A,
            face: FaceHandle::default(),
        };
        let face_b = FaceRef {
            operand: Operand::B,
            face: FaceHandle::default(),
        };
        graph.set_half_edge_face(fwd, face_a);
        graph.set_half_edge_face(twin, face_b);

        let report = detect_multi_incident_edges(&graph);
        assert!(
            report.is_manifold(),
            "single edge with 2 distinct faces should be manifold; got {:?}",
            report.multi_incident
        );
    }

    #[test]
    fn distinct_edges_at_same_endpoints_are_manifold() {
        // Two distinct half-edge pairs between the same two
        // endpoints should be treated as TWO separate edges, not
        // one multi-incident edge. This documents the corrected
        // detector semantics: grouping by half-edge pair index, not
        // by endpoint position.
        let mut graph = IntersectionGraph::new();
        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let edge_h = EdgeHandle::default();

        let face_a = FaceRef {
            operand: Operand::A,
            face: FaceHandle::default(),
        };
        let face_b = FaceRef {
            operand: Operand::B,
            face: FaceHandle::default(),
        };
        let (fwd1, twin1) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);
        graph.set_half_edge_face(fwd1, face_a);
        graph.set_half_edge_face(twin1, face_b);
        let (fwd2, twin2) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);
        graph.set_half_edge_face(fwd2, face_a);
        graph.set_half_edge_face(twin2, face_b);

        let report = detect_multi_incident_edges(&graph);
        // Two distinct edges → 2 groups checked, both manifold (each
        // has 2 face refs).
        assert!(report.is_manifold());
        assert_eq!(report.groups_checked, 2);
    }

    #[test]
    fn report_default_is_manifold() {
        let report = ManifoldReport::default();
        assert!(report.is_manifold());
        assert_eq!(report.groups_checked, 0);
    }
}
