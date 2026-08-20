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
/// Groups distinct half-edge PAIRS by their canonical endpoint pair
/// (start_vertex, end_vertex). A "pair" here is one half-edge plus
/// its twin — i.e. the underlying edge. We then count how many
/// such pairs share the same endpoint pair AND check whether their
/// face assignments produce a non-manifold result.
///
/// A manifold edge has exactly 1 or 2 unique face refs across the
/// pair (1 = free boundary on a sheet body; 2 = interior shared
/// edge between two faces). If multiple half-edge pairs share the
/// same endpoint pair, the union of all their face refs is what
/// matters: 3 or more distinct face refs means non-manifold after
/// dedup/merge during assembly.
///
/// Closed self-loops (start == end) are excluded from grouping
/// because each is a separate edge geometrically — multiple
/// intersection circles at the same vertex are NOT a multi-incident
/// edge.
pub(crate) fn detect_multi_incident_edges(graph: &IntersectionGraph) -> ManifoldReport {
    use std::collections::HashSet;

    // Group by (start_vertex, end_vertex) canonical pair, skipping
    // closed self-loops. The value is the set of distinct face
    // refs collected from all half-edge pairs in the group.
    let mut groups: HashMap<(usize, usize), HashSet<FaceRef>> = HashMap::new();
    let mut seen_pairs: HashSet<usize> = HashSet::new();

    for (he_idx, he) in graph.half_edges.iter().enumerate() {
        let twin_idx = he.twin;
        let pair_id = he_idx.min(twin_idx);
        if !seen_pairs.insert(pair_id) {
            continue;
        }
        let twin = &graph.half_edges[twin_idx];
        let a = he.origin;
        let b = twin.origin;
        if a == b {
            // Closed self-loop: not a candidate for multi-incidence.
            continue;
        }
        let key = if a <= b { (a, b) } else { (b, a) };
        let entry = groups.entry(key).or_default();
        if let Some(face) = he.face {
            entry.insert(face);
        }
        if let Some(face) = twin.face {
            entry.insert(face);
        }
    }

    let groups_checked = groups.len();
    let mut multi_incident = Vec::new();
    for (endpoints, face_set) in groups {
        if face_set.len() > 2 {
            let mut face_assignments: Vec<FaceRef> = face_set.into_iter().collect();
            face_assignments.sort_by_key(|f| (f.operand as u8, f.face));
            multi_incident.push(MultiIncidentEdge {
                endpoints,
                face_assignments,
            });
        }
    }
    // `groups` is a HashMap, so the loop above pushed entries in hash
    // order; the report Vec is destined for BoolDiagnostics::entries,
    // where order must be canonical (PK-300 hash-iteration audit).
    // Endpoints are already canonical min/max node pairs.
    multi_incident.sort_unstable_by_key(|e| e.endpoints);

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
    fn distinct_edges_at_same_endpoints_share_one_group() {
        // Two distinct half-edge pairs between the same two
        // endpoints collapse to a single endpoint-pair group. The
        // union of their face refs is {A, B} = 2 distinct, so the
        // group is manifold.
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
        // One endpoint-pair group with 2 unique face refs → manifold.
        assert!(report.is_manifold());
        assert_eq!(report.groups_checked, 1);
    }

    #[test]
    fn report_default_is_manifold() {
        let report = ManifoldReport::default();
        assert!(report.is_manifold());
        assert_eq!(report.groups_checked, 0);
    }
}
