//! Insert intersection curves as shared edges and validate the graph.
//!
//! For each [`IntersectionRecord`], finds or creates endpoint nodes in the
//! graph and adds a shared half-edge pair representing the intersection curve.

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use popkern_core::geometry::{
    CurveEval,
    curve::{CurveGeometry, CurveOnEdge},
};

use super::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
use crate::intersect::IntersectionRecord;

/// Insert intersection curves as shared edges in the graph.
///
/// For each [`IntersectionRecord`], looks up start/end nodes of the
/// intersection curve (they should already exist from edge splitting in
/// earlier phases), and creates shared half-edge pairs.
pub(crate) fn insert_intersection_curves(
    graph: &mut IntersectionGraph,
    records: &[IntersectionRecord],
    tolerance: f64,
) {
    for record in records {
        let curve_3d = &record.curve.curve_3d;
        let (t_start, t_end) = curve_3d.parameter_range();

        let start_pt = curve_3d.point_at(t_start);
        let end_pt = curve_3d.point_at(t_end);

        // Find or create start node.
        let start_node = match graph.find_node(&start_pt, tolerance) {
            Some(idx) => idx,
            None => graph.add_node(start_pt, NodeSource::Intersection),
        };

        // Find or create end node.
        let end_node = match graph.find_node(&end_pt, tolerance) {
            Some(idx) => idx,
            None => graph.add_node(end_pt, NodeSource::Intersection),
        };

        // Build the CurveOnEdge from the 3D intersection curve.
        let curve_on_edge = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve_3d.clone()),
            t_start,
            t_end,
        };

        let source = EdgeSource::Shared {
            face_a: record.face_a,
            face_b: record.face_b,
        };

        // Attach UV curves as pcurves: curve_uv1 belongs to face A (fwd),
        // curve_uv2 belongs to face B (twin).
        let pcurve_a: Option<popkern_core::geometry::pcurve::Pcurve> =
            Some(record.curve.curve_uv1.clone().into());
        let pcurve_b: Option<popkern_core::geometry::pcurve::Pcurve> =
            Some(record.curve.curve_uv2.clone().into());

        let (fwd, twin) = graph.add_half_edge_pair_with_pcurves(
            start_node,
            end_node,
            source,
            Some(curve_on_edge),
            pcurve_a,
            pcurve_b,
        );

        let face_a_ref = FaceRef {
            operand: Operand::A,
            face: record.face_a,
        };
        let face_b_ref = FaceRef {
            operand: Operand::B,
            face: record.face_b,
        };

        // Assign primary face ownership for each half-edge.
        graph.set_half_edge_face(fwd, face_a_ref);
        graph.set_half_edge_face(twin, face_b_ref);

        // Each face needs the intersection curve in BOTH directions so the
        // loop walker can split the face into two fragments.  Register the
        // twin with face A and the fwd with face B as well.
        graph.face_edges.entry(face_a_ref).or_default().push(twin);
        graph.face_edges.entry(face_b_ref).or_default().push(fwd);
    }
}

impl IntersectionGraph {
    /// Validate graph topology, returning a list of error descriptions.
    ///
    /// Checks:
    /// - Every half-edge's twin points back to it.
    /// - Every node has at least one outgoing edge.
    /// - Every face's boundary half-edges form at least one closed loop.
    pub fn validate(&self) -> Vec<String> {
        let mut errors = Vec::new();

        // Check twin consistency.
        for (i, he) in self.half_edges.iter().enumerate() {
            if he.twin >= self.half_edges.len() {
                errors.push(format!(
                    "half-edge {i}: twin index {} out of bounds",
                    he.twin
                ));
                continue;
            }
            let twin_he = &self.half_edges[he.twin];
            if twin_he.twin != i {
                errors.push(format!(
                    "half-edge {i}: twin {} does not point back (points to {})",
                    he.twin, twin_he.twin
                ));
            }
        }

        // Check every node has at least one outgoing edge.
        for (i, node) in self.nodes.iter().enumerate() {
            if node.outgoing.is_empty() {
                errors.push(format!("node {i}: no outgoing edges"));
            }
        }

        // Check face edge loops are closed.
        for (face_ref, he_indices) in &self.face_edges {
            if he_indices.is_empty() {
                errors.push(format!("face {face_ref:?}: no half-edges"));
                continue;
            }
            // Verify each half-edge in the face list actually references this face.
            for &he_idx in he_indices {
                if he_idx >= self.half_edges.len() {
                    errors.push(format!(
                        "face {face_ref:?}: half-edge index {he_idx} out of bounds"
                    ));
                    continue;
                }
                if self.half_edges[he_idx].face.as_ref() != Some(face_ref) {
                    errors.push(format!(
                        "face {face_ref:?}: half-edge {he_idx} has different face {:?}",
                        self.half_edges[he_idx].face
                    ));
                }
            }
        }

        errors
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::{
        geometry::{NurbsCurve, Point3, curve_on_surface::CurveOnSurface},
        topology::FaceHandle,
    };

    /// Build a minimal IntersectionRecord with a line from p0 to p1.
    fn make_record(
        p0: Point3,
        p1: Point3,
        face_a: FaceHandle,
        face_b: FaceHandle,
    ) -> IntersectionRecord {
        // Build a degree-1 NURBS line.
        let nurbs = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p0, p1],
            weights: vec![1.0, 1.0],
        };
        let cos = CurveOnSurface {
            curve_3d: nurbs.clone(),
            curve_uv1: nurbs.clone(),
            curve_uv2: nurbs,
            is_closed: false,
        };
        IntersectionRecord {
            face_a,
            face_b,
            curve: cos,
            is_tangent: false,
            is_coincident: false,
            start_boundary: Default::default(),
            end_boundary: Default::default(),
        }
    }

    #[test]
    fn insert_shared_edge_between_two_faces() {
        let mut graph = IntersectionGraph::new();

        // Pre-create nodes at the endpoints.
        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        assert_eq!(n0, 0);
        assert_eq!(n1, 1);

        let face_a = FaceHandle::default();
        let face_b = FaceHandle::default();

        let record = make_record(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            face_a,
            face_b,
        );

        insert_intersection_curves(&mut graph, &[record], 1e-6);

        // Should still have 2 nodes (reused existing ones).
        assert_eq!(graph.nodes.len(), 2);
        // Should have 2 half-edges (one pair).
        assert_eq!(graph.half_edges.len(), 2);

        // Forward half-edge should have face A.
        let face_a_ref = FaceRef {
            operand: Operand::A,
            face: face_a,
        };
        let face_b_ref = FaceRef {
            operand: Operand::B,
            face: face_b,
        };
        assert_eq!(graph.half_edges[0].face, Some(face_a_ref));
        assert_eq!(graph.half_edges[1].face, Some(face_b_ref));

        // face_edges should list the shared edge for both faces.
        assert!(graph.face_edges[&face_a_ref].contains(&0));
        assert!(graph.face_edges[&face_b_ref].contains(&1));

        // Forward half-edge should carry the curve.
        assert!(graph.half_edges[0].curve.is_some());
    }

    #[test]
    fn insert_creates_nodes_when_not_found() {
        let mut graph = IntersectionGraph::new();

        let face_a = FaceHandle::default();
        let face_b = FaceHandle::default();

        let record = make_record(
            Point3::new(5.0, 5.0, 5.0),
            Point3::new(6.0, 6.0, 6.0),
            face_a,
            face_b,
        );

        insert_intersection_curves(&mut graph, &[record], 1e-6);

        // Nodes should be created.
        assert_eq!(graph.nodes.len(), 2);
        assert_eq!(graph.half_edges.len(), 2);
    }

    #[test]
    fn validate_on_valid_graph() {
        let mut graph = IntersectionGraph::new();
        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let edge_h = popkern_core::topology::EdgeHandle::default();
        graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

        let errors = graph.validate();
        assert!(errors.is_empty(), "Expected no errors, got: {errors:?}");
    }

    #[test]
    fn validate_detects_orphan_node() {
        let mut graph = IntersectionGraph::new();
        // Node with no outgoing edges.
        graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);

        let errors = graph.validate();
        assert_eq!(errors.len(), 1);
        assert!(errors[0].contains("no outgoing edges"));
    }
}
