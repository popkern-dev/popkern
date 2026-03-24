//! Intersection graph data structures for the graph-based merge algorithm.
//!
//! The intersection graph represents the combined topology of two BRep operands
//! plus their intersection curves, using indexed arrays of nodes and half-edges.

use std::collections::HashMap;

use popkern_core::{
    geometry::{CurveOnEdge, Point3},
    topology::{EdgeHandle, FaceHandle, VertexHandle},
};

/// Identifies which boolean operand a piece of topology comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Operand {
    /// The first (left-hand) operand.
    A,
    /// The second (right-hand) operand.
    B,
}

/// A face qualified by which operand it belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct FaceRef {
    /// Which operand this face belongs to.
    pub operand: Operand,
    /// Handle to the face in its source BRep.
    pub face: FaceHandle,
}

/// Tracks the origin of a graph node.
#[derive(Debug, Clone)]
pub(crate) enum NodeSource {
    /// Copied from an original vertex in operand A.
    OriginalA(VertexHandle),
    /// Copied from an original vertex in operand B.
    OriginalB(VertexHandle),
    /// Created at a curve-curve or edge-face intersection point.
    Intersection,
}

/// Tracks the origin of a graph edge.
#[derive(Debug, Clone)]
pub(crate) enum EdgeSource {
    /// Copied from an original edge in operand A.
    OriginalA(EdgeHandle),
    /// Copied from an original edge in operand B.
    OriginalB(EdgeHandle),
    /// A shared edge lying on the intersection of two faces.
    Shared {
        /// Face from operand A.
        face_a: FaceHandle,
        /// Face from operand B.
        face_b: FaceHandle,
    },
    /// Created by splitting an original edge from operand A.
    SplitA {
        /// The original edge that was split.
        original: EdgeHandle,
    },
    /// Created by splitting an original edge from operand B.
    SplitB {
        /// The original edge that was split.
        original: EdgeHandle,
    },
}

/// A node (vertex) in the intersection graph.
#[derive(Debug)]
pub(crate) struct GraphNode {
    /// 3D position.
    pub position: Point3,
    /// UV coordinates of this node on each adjacent face.
    pub uv_coords: HashMap<FaceRef, (f64, f64)>,
    /// Where this node came from.
    pub source: NodeSource,
    /// Indices of outgoing half-edges in the graph's `half_edges` array.
    pub outgoing: Vec<usize>,
}

/// A directed half-edge in the intersection graph.
#[derive(Debug)]
pub(crate) struct GraphHalfEdge {
    /// Index of the origin node.
    pub origin: usize,
    /// Index of the twin (opposite-direction) half-edge.
    pub twin: usize,
    /// The face this half-edge borders, if assigned.
    pub face: Option<FaceRef>,
    /// Geometric curve along this half-edge, if any.
    pub curve: Option<CurveOnEdge>,
    /// Where this edge came from.
    pub source: EdgeSource,
}

/// The intersection graph: a combined representation of two operands'
/// topology plus intersection curves, ready for loop walking and face
/// reconstruction.
#[derive(Debug)]
pub(crate) struct IntersectionGraph {
    /// All nodes (vertices) in the graph.
    pub nodes: Vec<GraphNode>,
    /// All directed half-edges in the graph.
    pub half_edges: Vec<GraphHalfEdge>,
    /// Maps each face to the half-edge indices that border it.
    pub face_edges: HashMap<FaceRef, Vec<usize>>,
}

impl IntersectionGraph {
    /// Create an empty intersection graph.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            half_edges: Vec::new(),
            face_edges: HashMap::new(),
        }
    }

    /// Add a node at the given position with the given source.
    /// Returns the index of the new node.
    pub fn add_node(&mut self, position: Point3, source: NodeSource) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(GraphNode {
            position,
            uv_coords: HashMap::new(),
            source,
            outgoing: Vec::new(),
        });
        idx
    }

    /// Add a half-edge pair (forward + twin) between two nodes.
    ///
    /// Both half-edges share the same `source` and `curve`. The twin's curve
    /// is set to `None` (the forward half-edge owns the geometry).
    /// Returns `(forward_idx, twin_idx)`.
    pub fn add_half_edge_pair(
        &mut self,
        start: usize,
        end: usize,
        source: EdgeSource,
        curve: Option<CurveOnEdge>,
    ) -> (usize, usize) {
        let fwd_idx = self.half_edges.len();
        let twin_idx = fwd_idx + 1;

        self.half_edges.push(GraphHalfEdge {
            origin: start,
            twin: twin_idx,
            face: None,
            curve,
            source: source.clone(),
        });
        self.half_edges.push(GraphHalfEdge {
            origin: end,
            twin: fwd_idx,
            face: None,
            curve: None,
            source,
        });

        self.nodes[start].outgoing.push(fwd_idx);
        self.nodes[end].outgoing.push(twin_idx);

        (fwd_idx, twin_idx)
    }

    /// Find a node whose position is within `tolerance` of `position`.
    /// Returns `None` if no such node exists. O(N) scan.
    pub fn find_node(&self, position: &Point3, tolerance: f64) -> Option<usize> {
        self.nodes
            .iter()
            .enumerate()
            .find(|(_, n)| n.position.distance_to(position) < tolerance)
            .map(|(i, _)| i)
    }

    /// Assign a face to a half-edge and update the `face_edges` index.
    pub fn set_half_edge_face(&mut self, he_idx: usize, face: FaceRef) {
        self.half_edges[he_idx].face = Some(face);
        self.face_edges.entry(face).or_default().push(he_idx);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_graph() {
        let g = IntersectionGraph::new();
        assert!(g.nodes.is_empty());
        assert!(g.half_edges.is_empty());
        assert!(g.face_edges.is_empty());
    }

    #[test]
    fn add_node_returns_index() {
        let mut g = IntersectionGraph::new();
        let i0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let i1 = g.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(g.nodes.len(), 2);
    }

    #[test]
    fn add_half_edge_pair_wiring() {
        let mut g = IntersectionGraph::new();
        let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = g.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);

        let edge_h = popkern_core::topology::EdgeHandle::default();
        let (fwd, twin) = g.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

        // Twin pairing
        assert_eq!(g.half_edges[fwd].twin, twin);
        assert_eq!(g.half_edges[twin].twin, fwd);

        // Origins
        assert_eq!(g.half_edges[fwd].origin, n0);
        assert_eq!(g.half_edges[twin].origin, n1);

        // Outgoing lists
        assert!(g.nodes[n0].outgoing.contains(&fwd));
        assert!(g.nodes[n1].outgoing.contains(&twin));
    }

    #[test]
    fn find_node_within_tolerance() {
        let mut g = IntersectionGraph::new();
        g.add_node(Point3::new(1.0, 2.0, 3.0), NodeSource::Intersection);

        let found = g.find_node(&Point3::new(1.0, 2.0, 3.0 + 1e-8), 1e-6);
        assert_eq!(found, Some(0));

        let not_found = g.find_node(&Point3::new(10.0, 0.0, 0.0), 1e-6);
        assert!(not_found.is_none());
    }

    #[test]
    fn set_half_edge_face_updates_index() {
        let mut g = IntersectionGraph::new();
        let n0 = g.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = g.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let edge_h = popkern_core::topology::EdgeHandle::default();
        let (fwd, _twin) = g.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(edge_h), None);

        let face_ref = FaceRef {
            operand: Operand::A,
            face: popkern_core::topology::FaceHandle::default(),
        };
        g.set_half_edge_face(fwd, face_ref);

        assert_eq!(g.half_edges[fwd].face, Some(face_ref));
        assert_eq!(g.face_edges[&face_ref], vec![fwd]);
    }
}
