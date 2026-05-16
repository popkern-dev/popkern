//! Shared loop and loop-adjacency facts for BRep repair knowledge.

use crate::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Shared loop facts for repair diagnostics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BrepLoopFacts {
    /// Face that owns the loop.
    pub face: FaceHandle,
    /// Loop handle.
    pub loop_handle: LoopHandle,
    /// Loop role on the owning face.
    pub kind: BrepLoopKind,
    /// Ordered boundary half-edges observed while walking the loop.
    pub boundary_half_edges: Vec<HalfEdgeHandle>,
    /// Ordered boundary vertices observed while walking the loop.
    pub boundary_vertices: Vec<VertexHandle>,
    /// Whether a vertex repeated while walking the loop.
    pub repeated_vertex: bool,
    /// Whether a half-edge repeated while walking the loop.
    pub repeated_half_edge: bool,
}

impl BrepLoopFacts {
    /// Owning face.
    pub fn face(&self) -> FaceHandle {
        self.face
    }

    /// Loop handle.
    pub fn loop_handle(&self) -> LoopHandle {
        self.loop_handle
    }

    /// Loop kind.
    pub fn kind(&self) -> BrepLoopKind {
        self.kind
    }

    /// Boundary half-edges.
    pub fn boundary_half_edges(&self) -> &[HalfEdgeHandle] {
        &self.boundary_half_edges
    }

    /// Boundary vertices.
    pub fn boundary_vertices(&self) -> &[VertexHandle] {
        &self.boundary_vertices
    }

    /// Legacy-compatible compact diagnostic label.
    pub fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}/{:?} vertices={:?} repeated_vertex={} repeated_half_edge={}",
            self.face,
            self.loop_handle,
            self.kind,
            self.boundary_vertices,
            self.repeated_vertex,
            self.repeated_half_edge
        )
    }
}

/// Shared loop adjacency fact.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BrepLoopAdjacencyFacts {
    /// First face.
    pub face: FaceHandle,
    /// First loop.
    pub loop_handle: LoopHandle,
    /// Half-edge on the first loop.
    pub half_edge: HalfEdgeHandle,
    /// Adjacent face across the half-edge twin.
    pub adjacent_face: FaceHandle,
    /// Adjacent loop across the half-edge twin.
    pub adjacent_loop: LoopHandle,
    /// Adjacent half-edge.
    pub adjacent_half_edge: HalfEdgeHandle,
}

impl BrepLoopAdjacencyFacts {
    /// First face.
    pub fn face(&self) -> FaceHandle {
        self.face
    }

    /// First loop.
    pub fn loop_handle(&self) -> LoopHandle {
        self.loop_handle
    }

    /// Half-edge on the first loop.
    pub fn half_edge(&self) -> HalfEdgeHandle {
        self.half_edge
    }

    /// Adjacent face.
    pub fn adjacent_face(&self) -> FaceHandle {
        self.adjacent_face
    }

    /// Adjacent loop.
    pub fn adjacent_loop(&self) -> LoopHandle {
        self.adjacent_loop
    }

    /// Adjacent half-edge.
    pub fn adjacent_half_edge(&self) -> HalfEdgeHandle {
        self.adjacent_half_edge
    }

    /// Legacy-compatible compact diagnostic label.
    pub fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}@{:?}<->{:?}/{:?}@{:?}",
            self.face,
            self.loop_handle,
            self.half_edge,
            self.adjacent_face,
            self.adjacent_loop,
            self.adjacent_half_edge
        )
    }
}

/// Loop role on the owning face.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepLoopKind {
    /// Outer loop.
    Outer,
    /// Inner loop.
    Inner,
    /// Loop role could not be determined.
    Unknown,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn loop_kind_debug_labels_match_legacy_apex_labels() {
        assert_eq!(format!("{:?}", BrepLoopKind::Outer), "Outer");
        assert_eq!(format!("{:?}", BrepLoopKind::Inner), "Inner");
        assert_eq!(format!("{:?}", BrepLoopKind::Unknown), "Unknown");
    }

    #[test]
    fn loop_fact_label_contains_legacy_fields() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has face");
        let loop_handle = brep.faces[face].outer_loop;
        let half_edge = brep.loops[loop_handle].half_edge;
        let vertex = brep.half_edges[half_edge].origin;
        let facts = BrepLoopFacts {
            face,
            loop_handle,
            kind: BrepLoopKind::Outer,
            boundary_half_edges: vec![half_edge],
            boundary_vertices: vec![vertex],
            repeated_vertex: false,
            repeated_half_edge: false,
        };

        let label = facts.diagnostic_label();
        assert!(label.contains("Outer"), "{label}");
        assert!(label.contains("repeated_vertex=false"), "{label}");
        assert!(label.contains("repeated_half_edge=false"), "{label}");
    }
}
