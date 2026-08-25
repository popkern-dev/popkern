//! Handle mapping between outer and inner skin BReps.

use std::collections::HashMap;

use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Provenance of one join face constructed by the gap-join stage
/// (design §5 Stage 3a, PK-325). Handles reference the INPUT body for
/// sources and the produced skin for `face`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JoinRecord {
    /// What the join fills and where it came from.
    pub kind: JoinKind,
    /// The join face in the produced skin.
    pub face: FaceHandle,
}

/// The kind of gap a join face fills.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JoinKind {
    /// Rolling-ball strip along a gap edge (cylinder/torus/canal window).
    PipeFill {
        /// The input-body edge the strip rolls along.
        source_edge: EdgeHandle,
    },
    /// Sphere window at a gap corner.
    SphereFill {
        /// The input-body vertex at the sphere center, when the corner
        /// maps back to one (rim-split midpoints do not).
        source_vertex: Option<VertexHandle>,
    },
}

/// Bidirectional handle mapping between the outer BRep and the cloned inner skin.
///
/// All maps are outer-handle -> inner-handle.
#[derive(Debug, Clone)]
pub struct SkinMap {
    /// Vertex mapping: outer -> inner.
    pub vertices: HashMap<VertexHandle, VertexHandle>,
    /// Edge mapping: outer -> inner.
    pub edges: HashMap<EdgeHandle, EdgeHandle>,
    /// Half-edge mapping: outer -> inner.
    pub half_edges: HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    /// Face mapping: outer -> inner.
    pub faces: HashMap<FaceHandle, FaceHandle>,
}

impl SkinMap {
    /// Build a reverse face map (inner -> outer).
    pub fn reverse_face_map(&self) -> HashMap<FaceHandle, FaceHandle> {
        self.faces
            .iter()
            .map(|(&outer, &inner)| (inner, outer))
            .collect()
    }

    /// Build a reverse vertex map (inner -> outer).
    pub fn reverse_vertex_map(&self) -> HashMap<VertexHandle, VertexHandle> {
        self.vertices
            .iter()
            .map(|(&outer, &inner)| (inner, outer))
            .collect()
    }

    /// Build a reverse edge map (inner -> outer).
    pub fn reverse_edge_map(&self) -> HashMap<EdgeHandle, EdgeHandle> {
        self.edges
            .iter()
            .map(|(&outer, &inner)| (inner, outer))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_maps_are_correct() {
        let map = SkinMap {
            vertices: HashMap::new(),
            edges: HashMap::new(),
            half_edges: HashMap::new(),
            faces: HashMap::new(),
        };
        assert!(map.reverse_face_map().is_empty());
        assert!(map.reverse_vertex_map().is_empty());
        assert!(map.reverse_edge_map().is_empty());
    }
}
