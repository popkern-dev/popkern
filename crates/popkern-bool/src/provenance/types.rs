//! Provenance types for tracking input→output entity correspondence.
//!
//! After a boolean operation, each output entity (face, edge, vertex)
//! has a provenance that explains where it came from:
//! - KEPT: directly from an input body (unmodified)
//! - SPLIT: part of an input entity that was split by intersection
//! - CREATED: new entity at an intersection (e.g. intersection edge)
//! - MERGED: result of merging coincident entities

use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

/// Which input body an entity came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceBody {
    /// From the first operand (A).
    A,
    /// From the second operand (B).
    B,
    /// Created at the intersection of A and B.
    Intersection,
}

/// How an output face relates to input faces.
#[derive(Debug, Clone)]
pub enum FaceProvenance {
    /// Face kept from an input body (possibly trimmed).
    Kept {
        /// Which body it came from.
        source: SourceBody,
        /// The original face handle in the source body.
        original_face: FaceHandle,
    },
    /// Face created at an intersection (e.g. a cap face from subtraction).
    Created {
        /// Face from body A that contributed.
        face_a: FaceHandle,
        /// Face from body B that contributed.
        face_b: FaceHandle,
    },
}

/// How an output edge relates to input edges.
#[derive(Debug, Clone)]
pub enum EdgeProvenance {
    /// Edge kept from an input body.
    Original {
        /// Which body.
        source: SourceBody,
        /// The original edge handle.
        original_edge: EdgeHandle,
    },
    /// Edge created at a face-face intersection.
    Intersection {
        /// Face from body A.
        face_a: FaceHandle,
        /// Face from body B.
        face_b: FaceHandle,
    },
    /// Edge created by splitting an original edge at an intersection point.
    Split {
        /// Which body the original edge belongs to.
        source: SourceBody,
        /// The original edge that was split.
        original_edge: EdgeHandle,
        /// Parameter on the original edge where the split occurred.
        split_parameter: f64,
    },
}

/// How an output vertex relates to input vertices.
#[derive(Debug, Clone)]
pub enum VertexProvenance {
    /// Vertex kept from an input body.
    Original {
        /// Which body.
        source: SourceBody,
        /// The original vertex handle.
        original_vertex: VertexHandle,
    },
    /// Vertex created at an intersection curve endpoint.
    IntersectionEndpoint,
    /// Vertex created by splitting an edge at an intersection point.
    EdgeSplit {
        /// Which body the edge belongs to.
        source: SourceBody,
        /// The edge that was split.
        edge: EdgeHandle,
        /// Parameter on the edge.
        parameter: f64,
    },
    /// Vertex from merging two coincident vertices.
    Merged {
        /// The vertices that were merged.
        sources: Vec<(SourceBody, VertexHandle)>,
    },
}

/// Complete entity map for a boolean result.
///
/// Maps output entity handles to their provenance. Not all entities
/// will have provenance (e.g. entities created during cleanup).
#[derive(Debug, Clone, Default)]
pub struct EntityMap {
    /// Face provenance.
    pub faces: Vec<(FaceHandle, FaceProvenance)>,
    /// Edge provenance.
    pub edges: Vec<(EdgeHandle, EdgeProvenance)>,
    /// Vertex provenance.
    pub vertices: Vec<(VertexHandle, VertexProvenance)>,
}

impl EntityMap {
    /// Look up the provenance of an output face.
    pub fn face_provenance(&self, fh: FaceHandle) -> Option<&FaceProvenance> {
        self.faces.iter().find(|(h, _)| *h == fh).map(|(_, p)| p)
    }

    /// Look up the provenance of an output edge.
    pub fn edge_provenance(&self, eh: EdgeHandle) -> Option<&EdgeProvenance> {
        self.edges.iter().find(|(h, _)| *h == eh).map(|(_, p)| p)
    }

    /// Look up the provenance of an output vertex.
    pub fn vertex_provenance(&self, vh: VertexHandle) -> Option<&VertexProvenance> {
        self.vertices.iter().find(|(h, _)| *h == vh).map(|(_, p)| p)
    }

    /// Number of tracked entities.
    pub fn total(&self) -> usize {
        self.faces.len() + self.edges.len() + self.vertices.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entity_map_lookup() {
        let mut map = EntityMap::default();
        assert_eq!(map.total(), 0);

        let fh = FaceHandle::default();
        map.faces.push((
            fh,
            FaceProvenance::Kept {
                source: SourceBody::A,
                original_face: fh,
            },
        ));
        assert_eq!(map.total(), 1);
        assert!(map.face_provenance(fh).is_some());
    }

    #[test]
    fn empty_map_returns_none() {
        let map = EntityMap::default();
        assert!(map.face_provenance(FaceHandle::default()).is_none());
        assert!(map.edge_provenance(EdgeHandle::default()).is_none());
        assert!(map.vertex_provenance(VertexHandle::default()).is_none());
    }
}
