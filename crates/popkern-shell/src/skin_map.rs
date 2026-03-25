//! Handle mapping between outer and inner skin BReps.

use std::collections::HashMap;

use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Bidirectional handle mapping between the outer BRep and the cloned inner skin.
///
/// All maps are outer-handle → inner-handle.
#[derive(Debug, Clone)]
pub struct SkinMap {
    /// Vertex mapping: outer → inner.
    pub vertices: HashMap<VertexHandle, VertexHandle>,
    /// Edge mapping: outer → inner.
    pub edges: HashMap<EdgeHandle, EdgeHandle>,
    /// Half-edge mapping: outer → inner.
    pub half_edges: HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    /// Face mapping: outer → inner.
    pub faces: HashMap<FaceHandle, FaceHandle>,
}
