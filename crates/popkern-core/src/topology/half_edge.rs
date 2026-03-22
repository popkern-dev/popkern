use crate::topology::{EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};
use serde::{Deserialize, Serialize};

/// A directed half-edge in the BRep model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HalfEdge {
    /// Vertex at the start of this half-edge.
    pub origin: VertexHandle,
    /// The opposite-direction half-edge.
    pub twin: HalfEdgeHandle,
    /// Next half-edge in the face loop.
    pub next: HalfEdgeHandle,
    /// Previous half-edge in the face loop.
    pub prev: HalfEdgeHandle,
    /// Parent edge.
    pub edge: EdgeHandle,
    /// Parent loop.
    pub loop_: LoopHandle,
}
