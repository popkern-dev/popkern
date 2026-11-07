use crate::topology::{FaceHandle, HalfEdgeHandle, VertexHandle};
use serde::{Deserialize, Serialize};

/// A boundary loop of a face (outer boundary or inner hole).
///
/// A loop of no fins stands at a lone vertex, which it names: the empty
/// loop the first Euler step leaves (`mvfs`), and a vertex touching a
/// face's interior (a vertex loop, PK-1116), the face holding the
/// vertex as topology.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Loop {
    /// Any half-edge in this loop (null for a loop of no fins).
    pub half_edge: HalfEdgeHandle,
    /// Parent face.
    pub face: FaceHandle,
    /// The lone vertex a loop of no fins stands at; `None` for a loop
    /// of fins.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vertex: Option<VertexHandle>,
}

impl Loop {
    /// A loop of `face` starting at `half_edge` (null while its fins
    /// are being linked), naming no lone vertex.
    pub fn new(half_edge: HalfEdgeHandle, face: FaceHandle) -> Self {
        Self {
            half_edge,
            face,
            vertex: None,
        }
    }
}
