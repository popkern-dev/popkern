use crate::topology::{FaceHandle, HalfEdgeHandle};
use serde::{Deserialize, Serialize};

/// A boundary loop of a face (outer boundary or inner hole).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Loop {
    /// Any half-edge in this loop.
    pub half_edge: HalfEdgeHandle,
    /// Parent face.
    pub face: FaceHandle,
}
