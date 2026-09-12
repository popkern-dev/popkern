//! Result types for the face/edge split operators (split from
//! `split.rs` under the 500-line structural budget).

use crate::topology::{EdgeHandle, FaceHandle, VertexHandle};

/// Result of splitting a face along a trim curve.
#[derive(Debug, Clone)]
pub struct SplitFaceResult {
    /// One side of the split (the original face).
    pub face_a: FaceHandle,
    /// The other side of the split (the new face from mef).
    pub face_b: FaceHandle,
    /// The newly created edges along the trim curve.
    pub new_edges: Vec<EdgeHandle>,
}

/// Result of splitting an edge at a parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SplitEdgeResult {
    /// The newly created vertex at the split point.
    pub vertex: VertexHandle,
    /// The original edge (start -> new vertex).
    pub edge_a: EdgeHandle,
    /// The new edge (new vertex -> end).
    pub edge_b: EdgeHandle,
}
