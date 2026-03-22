use crate::geometry::Point3;
use crate::topology::HalfEdgeHandle;
use serde::{Deserialize, Serialize};

/// A topological vertex — a point in the BRep model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Vertex {
    /// The geometric position.
    pub position: Point3,
    /// Any outgoing half-edge (set by mev, invalid before first edge).
    pub outgoing: HalfEdgeHandle,
}
