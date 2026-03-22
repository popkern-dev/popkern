use crate::topology::HalfEdgeHandle;
use serde::{Deserialize, Serialize};

/// A topological edge — an undirected edge with two half-edges.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    /// One of the two half-edges.
    pub half_edge: HalfEdgeHandle,
    /// Geometric curve. Placeholder until geometry is wired in (Task 18).
    pub curve: Option<()>, // TODO: Replace with Option<CurveOnEdge>
}
