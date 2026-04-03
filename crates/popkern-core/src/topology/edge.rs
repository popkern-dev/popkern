use crate::geometry::CurveOnEdge;
use crate::naming::NameId;
use crate::topology::HalfEdgeHandle;
use serde::{Deserialize, Serialize};

/// A topological edge — an undirected edge with two half-edges.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    /// One of the two half-edges.
    pub half_edge: HalfEdgeHandle,
    /// Geometric curve.
    pub curve: Option<CurveOnEdge>,
    /// Stable topological name for this edge.
    #[serde(default)]
    pub name: Option<NameId>,
    /// Per-entity tolerance (accuracy of this edge's curve geometry).
    /// `None` means use the BRep-level default tolerance.
    #[serde(default)]
    pub tolerance: Option<f64>,
}
