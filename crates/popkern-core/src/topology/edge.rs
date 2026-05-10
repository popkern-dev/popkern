use crate::geometry::CurveOnEdge;
use crate::naming::NameId;
use crate::topology::{HalfEdgeHandle, TopologyRole};
use serde::{Deserialize, Serialize};

/// Which representation is authoritative for edge geometry.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CanonicalForm {
    /// 3D curve is authoritative (default).
    #[default]
    Canonical3D,
    /// Parametric curve is authoritative (from 2D sketch).
    CanonicalParametric,
}

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
    /// Which geometric representation is authoritative for this edge.
    #[serde(default)]
    pub canonical: CanonicalForm,
    /// Whether this edge belongs to production topology or retained construction residue.
    #[serde(default)]
    pub role: TopologyRole,
}
