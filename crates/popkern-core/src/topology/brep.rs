//! BRep — the top-level boundary representation container.

use crate::tolerance::Tolerance;
use crate::topology::*;
use serde::{Deserialize, Serialize};
use slotmap::SlotMap;

/// A boundary representation model.
///
/// Owns all topological entities in typed arenas (SlotMaps).
/// All topology modifications go through Euler operators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    /// All vertices.
    pub vertices: SlotMap<VertexHandle, Vertex>,
    /// All half-edges.
    pub half_edges: SlotMap<HalfEdgeHandle, HalfEdge>,
    /// All edges.
    pub edges: SlotMap<EdgeHandle, Edge>,
    /// All loops.
    pub loops: SlotMap<LoopHandle, Loop>,
    /// All faces.
    pub faces: SlotMap<FaceHandle, Face>,
    /// All shells.
    pub shells: SlotMap<ShellHandle, Shell>,
    /// All solids.
    pub solids: SlotMap<SolidHandle, Solid>,
    /// Geometric tolerance for this model.
    pub tolerance: Tolerance,
}

impl BRep {
    /// Create a new empty BRep with the given tolerance.
    pub fn new(tolerance: Tolerance) -> Self {
        Self {
            vertices: SlotMap::with_key(),
            half_edges: SlotMap::with_key(),
            edges: SlotMap::with_key(),
            loops: SlotMap::with_key(),
            faces: SlotMap::with_key(),
            shells: SlotMap::with_key(),
            solids: SlotMap::with_key(),
            tolerance,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_brep_is_empty() {
        let brep = BRep::new(Tolerance::default());
        assert_eq!(brep.vertices.len(), 0);
        assert_eq!(brep.half_edges.len(), 0);
        assert_eq!(brep.edges.len(), 0);
        assert_eq!(brep.loops.len(), 0);
        assert_eq!(brep.faces.len(), 0);
        assert_eq!(brep.shells.len(), 0);
        assert_eq!(brep.solids.len(), 0);
    }

    #[test]
    fn brep_has_tolerance() {
        let tol = Tolerance::default().with_linear(1e-3);
        let brep = BRep::new(tol);
        assert_eq!(brep.tolerance.linear, 1e-3);
    }

    #[test]
    fn brep_serializes_round_trip() {
        let brep = BRep::new(Tolerance::default());
        let json = serde_json::to_string(&brep).unwrap();
        let brep2: BRep = serde_json::from_str(&json).unwrap();
        assert_eq!(brep2.vertices.len(), 0);
    }
}
