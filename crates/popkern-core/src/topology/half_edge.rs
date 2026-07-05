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
    /// UV parameter curve on the parent face's surface.
    /// Encoded as `Point3 { x: u, y: v, z: 0 }`. `None` for edges
    /// without surface geometry. Stored as a [`Pcurve`] so analytic
    /// variants (linear, etc.) can be carried exactly.
    ///
    /// DIRECTION: the intended convention is traversal-aligned
    /// (pcurve start maps to this half-edge's `origin`), which
    /// `pcurve_recompute` produces. Primitive constructors predate
    /// the convention and may store edge-aligned pcurves, so
    /// consumers that need direction must verify via the
    /// origin/destination endpoint check until producers are
    /// normalized (PK-100).
    ///
    /// [`Pcurve`]: crate::geometry::pcurve::Pcurve
    pub pcurve: Option<crate::geometry::pcurve::Pcurve>,
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn half_edge_pcurve_is_populated_for_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        for (_, he) in &brep.half_edges {
            assert!(
                he.pcurve.is_some(),
                "Every cube half-edge should have a pcurve"
            );
        }
    }
}
