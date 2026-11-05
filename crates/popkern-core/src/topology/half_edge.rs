use crate::topology::{EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};
use serde::{Deserialize, Serialize};

/// A fin: one loop's directed use of an edge.
///
/// The fins of an edge form its radial ring (none for a wire edge, one
/// for a free edge, two for a manifold edge, more for a non-manifold
/// one). Read the ring through [`HalfEdge::mate`] and the
/// [`BRep`](crate::topology::BRep) radial accessors, never through the
/// link fields: their encoding is the ring's, not the manifold pair's.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HalfEdge {
    /// Vertex at the start of this half-edge.
    pub origin: VertexHandle,
    /// The opposite-direction half-edge of a manifold pair, null when
    /// the fin has no mate. Read it with [`HalfEdge::mate`] and write it
    /// with [`HalfEdge::set_mate`] / [`HalfEdge::clear_mate`].
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

impl HalfEdge {
    /// The other fin of this fin's edge when the edge has exactly two
    /// fins (a manifold edge, or a seam), `None` for a free edge, a
    /// non-manifold edge or a fin not yet linked.
    ///
    /// Code that requires a manifold edge treats `None` as its own
    /// typed failure (or uses [`BRep::require_mate`](crate::topology::BRep::require_mate));
    /// code walking a non-manifold edge uses the radial accessors.
    pub fn mate(&self) -> Option<HalfEdgeHandle> {
        use slotmap::Key;
        (!self.twin.is_null()).then_some(self.twin)
    }

    /// [`HalfEdge::mate`] for code that requires a manifold edge.
    ///
    /// # Errors
    ///
    /// [`CoreError::NoMate`](crate::error::CoreError::NoMate) naming
    /// this fin's edge when it has no mate.
    pub fn require_mate(&self) -> Result<HalfEdgeHandle, crate::error::CoreError> {
        self.mate()
            .ok_or(crate::error::CoreError::NoMate { edge: self.edge })
    }

    /// Make `other` this fin's mate (one direction: the caller links
    /// the other fin back, usually through
    /// [`BRep::set_mates`](crate::topology::BRep::set_mates)).
    pub fn set_mate(&mut self, other: HalfEdgeHandle) {
        self.twin = other;
    }

    /// Whether this fin lies in a loop (a face's use of its edge); the
    /// half-edges of a wire edge do not.
    pub fn has_loop(&self) -> bool {
        use slotmap::Key;
        !self.loop_.is_null()
    }

    /// Leave this fin alone on its edge (a free edge).
    pub fn clear_mate(&mut self) {
        self.twin = HalfEdgeHandle::default();
    }
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
