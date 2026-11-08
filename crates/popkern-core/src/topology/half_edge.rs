use crate::topology::{EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};
use serde::{Deserialize, Serialize};

/// A fin's neighbours in its edge's radial ring: `next` is the following
/// fin counter-clockwise about the edge's direction, `prev` the one
/// before. Both are the mate on a manifold edge and null for a fin alone
/// on its edge (a free edge).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct RadialLinks {
    /// The next fin around the edge.
    pub next: HalfEdgeHandle,
    /// The previous fin around the edge.
    pub prev: HalfEdgeHandle,
}

impl RadialLinks {
    /// A fin alone on its edge.
    pub fn alone() -> Self {
        Self::default()
    }

    /// One of the two fins of a manifold edge, `mate` the other (null:
    /// alone).
    pub fn pair(mate: HalfEdgeHandle) -> Self {
        Self {
            next: mate,
            prev: mate,
        }
    }
}

/// A fin: one loop's directed use of an edge.
///
/// The fins of an edge form its radial ring (none for a wire edge, one
/// for a free edge, two for a manifold edge, more for a non-manifold
/// one). Read the ring through [`HalfEdge::mate`] and the
/// [`BRep`](crate::topology::BRep) radial accessors, never through the
/// link fields: their encoding is the ring's, not the manifold pair's.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(from = "HalfEdgeRecord")]
pub struct HalfEdge {
    /// Vertex at the start of this half-edge.
    pub origin: VertexHandle,
    /// This fin's neighbours in its edge's radial ring. Read it with
    /// [`HalfEdge::mate`] or the [`BRep`](crate::topology::BRep) radial
    /// accessors; write it with
    /// [`BRep::link_ring`](crate::topology::BRep::link_ring),
    /// [`HalfEdge::set_mate`] or [`HalfEdge::clear_mate`].
    pub ring: RadialLinks,
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
        // A ring of two links each fin both ways to the other; a longer
        // ring has distinct neighbours, a lone fin none.
        (!self.ring.next.is_null() && self.ring.next == self.ring.prev).then_some(self.ring.next)
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
        self.ring = RadialLinks::pair(other);
    }

    /// Whether this fin lies in a loop (a face's use of its edge); the
    /// half-edges of a wire edge do not.
    pub fn has_loop(&self) -> bool {
        use slotmap::Key;
        !self.loop_.is_null()
    }

    /// Leave this fin alone on its edge (a free edge).
    pub fn clear_mate(&mut self) {
        self.ring = RadialLinks::alone();
    }
}

/// The serialized form of a fin: the ring links, or the `twin` of
/// bodies written before the radial ring (a manifold pair, read as a
/// ring of two).
#[derive(Deserialize)]
struct HalfEdgeRecord {
    origin: VertexHandle,
    #[serde(default)]
    ring: Option<RadialLinks>,
    #[serde(default)]
    twin: Option<HalfEdgeHandle>,
    next: HalfEdgeHandle,
    prev: HalfEdgeHandle,
    edge: EdgeHandle,
    loop_: LoopHandle,
    #[serde(default)]
    pcurve: Option<crate::geometry::pcurve::Pcurve>,
}

impl From<HalfEdgeRecord> for HalfEdge {
    fn from(r: HalfEdgeRecord) -> Self {
        let ring = r
            .ring
            .unwrap_or_else(|| RadialLinks::pair(r.twin.unwrap_or_default()));
        HalfEdge {
            origin: r.origin,
            ring,
            next: r.next,
            prev: r.prev,
            edge: r.edge,
            loop_: r.loop_,
            pcurve: r.pcurve,
        }
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
