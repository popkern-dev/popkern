//! Non-manifold radial-cycle helper for [`BRep`].
//!
//! [`BRep::radial_next`] walks around an edge honouring the
//! optional non-manifold radial cycle stored in
//! [`BRep::nonmanifold_radial`]. The cycle is populated by welders
//! that merge multiple co-incident edges into one — e.g. the
//! sheet-sheet-union welder in `popkern-bool` — so downstream code
//! can visit every face-use of a non-manifold edge without needing
//! to know whether the edge is 2-manifold or not.

use crate::topology::{BRep, HalfEdgeHandle};

impl BRep {
    /// Next half-edge around the parent edge, honouring non-manifold
    /// radial cycles.
    ///
    /// For 2-manifold edges the result is `he.twin` — the only
    /// other half-edge of the pair. For non-manifold edges
    /// (N > 2 half-edges sharing an edge, e.g. a transversal
    /// sheet-sheet crossing), [`BRep::nonmanifold_radial`] overrides
    /// the fallback and returns the next half-edge in the cyclic
    /// chain `h1 → h2 → ... → hN → h1`. A full traversal visits
    /// every half-edge of the group exactly once before returning
    /// to the starting handle.
    ///
    /// Returns `None` only if `he` is not a valid half-edge handle.
    pub fn radial_next(&self, he: HalfEdgeHandle) -> Option<HalfEdgeHandle> {
        if let Some(next) = self.nonmanifold_radial.get(&he) {
            return Some(*next);
        }
        self.half_edges.get(he).map(|h| h.twin)
    }
}
