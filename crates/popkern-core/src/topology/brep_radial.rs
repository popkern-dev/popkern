//! Radial-ring accessors for [`BRep`]: every read of how fins share an
//! edge goes through these, so callers state whether they need a
//! manifold pair ([`BRep::mate`], [`BRep::require_mate`]) or walk any
//! ring ([`BRep::radial_next`], [`BRep::radial_fins`]).
//!
//! An edge's fins form its radial ring: none for a wire edge, one for a
//! free (sheet) edge, two for a manifold edge or a seam, more for a
//! non-manifold edge. N > 2 rings are currently threaded through
//! [`BRep::nonmanifold_radial`] by welders that merge coincident edges
//! (the sheet-sheet union welder in `popkern-bool`).

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

impl BRep {
    /// The other fin of `he`'s edge when the edge has exactly two fins,
    /// `None` for a free edge, a non-manifold edge, or an unknown
    /// handle.
    pub fn mate(&self, he: HalfEdgeHandle) -> Option<HalfEdgeHandle> {
        if self.nonmanifold_radial.contains_key(&he) {
            return None;
        }
        self.half_edges.get(he)?.mate()
    }

    /// [`BRep::mate`] for code that requires a manifold edge.
    ///
    /// # Errors
    ///
    /// [`CoreError::NoMate`] naming the edge when it does not have
    /// exactly two fins; [`CoreError::InvalidHandle`] for an unknown
    /// fin.
    pub fn require_mate(&self, he: HalfEdgeHandle) -> Result<HalfEdgeHandle, CoreError> {
        if let Some(mate) = self.mate(he) {
            return Ok(mate);
        }
        let edge = self
            .half_edges
            .get(he)
            .ok_or_else(|| CoreError::InvalidHandle(format!("fin {he:?}")))?
            .edge;
        Err(CoreError::NoMate { edge })
    }

    /// Link `a` and `b` as the two fins of a manifold edge.
    pub fn set_mates(&mut self, a: HalfEdgeHandle, b: HalfEdgeHandle) {
        if let Some(fin) = self.half_edges.get_mut(a) {
            fin.set_mate(b);
        }
        if let Some(fin) = self.half_edges.get_mut(b) {
            fin.set_mate(a);
        }
    }

    /// The next fin around `he`'s edge: its mate on a manifold edge,
    /// the next fin of the ring on a non-manifold edge, `he` itself on
    /// a free edge. `None` only for an unknown handle.
    pub fn radial_next(&self, he: HalfEdgeHandle) -> Option<HalfEdgeHandle> {
        if let Some(next) = self.nonmanifold_radial.get(&he) {
            return Some(*next);
        }
        let fin = self.half_edges.get(he)?;
        Some(
            fin.mate()
                .filter(|m| self.half_edges.contains_key(*m))
                .unwrap_or(he),
        )
    }

    /// Every fin of `edge`'s ring, starting at the edge's first fin and
    /// in ring order. Empty for a wire edge or an unknown edge.
    pub fn radial_fins(&self, edge: EdgeHandle) -> Vec<HalfEdgeHandle> {
        let Some(first) = self
            .edges
            .get(edge)
            .map(|e| e.half_edge)
            .filter(|h| self.half_edges.get(*h).is_some_and(|f| f.has_loop()))
        else {
            return Vec::new();
        };
        let mut fins = vec![first];
        let mut cursor = first;
        // A ring never holds more fins than the body has.
        for _ in 0..self.half_edges.len() {
            match self.radial_next(cursor) {
                Some(next) if next != first => {
                    fins.push(next);
                    cursor = next;
                }
                _ => break,
            }
        }
        fins
    }

    /// How many fins `edge`'s ring holds (0 for a wire edge).
    pub fn edge_fin_count(&self, edge: EdgeHandle) -> usize {
        self.radial_fins(edge).len()
    }
}

impl BRep {
    /// The vertex a fin runs to: its mate's origin, or for a fin with
    /// no mate (a free edge, a non-manifold ring) its loop successor's
    /// origin. The mate is read first so the answer holds while loops
    /// are being wired (construction links mates before `next`).
    pub fn fin_end(&self, he: HalfEdgeHandle) -> Option<crate::topology::VertexHandle> {
        let fin = self.half_edges.get(he)?;
        if let Some(mate) = self.mate(he).and_then(|m| self.half_edges.get(m)) {
            return Some(mate.origin);
        }
        if !fin.has_loop() {
            return None;
        }
        self.half_edges.get(fin.next).map(|f| f.origin)
    }

    /// An edge's two vertices, in the direction of its first fin.
    pub fn edge_vertices(
        &self,
        edge: EdgeHandle,
    ) -> Option<(crate::topology::VertexHandle, crate::topology::VertexHandle)> {
        let he = self.edges.get(edge)?.half_edge;
        Some((self.half_edges.get(he)?.origin, self.fin_end(he)?))
    }
}

#[cfg(test)]
impl BRep {
    /// Test fixture: give `edge` a third fin (a copy of its first fin's
    /// loop use) and thread the three into one radial ring, so the edge
    /// reads as non-manifold everywhere the ring is asked. Returns the
    /// new fin. The body is not valid afterwards (the copy is in no
    /// loop cycle); tests use it to check that classifiers see N > 2.
    pub(crate) fn thread_test_third_fin(&mut self, edge: EdgeHandle) -> HalfEdgeHandle {
        let base = self.edges[edge].half_edge;
        let mate = self.mate(base).unwrap_or(base);
        let copy = self.half_edges[base].clone();
        let third = self.half_edges.insert(copy);
        self.nonmanifold_radial.insert(base, mate);
        self.nonmanifold_radial.insert(mate, third);
        self.nonmanifold_radial.insert(third, base);
        third
    }
}
