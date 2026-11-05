//! Radial-ring accessors for [`BRep`]: every read of how fins share an
//! edge goes through these, so callers state whether they need a
//! manifold pair ([`BRep::mate`], [`BRep::require_mate`]) or walk any
//! ring ([`BRep::radial_next`], [`BRep::radial_fins`]).
//!
//! An edge's fins form its radial ring: none for a wire edge, one for a
//! free (sheet) edge, two for a manifold edge or a seam, more for a
//! non-manifold edge. Each fin links to its ring neighbours
//! (`radial_next`, `radial_prev`); a lone fin has null links.

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

impl BRep {
    /// The other fin of `he`'s edge when the edge has exactly two fins,
    /// `None` for a free edge, a non-manifold edge, or an unknown
    /// handle.
    pub fn mate(&self, he: HalfEdgeHandle) -> Option<HalfEdgeHandle> {
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
        let fin = self.half_edges.get(he)?;
        Some(
            Some(fin.ring.next)
                .filter(|n| self.half_edges.contains_key(*n))
                .unwrap_or(he),
        )
    }

    /// The previous fin around `he`'s edge (`he` itself on a free
    /// edge). `None` only for an unknown handle.
    pub fn radial_prev(&self, he: HalfEdgeHandle) -> Option<HalfEdgeHandle> {
        let fin = self.half_edges.get(he)?;
        Some(
            Some(fin.ring.prev)
                .filter(|p| self.half_edges.contains_key(*p))
                .unwrap_or(he),
        )
    }

    /// Link `fins` into one radial ring in the given order (each fin's
    /// `radial_next` is the following one, cyclically). One fin is left
    /// alone (a free edge, null links), two become mates. The caller
    /// keeps the fins' `edge` fields and the edge's `half_edge` in step.
    pub fn link_ring(&mut self, fins: &[HalfEdgeHandle]) {
        let n = fins.len();
        for (i, &fin) in fins.iter().enumerate() {
            let (next, prev) = if n < 2 {
                (HalfEdgeHandle::default(), HalfEdgeHandle::default())
            } else {
                (fins[(i + 1) % n], fins[(i + n - 1) % n])
            };
            if let Some(f) = self.half_edges.get_mut(fin) {
                f.ring = crate::topology::RadialLinks { next, prev };
            }
        }
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
    /// Test fixture: add copies of `edge`'s first fin (its loop use)
    /// until its ring holds three fins, so the edge reads as
    /// non-manifold everywhere the ring is asked. Returns the last fin
    /// added. The body is not valid afterwards (the copy is in no
    /// loop cycle); tests use it to check that classifiers see N > 2.
    pub(crate) fn thread_test_third_fin(&mut self, edge: EdgeHandle) -> HalfEdgeHandle {
        let base = self.edges[edge].half_edge;
        let mut ring = self.radial_fins(edge);
        let mut last = base;
        while ring.len() < 3 {
            let copy = self.half_edges[base].clone();
            last = self.half_edges.insert(copy);
            ring.push(last);
        }
        self.link_ring(&ring);
        last
    }
}
