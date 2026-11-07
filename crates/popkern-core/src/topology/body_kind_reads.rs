//! The reads a body's kind is decided from, without allocating per
//! edge or vertex: every operation's entry gate reads the kind
//! ([`BRep::require_kind`]), so it costs a pass over the arenas and no
//! more (a 255-face slotted disc read its kind in 0.88 ms through
//! per-edge ring vectors and a hashed fan union, five whole-body
//! validations' worth, at every boolean).

use slotmap::SecondaryMap;

use super::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

/// An edge's face uses, counted to three.
#[derive(Debug, Clone, Copy)]
pub(crate) struct RingSummary {
    /// Whether the edge and its first fin are in the body.
    pub present: bool,
    /// Face uses on the ring (fins in a loop), counted up to 3.
    pub count: usize,
    /// The first two face uses, in ring order.
    pub fins: [HalfEdgeHandle; 2],
}

impl BRep {
    /// `edge`'s ring as [`BRep::radial_fins`] reads it, without
    /// collecting it: the walk stops at the third face use.
    pub(crate) fn ring_summary(&self, edge: EdgeHandle) -> RingSummary {
        let mut out = RingSummary {
            present: false,
            count: 0,
            fins: [HalfEdgeHandle::default(); 2],
        };
        let Some(first) = self
            .edges
            .get(edge)
            .map(|e| e.half_edge)
            .filter(|h| self.half_edges.contains_key(*h))
        else {
            return out;
        };
        out.present = true;
        let mut cursor = first;
        for _ in 0..self.half_edges.len() {
            if self.half_edges.get(cursor).is_some_and(|f| f.has_loop()) {
                if out.count < 2 {
                    out.fins[out.count] = cursor;
                }
                out.count += 1;
                if out.count > 2 {
                    break;
                }
            }
            match self.radial_next(cursor) {
                Some(next) if next != first && next != cursor => cursor = next,
                _ => break,
            }
        }
        out
    }

    /// Whether the face uses about every vertex form one fan, for a
    /// body whose rings hold at most two face uses. A fan is the face
    /// uses leaving the vertex joined through loops (a use and its
    /// loop predecessor's edge) and edge rings, as
    /// [`BRep::vertex_uses_all`] reads them: here one walk from some use
    /// of each vertex must reach every use leaving it, and a vertex a
    /// loop of no fins stands at has no other use.
    pub(super) fn one_fan_per_vertex(&self) -> bool {
        let mut leaving: SecondaryMap<VertexHandle, (usize, HalfEdgeHandle)> = SecondaryMap::new();
        for (h, fin) in &self.half_edges {
            if !self.loops.contains_key(fin.loop_) {
                continue;
            }
            match leaving.get_mut(fin.origin) {
                Some(entry) => entry.0 += 1,
                None => {
                    leaving.insert(fin.origin, (1, h));
                }
            }
        }
        // A vertex loop is a use of its own: its vertex has no other.
        let mut lone: SecondaryMap<VertexHandle, ()> = SecondaryMap::new();
        for (_, l) in &self.loops {
            if self.half_edges.contains_key(l.half_edge) {
                continue;
            }
            if let Some(v) = l.vertex
                && (leaving.contains_key(v) || lone.insert(v, ()).is_some())
            {
                return false;
            }
        }
        let mut seen: SecondaryMap<HalfEdgeHandle, ()> = SecondaryMap::new();
        let mut stack = Vec::new();
        for (vertex, &(count, start)) in &leaving {
            // The face use about `vertex` a ring fin stands for: itself
            // when it leaves the vertex, its loop successor when it
            // arrives.
            let leaving_use = |r: HalfEdgeHandle| -> Option<HalfEdgeHandle> {
                let fin = self.half_edges.get(r)?;
                let u = if fin.origin == vertex { r } else { fin.next };
                let use_ = self.half_edges.get(u)?;
                (use_.origin == vertex && self.loops.contains_key(use_.loop_)).then_some(u)
            };
            stack.clear();
            stack.push(start);
            seen.insert(start, ());
            let mut reached = 1;
            while let Some(h) = stack.pop() {
                let Some(fin) = self.half_edges.get(h) else {
                    continue;
                };
                let prev_edge = self.half_edges.get(fin.prev).map(|p| p.edge);
                for edge in [Some(fin.edge), prev_edge].into_iter().flatten() {
                    let ring = self.ring_summary(edge);
                    for &r in &ring.fins[..ring.count.min(2)] {
                        if let Some(u) = leaving_use(r)
                            && seen.insert(u, ()).is_none()
                        {
                            reached += 1;
                            stack.push(u);
                        }
                    }
                }
            }
            if reached != count {
                return false;
            }
        }
        true
    }
}
