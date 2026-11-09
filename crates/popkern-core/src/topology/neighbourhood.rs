//! What a body is like around one of its vertices or edges: one
//! authority on manifoldness (PK-1129).
//!
//! A body is a manifold at a vertex when the faces around the vertex
//! close up into one disc: the half-edges leaving it form one cycle
//! under `h -> mate(h).next`. Several cycles are several cones of
//! material meeting at the point, non-manifold. A face may appear more
//! than once in that one cycle (a face touching itself at the vertex, a
//! pinch use: the top face of a plate a spherical dimple touches at its
//! edge), which is still a manifold. An open cycle is the vertex on a
//! sheet's free boundary. Whether two of a body's vertices happen to
//! coincide says nothing here: that is identity, not manifoldness.
//!
//! At an edge: a free edge has one fin, a manifold edge two, an edge of
//! a general body's non-manifold ring more.

use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// The body around a vertex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VertexNeighbourhood {
    /// One closed fan: manifold. The faces it holds more than once.
    Disc {
        /// Faces met more than once round the vertex (pinch uses).
        pinch_uses: Vec<FaceHandle>,
    },
    /// One open fan: the vertex on a sheet's free boundary.
    HalfDisc {
        /// Faces met more than once.
        pinch_uses: Vec<FaceHandle>,
    },
    /// Several fans: cones of material meeting at the point.
    Fans {
        /// How many.
        count: usize,
    },
    /// No face uses the vertex (an acorn, a wire's end, a vertex loop).
    Lower,
}

impl VertexNeighbourhood {
    /// Whether the body is a manifold (with or without boundary) there.
    pub fn is_manifold(&self) -> bool {
        matches!(self, Self::Disc { .. } | Self::HalfDisc { .. })
    }
}

/// The body around an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeNeighbourhood {
    /// No face uses it (a wire edge).
    Wire,
    /// One fin: a sheet's free edge.
    Free,
    /// Two fins: manifold.
    Manifold,
    /// A ring of more than two fins: non-manifold.
    Ring {
        /// How many.
        fins: usize,
    },
}

impl BRep {
    /// The body around vertex `v` (see the module docs).
    pub fn vertex_neighbourhood(&self, v: VertexHandle) -> VertexNeighbourhood {
        let mut leaving: Vec<HalfEdgeHandle> = self
            .half_edges
            .iter()
            .filter(|(_, h)| h.origin == v)
            .map(|(k, _)| k)
            .collect();
        leaving.sort();
        if leaving.is_empty() {
            return VertexNeighbourhood::Lower;
        }
        let mut seen: Vec<HalfEdgeHandle> = Vec::with_capacity(leaving.len());
        let mut fans: Vec<(Vec<FaceHandle>, bool)> = Vec::new();
        for &start in &leaving {
            if seen.contains(&start) {
                continue;
            }
            // Walk the fan both ways from `start`: forward under
            // h -> mate(h).next, backward under h -> mate(prev(h)).
            let mut faces = Vec::new();
            let mut open = false;
            let mut cur = start;
            let mut closed = false;
            for _ in 0..=leaving.len() {
                seen.push(cur);
                if let Some(f) = self.face_of(cur) {
                    faces.push(f);
                }
                match self.fan_next(cur, v) {
                    Some(n) if n == start => {
                        closed = true;
                        break;
                    }
                    Some(n) if !seen.contains(&n) => cur = n,
                    Some(_) => break,
                    None => {
                        open = true;
                        break;
                    }
                }
            }
            if !closed {
                // The other way round to the fan's other end.
                let mut back = start;
                for _ in 0..=leaving.len() {
                    match self.fan_prev(back, v) {
                        Some(p) if !seen.contains(&p) => {
                            seen.push(p);
                            if let Some(f) = self.face_of(p) {
                                faces.push(f);
                            }
                            back = p;
                        }
                        Some(_) => break,
                        None => {
                            open = true;
                            break;
                        }
                    }
                }
            }
            fans.push((faces, open));
        }
        if fans.len() > 1 {
            return VertexNeighbourhood::Fans { count: fans.len() };
        }
        let (faces, open) = fans.swap_remove(0);
        let mut sorted = faces.clone();
        sorted.sort();
        let mut pinch_uses: Vec<FaceHandle> = sorted
            .windows(2)
            .filter(|w| w[0] == w[1])
            .map(|w| w[0])
            .collect();
        pinch_uses.dedup();
        if open {
            VertexNeighbourhood::HalfDisc { pinch_uses }
        } else {
            VertexNeighbourhood::Disc { pinch_uses }
        }
    }

    /// The body around edge `e`.
    pub fn edge_neighbourhood(&self, e: EdgeHandle) -> EdgeNeighbourhood {
        match self.ring_summary(e).count {
            0 => EdgeNeighbourhood::Wire,
            1 => EdgeNeighbourhood::Free,
            2 => EdgeNeighbourhood::Manifold,
            n => EdgeNeighbourhood::Ring { fins: n },
        }
    }

    fn face_of(&self, h: HalfEdgeHandle) -> Option<FaceHandle> {
        let l = self.half_edges.get(h)?.loop_;
        self.loops.get(l).map(|l| l.face)
    }

    /// The next half-edge leaving `v` round the fan: across `h`'s edge,
    /// then on along the other face's loop.
    fn fan_next(&self, h: HalfEdgeHandle, v: VertexHandle) -> Option<HalfEdgeHandle> {
        let m = self.mate(h)?;
        let n = self.half_edges.get(m)?.next;
        (self.half_edges.get(n)?.origin == v).then_some(n)
    }

    /// The previous one: the half-edge arriving at `v` before `h` in its
    /// loop, across its edge.
    fn fan_prev(&self, h: HalfEdgeHandle, v: VertexHandle) -> Option<HalfEdgeHandle> {
        let l = self.half_edges.get(h)?.loop_;
        let start = self.loops.get(l)?.half_edge;
        // The loop's half-edge whose next is `h`.
        let mut p = start;
        for _ in 0..self.half_edges.len() {
            let he = self.half_edges.get(p)?;
            if he.next == h {
                break;
            }
            p = he.next;
            if p == start {
                return None;
            }
        }
        let m = self.mate(p)?;
        (self.half_edges.get(m)?.origin == v).then_some(m)
    }
}

#[cfg(test)]
#[path = "neighbourhood_tests.rs"]
mod tests;
