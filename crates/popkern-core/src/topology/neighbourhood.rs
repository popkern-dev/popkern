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

/// The faces met in more than one run round a fan (`faces` in fan
/// order, cyclic when the fan closes). Consecutive uses of one face are
/// one run: the fan crossed a seam edge, the face continuing into
/// itself (a cylinder's lateral face either side of its seam at a rim
/// vertex, a sphere's at its pole), which is no pinch. A pinch is the
/// face coming back after others (a plate's top face either side of a
/// dimple touching its edge).
fn pinch_uses_of(faces: &[FaceHandle], open: bool) -> Vec<FaceHandle> {
    let n = faces.len();
    let mut runs: Vec<FaceHandle> = (0..n)
        .filter(|&i| i == 0 || faces[i] != faces[i - 1])
        .map(|i| faces[i])
        .collect();
    // A closed fan's last run continues into its first.
    if !open && runs.len() > 1 && runs.first() == runs.last() {
        runs.pop();
    }
    let mut sorted = runs;
    sorted.sort();
    let mut pinch: Vec<FaceHandle> = sorted
        .windows(2)
        .filter(|w| w[0] == w[1])
        .map(|w| w[0])
        .collect();
    pinch.dedup();
    pinch
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
        self.neighbourhood_of(v, &leaving)
    }

    /// The body around every vertex, in vertex order: one pass over the
    /// half-edges, where asking per vertex scans them all each time.
    pub fn vertex_neighbourhoods(&self) -> Vec<(VertexHandle, VertexNeighbourhood)> {
        let mut leaving: Vec<(VertexHandle, HalfEdgeHandle)> = self
            .half_edges
            .iter()
            .map(|(h, he)| (he.origin, h))
            .collect();
        leaving.sort();
        let mut out = Vec::with_capacity(self.vertices.len());
        let mut runs = leaving.chunk_by(|a, b| a.0 == b.0).peekable();
        for v in self.vertices.keys() {
            // Runs of a vertex no longer in the body are skipped.
            while runs.peek().is_some_and(|run| run[0].0 < v) {
                runs.next();
            }
            let hs: Vec<HalfEdgeHandle> = match runs.peek() {
                Some(run) if run[0].0 == v => runs
                    .next()
                    .map(|run| run.iter().map(|&(_, h)| h).collect())
                    .unwrap_or_default(),
                _ => Vec::new(),
            };
            out.push((v, self.neighbourhood_of(v, &hs)));
        }
        out
    }

    /// The neighbourhood of `v` from `leaving`, its leaving half-edges
    /// in handle order.
    fn neighbourhood_of(&self, v: VertexHandle, leaving: &[HalfEdgeHandle]) -> VertexNeighbourhood {
        if leaving.is_empty() {
            return VertexNeighbourhood::Lower;
        }
        let mut seen: Vec<HalfEdgeHandle> = Vec::with_capacity(leaving.len());
        let mut fans: Vec<(Vec<FaceHandle>, bool)> = Vec::new();
        for &start in leaving {
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
                // The other way round to the fan's other end, kept in
                // fan order (before the start).
                let mut back = start;
                let mut before = Vec::new();
                for _ in 0..=leaving.len() {
                    match self.fan_prev(back, v) {
                        Some(p) if !seen.contains(&p) => {
                            seen.push(p);
                            if let Some(f) = self.face_of(p) {
                                before.push(f);
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
                before.reverse();
                before.extend(faces);
                faces = before;
            }
            fans.push((faces, open));
        }
        if fans.len() > 1 {
            return VertexNeighbourhood::Fans { count: fans.len() };
        }
        let (faces, open) = fans.swap_remove(0);
        let pinch_uses = pinch_uses_of(&faces, open);
        if open {
            VertexNeighbourhood::HalfDisc { pinch_uses }
        } else {
            VertexNeighbourhood::Disc { pinch_uses }
        }
    }

    /// Every pinch use of the body, in vertex order: a face meeting a
    /// vertex more than once round it.
    pub fn pinch_uses(&self) -> Vec<(FaceHandle, VertexHandle)> {
        self.vertex_neighbourhoods()
            .into_iter()
            .flat_map(|(v, n)| match n {
                VertexNeighbourhood::Disc { pinch_uses }
                | VertexNeighbourhood::HalfDisc { pinch_uses } => {
                    pinch_uses.into_iter().map(|f| (f, v)).collect()
                }
                _ => Vec::new(),
            })
            .collect()
    }

    /// Decline `operation` at the body's first pinch use, or the first
    /// at one of `at` when given: for operations that do not handle a
    /// face touching itself at a vertex yet.
    ///
    /// # Errors
    ///
    /// [`CoreError::PinchUse`] naming the face and the vertex.
    pub fn require_no_pinch_use(
        &self,
        operation: &'static str,
        at: Option<&[VertexHandle]>,
    ) -> Result<(), crate::error::CoreError> {
        let first = self
            .pinch_uses()
            .into_iter()
            .find(|(_, v)| at.is_none_or(|at| at.contains(v)));
        match first {
            Some((face, vertex)) => Err(crate::error::CoreError::PinchUse {
                operation,
                face,
                vertex,
            }),
            None => Ok(()),
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

    /// The previous one: across the edge of the half-edge arriving at
    /// `v` before `h` in its loop.
    fn fan_prev(&self, h: HalfEdgeHandle, v: VertexHandle) -> Option<HalfEdgeHandle> {
        let p = self.half_edges.get(h)?.prev;
        let m = self.mate(p)?;
        (self.half_edges.get(m)?.origin == v).then_some(m)
    }
}

#[cfg(test)]
#[path = "neighbourhood_tests.rs"]
mod tests;
