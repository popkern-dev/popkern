//! A face's loops joined where they share a vertex, and parted again:
//! the Euler operators of a pinch use (PK-1129).
//!
//! A face can touch itself at a vertex: the band of a cylinder between
//! the two lobes of Viviani's figure eight, which touch at its double
//! point, or a plate's top face where a dimple's rim touches its edge.
//! Its interior is still a disc, so its boundary is ONE loop passing
//! through the vertex twice. Carried as two loops sharing the vertex
//! (an outer and an inner), the face counts as a disc with a hole and
//! the body's Euler–Poincaré sum is off by one (measured: a Viviani
//! common part read 4 for 2).
//!
//! [`BRep::join_loops_at_vertex`] splices two loops of one face at a
//! vertex both pass through into one loop visiting it twice: at the
//! vertex, each loop's arriving fin is relinked to the other loop's
//! leaving fin. Counts: L − 1 (one ring fewer on the face); V, E, F
//! unchanged; V − E + F − (L − F) rises by one, as the face's hole is
//! gone. [`BRep::split_loop_at_vertex`] is its inverse: a loop through a
//! vertex twice parts into two there, the new one an inner loop.

use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, Loop, LoopHandle, VertexHandle};

use super::EulerError;

impl BRep {
    /// Join loops `keep` and `absorb` of `face` at `vertex`, which each
    /// passes through once: `keep` becomes one loop through it twice and
    /// `absorb` is removed (the face's outer loop is kept as outer).
    ///
    /// # Errors
    ///
    /// [`EulerError::NotJoinable`] before anything changes when the
    /// loops are not two distinct boundary loops of `face` or do not
    /// each pass through `vertex` exactly once.
    pub fn join_loops_at_vertex(
        &mut self,
        face: FaceHandle,
        keep: LoopHandle,
        absorb: LoopHandle,
        vertex: VertexHandle,
    ) -> Result<(), EulerError> {
        let decline = |reason: &str| EulerError::NotJoinable {
            face,
            reason: reason.into(),
        };
        let f = self.faces.get(face).ok_or(decline("no such face"))?;
        let on_face = |l: LoopHandle| f.outer_loop == l || f.inner_loops.contains(&l);
        if keep == absorb || !on_face(keep) || !on_face(absorb) {
            return Err(decline("the loops are not two distinct loops of the face"));
        }
        let (Some(ha), Some(hb)) = (
            self.leaving_once(keep, vertex),
            self.leaving_once(absorb, vertex),
        ) else {
            return Err(decline(
                "a loop does not pass through the vertex exactly once",
            ));
        };
        let (pa, pb) = (self.half_edges[ha].prev, self.half_edges[hb].prev);
        // Every fin of `absorb` joins `keep`.
        for h in self.loop_fins(absorb) {
            self.half_edges[h].loop_ = keep;
        }
        // At the vertex: each arrival continues into the other's
        // departure.
        self.half_edges[pa].next = hb;
        self.half_edges[hb].prev = pa;
        self.half_edges[pb].next = ha;
        self.half_edges[ha].prev = pb;
        let outer_was_absorbed = self.faces[face].outer_loop == absorb;
        let f = &mut self.faces[face];
        f.inner_loops.retain(|&l| l != absorb && l != keep);
        if outer_was_absorbed || f.outer_loop == keep {
            f.outer_loop = keep;
        } else {
            f.inner_loops.push(keep);
        }
        self.loops.remove(absorb);
        Ok(())
    }

    /// Part loop `l` of `face`, which passes through `vertex` twice, into
    /// two loops there; returns the new one, added as an inner loop.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotJoinable`] before anything changes when `l` is
    /// not a loop of `face` or does not pass through `vertex` twice.
    pub fn split_loop_at_vertex(
        &mut self,
        face: FaceHandle,
        l: LoopHandle,
        vertex: VertexHandle,
    ) -> Result<LoopHandle, EulerError> {
        let decline = |reason: &str| EulerError::NotJoinable {
            face,
            reason: reason.into(),
        };
        let f = self.faces.get(face).ok_or(decline("no such face"))?;
        if f.outer_loop != l && !f.inner_loops.contains(&l) {
            return Err(decline("the loop is not the face's"));
        }
        let leaving: Vec<HalfEdgeHandle> = self
            .loop_fins(l)
            .into_iter()
            .filter(|&h| self.half_edges[h].origin == vertex)
            .collect();
        let [ha, hb] = leaving[..] else {
            return Err(decline("the loop does not pass through the vertex twice"));
        };
        let (pa, pb) = (self.half_edges[ha].prev, self.half_edges[hb].prev);
        // Undo the join: each arrival continues into its own run again.
        self.half_edges[pa].next = hb;
        self.half_edges[hb].prev = pa;
        self.half_edges[pb].next = ha;
        self.half_edges[ha].prev = pb;
        // The run from `ha` is the one `l` keeps if it holds `l`'s
        // start; the other run becomes the new loop.
        let start = self.loops[l].half_edge;
        let run_a = self.loop_fins_from(ha);
        let (kept_start, moved) = if run_a.contains(&start) {
            (ha, self.loop_fins_from(hb))
        } else {
            (hb, run_a)
        };
        let new = self.loops.insert(Loop::new(moved[0], face));
        for h in moved {
            self.half_edges[h].loop_ = new;
        }
        self.loops[l].half_edge = kept_start;
        self.faces[face].inner_loops.push(new);
        Ok(new)
    }

    /// Join, on every face, any two loops passing once each through a
    /// vertex they share, until none do; returns the joins made. A face
    /// whose loops touch at a vertex is one disc there, its boundary one
    /// loop through the vertex twice (see the module docs).
    pub fn join_loops_sharing_vertices(&mut self) -> usize {
        let mut joined = 0;
        let faces: Vec<FaceHandle> = self.faces.keys().collect();
        for face in faces {
            while let Some((keep, absorb, v)) = self.first_shared_vertex(face) {
                if self.join_loops_at_vertex(face, keep, absorb, v).is_err() {
                    break;
                }
                joined += 1;
            }
        }
        joined
    }

    /// Two loops of `face` and a vertex each passes through once: the
    /// first pair in loop order, and the first such vertex in the first
    /// loop's order. Each loop is walked once and its vertices indexed,
    /// so a face of many holes costs its fins, not their pairs (a disc
    /// top of 64 loops took 0.1 ms a boolean, PK-1161).
    fn first_shared_vertex(
        &self,
        face: FaceHandle,
    ) -> Option<(LoopHandle, LoopHandle, VertexHandle)> {
        let f = self.faces.get(face)?;
        if f.inner_loops.is_empty() {
            return None;
        }
        let loops: Vec<LoopHandle> = std::iter::once(f.outer_loop)
            .chain(f.inner_loops.iter().copied())
            .collect();
        let origins: Vec<Vec<VertexHandle>> = loops
            .iter()
            .map(|&l| {
                self.loop_fins(l)
                    .into_iter()
                    .map(|h| self.half_edges[h].origin)
                    .collect()
            })
            .collect();
        // Each vertex's passes per loop, as (loop index, passes).
        let mut passes: std::collections::HashMap<VertexHandle, Vec<(usize, usize)>> =
            std::collections::HashMap::new();
        for (i, vs) in origins.iter().enumerate() {
            for &v in vs {
                let entry = passes.entry(v).or_default();
                match entry.last_mut() {
                    Some((k, n)) if *k == i => *n += 1,
                    _ => entry.push((i, 1)),
                }
            }
        }
        let once_in = |v: &VertexHandle, i: usize| passes[v].iter().any(|&(k, n)| k == i && n == 1);
        for (i, va) in origins.iter().enumerate() {
            // The nearest later loop through a vertex of this one once
            // each, then that vertex first in this loop's order.
            let found = va
                .iter()
                .enumerate()
                .filter(|(_, v)| once_in(v, i))
                .flat_map(|(pos, v)| {
                    passes[v]
                        .iter()
                        .filter(move |&&(k, n)| k > i && n == 1)
                        .map(move |&(k, _)| (k, pos, *v))
                })
                .min_by_key(|&(k, pos, _)| (k, pos));
            if let Some((k, _, v)) = found {
                return Some((loops[i], loops[k], v));
            }
        }
        None
    }

    /// The fin of loop `l` leaving `vertex`, when exactly one does.
    fn leaving_once(&self, l: LoopHandle, vertex: VertexHandle) -> Option<HalfEdgeHandle> {
        let leaving: Vec<HalfEdgeHandle> = self
            .loop_fins(l)
            .into_iter()
            .filter(|&h| self.half_edges[h].origin == vertex)
            .collect();
        match leaving[..] {
            [h] => Some(h),
            _ => None,
        }
    }

    /// Loop `l`'s fins in order from its start.
    fn loop_fins(&self, l: LoopHandle) -> Vec<HalfEdgeHandle> {
        self.loops
            .get(l)
            .map_or_else(Vec::new, |lp| self.loop_fins_from(lp.half_edge))
    }

    /// The fins of the cycle under `next` from `start`.
    fn loop_fins_from(&self, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
        let mut out = Vec::new();
        let mut h = start;
        for _ in 0..=self.half_edges.len() {
            if !self.half_edges.contains_key(h) {
                break;
            }
            out.push(h);
            h = self.half_edges[h].next;
            if h == start {
                break;
            }
        }
        out
    }
}

#[cfg(test)]
#[path = "pinch_loop_tests.rs"]
mod tests;
