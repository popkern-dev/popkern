//! make_face_from_wire_cycle / kill_face_to_wire_cycle: a lamina face
//! spanning a closed cycle of wire edges, and back.
//!
//! The face lies in the region its wires did, both sides in their shell
//! (a sheet piece in the void, a lamina inside a solid). Each edge
//! keeps the fin the cycle runs along, now the face's use and alone on
//! its ring (a free edge); its mate goes. Geometry is the caller's: the
//! face takes the surface given, and pcurves are authored afterwards
//! (`BRep::recompute_face_pcurves`), as for [`BRep::closed_lamina`].
//!
//! [`BRep::closed_lamina`]: crate::topology::BRep::closed_lamina

use std::collections::BTreeSet;

use slotmap::Key;

use super::EulerError;
use crate::geometry::surface::SurfaceOnFace;
use crate::topology::{
    BRep, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, RadialLinks,
};

/// What [`BRep::make_face_from_wire_cycle`] made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WireFaceResult {
    /// The face.
    pub face: FaceHandle,
    /// Its loop, running the cycle's fins in order.
    pub loop_: LoopHandle,
}

impl BRep {
    /// Make a face bounded by the closed cycle of wire-edge fins
    /// `cycle` (each leaving the previous one's end, the last returning
    /// to the first's origin), carrying `surface`. Counts: faces +1,
    /// loops +1, fins −n, wire edges −n.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotAWireCycle`] when the fins are not a closed
    /// cycle of distinct wire edges held by one shell.
    pub fn make_face_from_wire_cycle(
        &mut self,
        cycle: &[HalfEdgeHandle],
        surface: Option<SurfaceOnFace>,
    ) -> Result<WireFaceResult, EulerError> {
        let shell = self.check_wire_cycle(cycle)?;
        let face = self.faces.insert(Face {
            outer_loop: LoopHandle::default(),
            inner_loops: Vec::new(),
            back: shell,
            front: shell,
            surface,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let loop_ = self.loops.insert(Loop {
            half_edge: cycle[0],
            face,
        });
        self.faces[face].outer_loop = loop_;
        let n = cycle.len();
        for (i, &fin) in cycle.iter().enumerate() {
            if let Some(mate) = self.mate(fin) {
                self.remove_half_edge_reanchored(mate);
            }
            let edge = self.half_edges[fin].edge;
            let h = &mut self.half_edges[fin];
            h.next = cycle[(i + 1) % n];
            h.prev = cycle[(i + n - 1) % n];
            h.loop_ = loop_;
            h.ring = RadialLinks::alone();
            self.edges[edge].half_edge = fin;
        }
        let edges: BTreeSet<EdgeHandle> = cycle.iter().map(|&f| self.half_edges[f].edge).collect();
        let s = &mut self.shells[shell];
        s.wire_edges.retain(|e| !edges.contains(e));
        if !self.faces.contains_key(self.shells[shell].face) {
            self.shells[shell].face = face;
        }
        Ok(WireFaceResult { face, loop_ })
    }

    /// Turn lamina `face` back into the wire cycle it spans, the inverse
    /// of [`BRep::make_face_from_wire_cycle`]: every edge of its one
    /// loop must be free (the face its only use) and both its sides in
    /// one shell. Each edge regains a mate and is held as a wire by that
    /// shell; the face, its loop and the fins' pcurves go. Returns the
    /// edges in loop order. Counts: faces −1, loops −1, fins +n, wire
    /// edges +n.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotALaminaFace`] otherwise.
    pub fn kill_face_to_wire_cycle(
        &mut self,
        face: FaceHandle,
    ) -> Result<Vec<EdgeHandle>, EulerError> {
        let data = self
            .faces
            .get(face)
            .ok_or(EulerError::NotALaminaFace { face })?;
        let (shell, loop_) = (data.back, data.outer_loop);
        if data.front != shell || !data.inner_loops.is_empty() || !self.shells.contains_key(shell) {
            return Err(EulerError::NotALaminaFace { face });
        }
        let fins = self.loop_cycle(loop_);
        let free = !fins.is_empty()
            && fins.iter().all(|&f| {
                self.half_edges[f].ring.next.is_null() && self.half_edges[f].ring.prev.is_null()
            });
        if !free {
            return Err(EulerError::NotALaminaFace { face });
        }
        let ends: Vec<_> = fins
            .iter()
            .map(|&f| self.half_edges[self.half_edges[f].next].origin)
            .collect();
        let mut edges = Vec::with_capacity(fins.len());
        for (&fin, end) in fins.iter().zip(ends) {
            let edge = self.half_edges[fin].edge;
            let mate = self.half_edges.insert(HalfEdge {
                origin: end,
                ring: RadialLinks::alone(),
                next: HalfEdgeHandle::default(),
                prev: HalfEdgeHandle::default(),
                edge,
                loop_: LoopHandle::default(),
                pcurve: None,
            });
            let h = &mut self.half_edges[fin];
            h.next = HalfEdgeHandle::default();
            h.prev = HalfEdgeHandle::default();
            h.loop_ = LoopHandle::default();
            h.pcurve = None;
            self.set_mates(fin, mate);
            self.face_use_pcurves
                .remove(&crate::topology::FaceUseTrimKey {
                    face,
                    half_edge: fin,
                });
            edges.push(edge);
        }
        self.loops.remove(loop_);
        self.faces.remove(face);
        self.shells[shell].wire_edges.extend(edges.iter().copied());
        if self.shells[shell].face == face {
            let next = self
                .faces
                .iter()
                .find(|(_, f)| f.back == shell || f.front == shell)
                .map(|(h, _)| h)
                .unwrap_or_default();
            self.shells[shell].face = next;
        }
        Ok(edges)
    }

    /// The shell holding every edge of `cycle`, when the fins are a
    /// closed cycle of distinct held wire edges.
    fn check_wire_cycle(
        &self,
        cycle: &[HalfEdgeHandle],
    ) -> Result<crate::topology::ShellHandle, EulerError> {
        let bad = |reason: &'static str| EulerError::NotAWireCycle { reason };
        let mut shell = None;
        let mut edges = BTreeSet::new();
        for (i, &fin) in cycle.iter().enumerate() {
            let h = self
                .half_edges
                .get(fin)
                .ok_or(bad("a fin is not in the body"))?;
            if !self.is_wire_edge(h.edge) || !edges.insert(h.edge) {
                return Err(bad("a fin's edge is not a distinct wire edge"));
            }
            let held = self
                .wire_edge_shell(h.edge)
                .ok_or(bad("a wire edge is held by no shell"))?;
            if *shell.get_or_insert(held) != held {
                return Err(bad("the edges are held by different shells"));
            }
            let next = cycle[(i + 1) % cycle.len()];
            let end = self.fin_end(fin).ok_or(bad("a fin has no end"))?;
            if self.half_edges.get(next).map(|n| n.origin) != Some(end) {
                return Err(bad("a fin does not leave the previous fin's end"));
            }
        }
        shell.ok_or(bad("the cycle is empty"))
    }

    /// The fins of loop `l` in order.
    pub(super) fn loop_cycle(&self, l: LoopHandle) -> Vec<HalfEdgeHandle> {
        let Some(start) = self.loops.get(l).map(|l| l.half_edge) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let mut cursor = start;
        for _ in 0..self.half_edges.len() {
            if self.half_edges.get(cursor).is_none() {
                return Vec::new();
            }
            out.push(cursor);
            cursor = self.half_edges[cursor].next;
            if cursor == start {
                return out;
            }
        }
        Vec::new()
    }
}
