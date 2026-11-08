//! make_wire_edge / kill_wire_edge: an edge on no face, held by a shell
//! of the region it lies in.
//!
//! A shell is one connected piece of a region's boundary, so a wire
//! joins the shells its ends lie in: an acorn end is absorbed (the
//! vertex is now the wire's), a wire-only shell at the other end merges
//! into the first, and a wire standing on a face joins that face's
//! shell. Two shells holding face sides are not joined by a wire (their
//! faces stay separate surfaces, `EulerError::WireJoinsFaceShells`).
//! Killing a wire edge splits its shell where that disconnects it and
//! leaves an end no other edge reaches as an acorn, so a vertex is never
//! orphaned.

use std::collections::{BTreeSet, HashSet};

use super::EulerError;
use crate::topology::{
    BRep, Edge, EdgeHandle, HalfEdge, HalfEdgeHandle, LoopHandle, RadialLinks, RegionHandle, Shell,
    ShellHandle, VertexHandle,
};

/// What [`BRep::make_wire_edge`] made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WireEdgeResult {
    /// The edge, its first fin running `v0` → `v1`.
    pub edge: EdgeHandle,
    /// The fin leaving `v0`.
    pub fin: HalfEdgeHandle,
    /// The fin leaving `v1`.
    pub mate: HalfEdgeHandle,
    /// The shell holding the edge.
    pub shell: ShellHandle,
}

/// What [`BRep::kill_wire_edge`] changed besides removing the edge.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct WireEdgeKill {
    /// Ends no edge reaches any more, each now an acorn in its shell.
    pub acorns: Vec<(VertexHandle, ShellHandle)>,
    /// Shells split off where removing the edge disconnected its shell.
    pub split: Vec<ShellHandle>,
}

/// Where a vertex lies in a region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Holder {
    /// An acorn shell of the region holds it.
    Acorn(ShellHandle),
    /// A shell of the region reaches it through wire edges or faces.
    Shell(ShellHandle),
    /// Nothing uses the vertex.
    Orphan,
}

impl Holder {
    fn shell(self) -> Option<ShellHandle> {
        match self {
            Holder::Acorn(s) | Holder::Shell(s) => Some(s),
            Holder::Orphan => None,
        }
    }
}

impl BRep {
    /// Make a wire edge from `v0` to `v1` in `region`. Counts: edges
    /// +1, fins +2, wire edges +1; an acorn end absorbed: acorns −1
    /// (and shells −1 unless its shell holds the edge); a second
    /// wire-only shell merged: shells −1. The edge may close a cycle.
    ///
    /// # Errors
    ///
    /// Unknown handles, [`EulerError::DegenerateWireEdge`] for
    /// `v0 == v1`, [`EulerError::VertexNotInRegion`] for an end whose
    /// uses lie in other regions, [`EulerError::WireJoinsFaceShells`].
    pub fn make_wire_edge(
        &mut self,
        region: RegionHandle,
        v0: VertexHandle,
        v1: VertexHandle,
    ) -> Result<WireEdgeResult, EulerError> {
        if !self.regions.contains_key(region) {
            return Err(EulerError::UnknownRegion(region));
        }
        for v in [v0, v1] {
            if !self.vertices.contains_key(v) {
                return Err(EulerError::UnknownVertex(v));
            }
        }
        if v0 == v1 {
            return Err(EulerError::DegenerateWireEdge { vertex: v0 });
        }
        let (h0, h1) = (self.holder(v0, region)?, self.holder(v1, region)?);
        let has_faces = |s: Option<ShellHandle>| s.is_some_and(|s| self.shell_has_faces(s));
        if let (Some(a), Some(b)) = (h0.shell(), h1.shell())
            && a != b
            && has_faces(Some(a))
            && has_faces(Some(b))
        {
            return Err(EulerError::WireJoinsFaceShells(a, b));
        }
        // The edge's shell: a face shell when an end has one, else
        // either end's, else a new one.
        let (keep, other) = if has_faces(h1.shell()) {
            (h1, h0)
        } else {
            (h0, h1)
        };
        let target = match keep.shell().or(other.shell()) {
            Some(s) => s,
            None => self.shells.insert(Shell::new(Default::default(), region)),
        };
        for holder in [keep, other] {
            match holder {
                Holder::Acorn(s) if s == target => self.shells[s].acorn = None,
                Holder::Acorn(s) => {
                    self.shells.remove(s);
                }
                Holder::Shell(s) if s != target => {
                    let moved = std::mem::take(&mut self.shells[s].wire_edges);
                    self.shells[target].wire_edges.extend(moved);
                    self.shells.remove(s);
                }
                _ => {}
            }
        }
        let (edge, fin, mate) = self.insert_wire_fins(v0, v1);
        self.shells[target].wire_edges.push(edge);
        Ok(WireEdgeResult {
            edge,
            fin,
            mate,
            shell: target,
        })
    }

    /// Remove wire edge `edge`, the inverse of [`BRep::make_wire_edge`]:
    /// its shell splits where that disconnects it, and an end no edge
    /// reaches any more becomes an acorn (reusing the shell when the
    /// shell held only this edge). Counts: edges −1, fins −2, wire edges
    /// −1, plus one acorn per lone end and one shell per split piece or
    /// further acorn.
    ///
    /// # Errors
    ///
    /// [`EulerError::UnknownEdge`], [`EulerError::NotAWireEdge`].
    pub fn kill_wire_edge(&mut self, edge: EdgeHandle) -> Result<WireEdgeKill, EulerError> {
        if !self.edges.contains_key(edge) {
            return Err(EulerError::UnknownEdge(edge));
        }
        let shell = self
            .wire_edge_shell(edge)
            .filter(|_| self.is_wire_edge(edge))
            .ok_or(EulerError::NotAWireEdge { edge })?;
        let (a, b) = self
            .edge_vertices(edge)
            .ok_or(EulerError::NotAWireEdge { edge })?;
        let region = self.shells[shell].region;
        let fins: Vec<HalfEdgeHandle> = self
            .half_edges
            .iter()
            .filter(|(_, h)| h.edge == edge)
            .map(|(k, _)| k)
            .collect();
        for f in fins {
            self.remove_half_edge_reanchored(f);
        }
        self.edges.remove(edge);
        self.shells[shell].wire_edges.retain(|&e| e != edge);
        for v in [a, b] {
            self.reanchor_vertex(v);
        }
        let mut out = WireEdgeKill {
            split: self.split_wire_shell(shell),
            ..Default::default()
        };
        let mut ends = vec![a];
        if b != a {
            ends.push(b);
        }
        for v in ends {
            if self.half_edges.values().any(|h| h.origin == v) {
                continue;
            }
            let reuse = self.shells.get(shell).is_some_and(|s| {
                s.wire_edges.is_empty() && s.acorn.is_none() && !self.shell_has_faces(shell)
            });
            let held = if reuse {
                self.shells[shell].acorn = Some(v);
                shell
            } else {
                let mut s = Shell::new(Default::default(), region);
                s.acorn = Some(v);
                self.shells.insert(s)
            };
            out.acorns.push((v, held));
        }
        if self.shells.get(shell).is_some_and(|s| {
            s.wire_edges.is_empty() && s.acorn.is_none() && !self.shell_has_faces(shell)
        }) {
            self.shells.remove(shell);
        }
        Ok(out)
    }

    /// Where `vertex` lies in `region`.
    fn holder(&self, vertex: VertexHandle, region: RegionHandle) -> Result<Holder, EulerError> {
        if let Some(s) = self.acorn_shell(vertex) {
            return if self.shells[s].region == region {
                Ok(Holder::Acorn(s))
            } else {
                Err(EulerError::VertexNotInRegion { vertex, region })
            };
        }
        let mut used = false;
        let mut shells = BTreeSet::new();
        for fin in self.half_edges.values().filter(|h| h.origin == vertex) {
            used = true;
            if let Some(s) = self.wire_edge_shell(fin.edge) {
                shells.insert(s);
            } else if let Some(face) = self.loops.get(fin.loop_).map(|l| l.face)
                && let Some(f) = self.faces.get(face)
            {
                shells.insert(f.back);
                shells.insert(f.front);
            }
        }
        // A fin arriving at the vertex on a face (a free edge's lone
        // fin) also places it: its loop successor leaves the vertex, so
        // the scan above already saw that face.
        let in_region = shells
            .into_iter()
            .find(|&s| self.shells.get(s).is_some_and(|s| s.region == region));
        match (in_region, used) {
            (Some(s), _) => Ok(Holder::Shell(s)),
            (None, true) => Err(EulerError::VertexNotInRegion { vertex, region }),
            (None, false) => Ok(Holder::Orphan),
        }
    }

    /// Whether some face has a side in `shell`.
    pub(crate) fn shell_has_faces(&self, shell: ShellHandle) -> bool {
        self.faces
            .values()
            .any(|f| f.back == shell || f.front == shell)
    }

    /// A wire edge's two fins, mated, in no loop.
    fn insert_wire_fins(
        &mut self,
        v0: VertexHandle,
        v1: VertexHandle,
    ) -> (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle) {
        let edge = self.insert_edge_tracked(Edge {
            half_edge: HalfEdgeHandle::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        let fin = |origin| HalfEdge {
            origin,
            ring: RadialLinks::alone(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge,
            loop_: LoopHandle::default(),
            pcurve: None,
        };
        let a = self.half_edges.insert(fin(v0));
        let b = self.half_edges.insert(fin(v1));
        self.set_mates(a, b);
        self.edges[edge].half_edge = a;
        for (v, h) in [(v0, a), (v1, b)] {
            if !self.half_edges.contains_key(self.vertices[v].outgoing) {
                self.vertices[v].outgoing = h;
            }
        }
        (edge, a, b)
    }

    /// Point `vertex`'s outgoing fin at a live fin leaving it, or null.
    fn reanchor_vertex(&mut self, vertex: VertexHandle) {
        let Some(v) = self.vertices.get(vertex) else {
            return;
        };
        if self.half_edges.contains_key(v.outgoing) {
            return;
        }
        let next = self
            .half_edges
            .iter()
            .find(|(_, h)| h.origin == vertex)
            .map(|(k, _)| k)
            .unwrap_or_default();
        self.vertices[vertex].outgoing = next;
    }

    /// Split `shell`'s wire edges into the pieces still connected to it:
    /// with face sides, the wires reaching a vertex of its faces stay
    /// and every other connected set of wires gets a shell of its own;
    /// without, the first connected set stays. Returns the new shells.
    pub(super) fn split_wire_shell(&mut self, shell: ShellHandle) -> Vec<ShellHandle> {
        let Some(s) = self.shells.get(shell) else {
            return Vec::new();
        };
        let (region, wires) = (s.region, s.wire_edges.clone());
        if wires.is_empty() {
            return Vec::new();
        }
        let face_vertices: HashSet<VertexHandle> = self
            .faces
            .values()
            .filter(|f| f.back == shell || f.front == shell)
            .flat_map(|f| std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()))
            .flat_map(|l| self.loop_origins(l))
            .collect();
        let mut pieces = self.wire_components(&wires);
        let touches_faces = |piece: &[EdgeHandle], brep: &BRep| {
            piece.iter().any(|&e| {
                brep.edge_vertices(e)
                    .is_some_and(|(a, b)| face_vertices.contains(&a) || face_vertices.contains(&b))
            })
        };
        let stay = if face_vertices.is_empty() {
            Some(0)
        } else {
            pieces.iter().position(|p| touches_faces(p, self))
        };
        let mut kept = Vec::new();
        let mut made = Vec::new();
        for (i, piece) in pieces.drain(..).enumerate() {
            let keep =
                Some(i) == stay || (!face_vertices.is_empty() && touches_faces(&piece, self));
            if keep {
                kept.extend(piece);
            } else {
                let mut s = Shell::new(Default::default(), region);
                s.wire_edges = piece;
                made.push(self.shells.insert(s));
            }
        }
        self.shells[shell].wire_edges = kept;
        made
    }

    /// The origins of loop `l`'s fins.
    fn loop_origins(&self, l: LoopHandle) -> Vec<VertexHandle> {
        let Some(start) = self.loops.get(l).map(|l| l.half_edge) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let mut cursor = start;
        for _ in 0..self.half_edges.len() {
            let Some(h) = self.half_edges.get(cursor) else {
                break;
            };
            out.push(h.origin);
            cursor = h.next;
            if cursor == start {
                break;
            }
        }
        out
    }
}
