//! Wire edges, acorns and vertex uses (PK-1102 S3c).
//!
//! A wire edge (an edge on no face, [`BRep::is_wire_edge`]) and an
//! isolated vertex (an acorn) lie in a region like a sheet does, so each
//! is held by a shell of that region: [`Shell::wire_edges`] and
//! [`Shell::acorn`]. A wire body is a shell of the outside holding its
//! edges, one shell per connected wire. Producers do not place wires:
//! the builder settles every wire edge no shell holds into the outside
//! when the body is built ([`BRep::settle_wire_edges`]), as it settles
//! open shells. An acorn is made only on purpose
//! ([`BRep::add_acorn_shell`]): a vertex no edge reaches is otherwise an
//! orphan, which validation reports.
//!
//! A vertex's uses ([`VertexUse`]) are derived, not stored: one fan per
//! set of face uses about it connected through loops and edge rings,
//! one use per wire edge ending there, one per acorn shell holding it,
//! one per loop of no fins standing at it (a vertex on a face's
//! interior).
//! A vertex with more than one use is non-manifold (two cones touching
//! at their apexes, a wire standing on a face).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use slotmap::Key;

use super::{
    BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, RegionHandle, Shell, ShellHandle, VertexHandle,
};

/// One use of a vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum VertexUse {
    /// A fan of face uses about the vertex, named by one of its fins
    /// leaving the vertex.
    Fan(HalfEdgeHandle),
    /// A wire edge ending at the vertex.
    Wire(EdgeHandle),
    /// The acorn shell holding the vertex.
    Isolated(ShellHandle),
    /// A loop of no fins standing at the vertex: the vertex touches its
    /// face's interior (PK-1116), or is the first Euler step's.
    Loop(LoopHandle),
}

impl BRep {
    /// The shell holding wire edge `edge`, if any shell does.
    pub fn wire_edge_shell(&self, edge: EdgeHandle) -> Option<ShellHandle> {
        self.shells
            .iter()
            .find(|(_, s)| s.wire_edges.contains(&edge))
            .map(|(h, _)| h)
    }

    /// The acorn shells and the vertex each holds, in arena order.
    pub fn acorns(&self) -> Vec<(ShellHandle, VertexHandle)> {
        self.shells
            .iter()
            .filter_map(|(h, s)| s.acorn.map(|v| (h, v)))
            .collect()
    }

    /// The loops of no fins and the vertex each stands at, in arena
    /// order.
    pub fn vertex_loops(&self) -> Vec<(LoopHandle, VertexHandle)> {
        self.loops
            .iter()
            .filter(|(_, l)| !self.half_edges.contains_key(l.half_edge))
            .filter_map(|(h, l)| l.vertex.map(|v| (h, v)))
            .collect()
    }

    /// `face`'s vertex loops: the loops of no fins naming it that are
    /// none of its boundary loops (a vertex on its interior, PK-1116).
    pub fn face_vertex_loops(&self, face: super::FaceHandle) -> Vec<LoopHandle> {
        let outer = self.faces.get(face).map(|f| f.outer_loop);
        self.loops
            .iter()
            .filter(|&(h, l)| {
                l.face == face
                    && Some(h) != outer
                    && l.vertex.is_some()
                    && !self.half_edges.contains_key(l.half_edge)
            })
            .map(|(h, _)| h)
            .collect()
    }

    /// Every vertex's uses (see the module docs), in vertex order. A
    /// vertex with no use (an orphan) maps to an empty list.
    pub fn vertex_uses_all(&self) -> BTreeMap<VertexHandle, Vec<VertexUse>> {
        let mut out: BTreeMap<VertexHandle, Vec<VertexUse>> =
            self.vertices.keys().map(|v| (v, Vec::new())).collect();
        // Face uses leaving each vertex.
        let mut leaving: HashMap<VertexHandle, Vec<HalfEdgeHandle>> = HashMap::new();
        for (h, fin) in &self.half_edges {
            if self.loops.contains_key(fin.loop_) {
                leaving.entry(fin.origin).or_default().push(h);
            }
        }
        for (v, fins) in leaving {
            for fan in self.fans_about(v, &fins) {
                out.entry(v).or_default().push(VertexUse::Fan(fan));
            }
        }
        for edge in self.wire_edges() {
            let Some((a, b)) = self.edge_vertices(edge) else {
                continue;
            };
            out.entry(a).or_default().push(VertexUse::Wire(edge));
            if b != a {
                out.entry(b).or_default().push(VertexUse::Wire(edge));
            }
        }
        for (shell, v) in self.acorns() {
            out.entry(v).or_default().push(VertexUse::Isolated(shell));
        }
        for (l, v) in self.vertex_loops() {
            out.entry(v).or_default().push(VertexUse::Loop(l));
        }
        for uses in out.values_mut() {
            uses.sort_unstable();
        }
        out
    }

    /// `vertex`'s uses (see the module docs). Scans the body; ask
    /// [`BRep::vertex_uses_all`] for every vertex at once.
    pub fn vertex_uses(&self, vertex: VertexHandle) -> Vec<VertexUse> {
        let fins = self.face_uses_leaving(vertex);
        let mut uses: Vec<VertexUse> = self
            .fans_about(vertex, &fins)
            .into_iter()
            .map(VertexUse::Fan)
            .collect();
        uses.extend(
            self.wire_edges()
                .into_iter()
                .filter(|&e| {
                    self.edge_vertices(e)
                        .is_some_and(|(a, b)| a == vertex || b == vertex)
                })
                .map(VertexUse::Wire),
        );
        uses.extend(
            self.acorns()
                .into_iter()
                .filter(|&(_, v)| v == vertex)
                .map(|(s, _)| VertexUse::Isolated(s)),
        );
        uses.extend(
            self.vertex_loops()
                .into_iter()
                .filter(|&(_, v)| v == vertex)
                .map(|(l, _)| VertexUse::Loop(l)),
        );
        uses.sort_unstable();
        uses
    }

    /// Whether `vertex` has more than one use: the body is non-manifold
    /// there.
    pub fn is_nonmanifold_vertex(&self, vertex: VertexHandle) -> bool {
        self.vertex_uses(vertex).len() > 1
    }

    /// The fans the face uses `fins` leaving `vertex` fall into, each
    /// named by its least fin.
    fn fans_about(&self, vertex: VertexHandle, fins: &[HalfEdgeHandle]) -> Vec<HalfEdgeHandle> {
        self.fan_groups(vertex, fins)
            .into_iter()
            .filter_map(|g| g.first().copied())
            .collect()
    }

    /// The face uses (fins in a loop) leaving `vertex`, in arena order.
    pub(crate) fn face_uses_leaving(&self, vertex: VertexHandle) -> Vec<HalfEdgeHandle> {
        self.half_edges
            .iter()
            .filter(|(_, f)| f.origin == vertex && self.loops.contains_key(f.loop_))
            .map(|(h, _)| h)
            .collect()
    }

    /// The face use about `vertex` ring fin `r` stands for at the end of
    /// its edge where fin `reference` leaves the vertex (`leaves`) or
    /// arrives at it: `r` itself when it leaves there, its loop successor
    /// when it arrives. On an open edge the end is the vertex's; on a
    /// closed edge both ends are, and a fin leaves at the reference's end
    /// when it runs the edge the reference's way (its mate in a ring of
    /// two arrives there). Read by origin alone, every fin of a closed
    /// edge left at both ends, and the two lobes of a figure eight
    /// meeting at its double point read as one fan though the
    /// neighbourhood walk ([`BRep::vertex_neighbourhood`]) finds two.
    pub(crate) fn ring_use_at(
        &self,
        vertex: VertexHandle,
        r: HalfEdgeHandle,
        reference: HalfEdgeHandle,
        leaves: bool,
    ) -> Option<HalfEdgeHandle> {
        let fin = self.half_edges.get(r)?;
        // Closed: the fin ends where it starts (its successor leaves its
        // origin) round a curve of length, its middle off the vertex. A
        // degenerate edge (a pole's) leaves and arrives through no arc,
        // read by origin: a pole's fan of collapsed fins is one fan.
        let closed = self.half_edges.get(fin.next).map(|n| n.origin) == Some(fin.origin)
            && self.closed_round_curve(fin.edge, fin.origin);
        let leaves_here = if closed {
            (self.fin_sense(r)? == self.fin_sense(reference)?) == leaves
        } else {
            fin.origin == vertex
        };
        Some(if leaves_here { r } else { fin.next })
    }

    /// Whether `edge`, whose ends are both `vertex`, runs round a curve
    /// whose middle lies farther than the linear tolerance from it (a
    /// circle's, not a pole's collapsed edge).
    fn closed_round_curve(&self, edge: EdgeHandle, vertex: VertexHandle) -> bool {
        use crate::geometry::curve::CurveEval;
        let (Some(c), Some(v)) = (
            self.edges.get(edge).and_then(|e| e.curve.as_ref()),
            self.vertices.get(vertex),
        ) else {
            return false;
        };
        c.geometry
            .point_at(0.5 * (c.t_start + c.t_end))
            .distance_to(&v.position)
            > self.tolerance.linear()
    }

    /// The face uses `fins` leaving `vertex`, grouped into fans: two
    /// fins are in one fan when one follows the other's loop predecessor
    /// across an edge ring, or they share an edge's ring. Each group is
    /// sorted, the groups by their least fin.
    pub(crate) fn fan_groups(
        &self,
        vertex: VertexHandle,
        fins: &[HalfEdgeHandle],
    ) -> Vec<Vec<HalfEdgeHandle>> {
        let index: HashMap<HalfEdgeHandle, usize> =
            fins.iter().enumerate().map(|(i, &h)| (h, i)).collect();
        let mut parent: Vec<usize> = (0..fins.len()).collect();
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let leaving_use = |r: HalfEdgeHandle, reference: HalfEdgeHandle, leaves: bool| {
            index
                .get(&self.ring_use_at(vertex, r, reference, leaves)?)
                .copied()
        };
        for (i, &h) in fins.iter().enumerate() {
            let Some(fin) = self.half_edges.get(h) else {
                continue;
            };
            let along = [(fin.edge, h, true)].into_iter().chain(
                self.half_edges
                    .get(fin.prev)
                    .map(|p| (p.edge, fin.prev, false)),
            );
            for (edge, reference, leaves) in along {
                for r in self.radial_fins(edge) {
                    if let Some(j) = leaving_use(r, reference, leaves) {
                        let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                        parent[a] = b;
                    }
                }
            }
        }
        let mut fans: BTreeMap<usize, Vec<HalfEdgeHandle>> = BTreeMap::new();
        for (i, &h) in fins.iter().enumerate() {
            let r = root(&mut parent, i);
            fans.entry(r).or_default().push(h);
        }
        let mut out: Vec<Vec<HalfEdgeHandle>> = fans
            .into_values()
            .map(|mut g| {
                g.sort_unstable();
                g
            })
            .collect();
        out.sort_unstable();
        out
    }

    /// Hold every wire edge no shell holds in a new shell of the
    /// outside, one shell per set of such edges connected through their
    /// vertices, and drop held entries that are no longer wire edges (an
    /// edge removed or given a face); a shell left holding nothing that
    /// held wire edges goes. Called where a producer has finished a body
    /// that may carry wires (the builder, old bodies loading). Returns
    /// the shells made.
    pub fn settle_wire_edges(&mut self) -> Vec<ShellHandle> {
        let wires: BTreeSet<EdgeHandle> = self.wire_edges().into_iter().collect();
        let mut emptied = Vec::new();
        let mut held = BTreeSet::new();
        for (h, shell) in self.shells.iter_mut() {
            if shell.wire_edges.is_empty() {
                continue;
            }
            shell
                .wire_edges
                .retain(|e| wires.contains(e) && held.insert(*e));
            if shell.wire_edges.is_empty() {
                emptied.push(h);
            }
        }
        for h in emptied {
            let faceless = !self.faces.values().any(|f| f.back == h || f.front == h);
            if faceless && self.shells[h].acorn.is_none() {
                self.shells.remove(h);
            }
        }
        let loose: Vec<EdgeHandle> = wires.difference(&held).copied().collect();
        let mut made = Vec::new();
        for component in self.wire_components(&loose) {
            let mut shell = Shell::new(Default::default(), self.outside);
            shell.wire_edges = component;
            made.push(self.shells.insert(shell));
        }
        made
    }

    /// `edges` split into sets connected through shared vertices, each
    /// in edge order.
    pub(crate) fn wire_components(&self, edges: &[EdgeHandle]) -> Vec<Vec<EdgeHandle>> {
        let mut parent: Vec<usize> = (0..edges.len()).collect();
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let mut at: HashMap<VertexHandle, usize> = HashMap::new();
        for (i, &e) in edges.iter().enumerate() {
            let Some((a, b)) = self.edge_vertices(e) else {
                continue;
            };
            for v in [a, b] {
                match at.get(&v) {
                    Some(&j) => {
                        let (x, y) = (root(&mut parent, i), root(&mut parent, j));
                        parent[x] = y;
                    }
                    None => {
                        at.insert(v, i);
                    }
                }
            }
        }
        let mut groups: BTreeMap<usize, Vec<EdgeHandle>> = BTreeMap::new();
        for (i, &e) in edges.iter().enumerate() {
            let r = root(&mut parent, i);
            groups.entry(r).or_default().push(e);
        }
        let mut out: Vec<Vec<EdgeHandle>> = groups.into_values().collect();
        out.sort_unstable();
        out
    }

    /// Hold `vertex` as an isolated vertex (an acorn) in a new shell of
    /// `region`. Returns the shell, or `None` when the region or vertex
    /// is unknown, or an edge or another shell already uses the vertex.
    pub fn add_acorn_shell(
        &mut self,
        region: RegionHandle,
        vertex: VertexHandle,
    ) -> Option<ShellHandle> {
        if !self.regions.contains_key(region) || !self.vertices.contains_key(vertex) {
            return None;
        }
        let used = self.half_edges.values().any(|f| f.origin == vertex)
            || self.shells.values().any(|s| s.acorn == Some(vertex));
        if used {
            return None;
        }
        let mut shell = Shell::new(Default::default(), region);
        shell.acorn = Some(vertex);
        let shell = self.shells.insert(shell);
        self.vertices[vertex].outgoing = HalfEdgeHandle::null();
        Some(shell)
    }
}

#[cfg(test)]
#[path = "brep_wires_tests.rs"]
mod tests;
