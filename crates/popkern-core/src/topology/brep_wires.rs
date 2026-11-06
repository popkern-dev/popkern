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
//! one use per wire edge ending there, one per acorn shell holding it.
//! A vertex with more than one use is non-manifold (two cones touching
//! at their apexes, a wire standing on a face).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use slotmap::Key;

use super::{BRep, EdgeHandle, HalfEdgeHandle, RegionHandle, Shell, ShellHandle, VertexHandle};

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
        for uses in out.values_mut() {
            uses.sort_unstable();
        }
        out
    }

    /// `vertex`'s uses (see the module docs). Scans the body; ask
    /// [`BRep::vertex_uses_all`] for every vertex at once.
    pub fn vertex_uses(&self, vertex: VertexHandle) -> Vec<VertexUse> {
        let fins: Vec<HalfEdgeHandle> = self
            .half_edges
            .iter()
            .filter(|(_, f)| f.origin == vertex && self.loops.contains_key(f.loop_))
            .map(|(h, _)| h)
            .collect();
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
        uses.sort_unstable();
        uses
    }

    /// Whether `vertex` has more than one use: the body is non-manifold
    /// there.
    pub fn is_nonmanifold_vertex(&self, vertex: VertexHandle) -> bool {
        self.vertex_uses(vertex).len() > 1
    }

    /// The fans the face uses `fins` leaving `vertex` fall into, each
    /// named by its least fin. Two fins are in one fan when one follows
    /// the other's loop predecessor across an edge ring, or they share
    /// an edge's ring.
    fn fans_about(&self, vertex: VertexHandle, fins: &[HalfEdgeHandle]) -> Vec<HalfEdgeHandle> {
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
        // The face use about `vertex` a ring fin `r` stands for: itself
        // when it leaves the vertex, its loop successor when it arrives.
        let leaving_use = |r: HalfEdgeHandle| -> Option<usize> {
            let fin = self.half_edges.get(r)?;
            if fin.origin == vertex {
                index.get(&r).copied()
            } else {
                index.get(&fin.next).copied()
            }
        };
        for (i, &h) in fins.iter().enumerate() {
            let Some(fin) = self.half_edges.get(h) else {
                continue;
            };
            let prev_edge = self.half_edges.get(fin.prev).map(|p| p.edge);
            for edge in [Some(fin.edge), prev_edge].into_iter().flatten() {
                for r in self.radial_fins(edge) {
                    if let Some(j) = leaving_use(r) {
                        let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                        parent[a] = b;
                    }
                }
            }
        }
        let mut fans: BTreeMap<usize, HalfEdgeHandle> = BTreeMap::new();
        for (i, &h) in fins.iter().enumerate() {
            let r = root(&mut parent, i);
            let least = fans.entry(r).or_insert(h);
            if h < *least {
                *least = h;
            }
        }
        let mut out: Vec<HalfEdgeHandle> = fans.into_values().collect();
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
    fn wire_components(&self, edges: &[EdgeHandle]) -> Vec<Vec<EdgeHandle>> {
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
