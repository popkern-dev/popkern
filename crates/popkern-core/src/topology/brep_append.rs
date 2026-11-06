//! Appending one body to another (PK-1102): every entity of the
//! appended body copied into fresh handles, its links remapped, no
//! welding. The two bodies' regions stay apart; they share nothing.

use std::collections::HashMap;

use slotmap::{Key, SlotMap};

use super::{
    BRep, ContactEntity, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, ShellHandle,
    SolidHandle, VertexHandle,
};

/// Where each entity of an appended body went ([`BRep::append_body`]).
#[derive(Debug, Clone, Default)]
pub struct AppendMap {
    /// Vertices.
    pub vertices: HashMap<VertexHandle, VertexHandle>,
    /// Edges.
    pub edges: HashMap<EdgeHandle, EdgeHandle>,
    /// Fins.
    pub half_edges: HashMap<HalfEdgeHandle, HalfEdgeHandle>,
    /// Loops.
    pub loops: HashMap<LoopHandle, LoopHandle>,
    /// Faces.
    pub faces: HashMap<FaceHandle, FaceHandle>,
    /// Shells.
    pub shells: HashMap<ShellHandle, ShellHandle>,
    /// Regions.
    pub regions: HashMap<SolidHandle, SolidHandle>,
}

/// `h` through `map`; a handle the appended body did not hold (a null
/// or dangling link of a body under construction) stays null, so a
/// broken input is carried as broken rather than pointing into `self`.
fn remap<K: Key>(map: &HashMap<K, K>, h: K) -> K {
    map.get(&h).copied().unwrap_or_default()
}

/// Insert a clone of every value of `from` into `into`, recording
/// where each went.
fn copy_arena<K: Key, V: Clone>(
    from: &SlotMap<K, V>,
    into: &mut SlotMap<K, V>,
    map: &mut HashMap<K, K>,
) {
    for (h, value) in from {
        map.insert(h, into.insert(value.clone()));
    }
}

impl BRep {
    /// Append `other` to this body: every vertex, edge, fin, loop,
    /// face, shell and region of `other` is copied into fresh handles
    /// with its geometry, names, tolerances and roles, every link
    /// (loop cycles, edge rings, anchors, ownership) remapped, and its
    /// recorded contacts carried. Nothing is welded: the result holds
    /// both bodies side by side. Returns where each entity went.
    pub fn append_body(&mut self, other: &BRep) -> AppendMap {
        let mut m = AppendMap::default();
        copy_arena(&other.vertices, &mut self.vertices, &mut m.vertices);
        for (h, edge) in &other.edges {
            let new = self.insert_edge_tracked(edge.clone());
            m.edges.insert(h, new);
        }
        copy_arena(&other.half_edges, &mut self.half_edges, &mut m.half_edges);
        copy_arena(&other.loops, &mut self.loops, &mut m.loops);
        copy_arena(&other.faces, &mut self.faces, &mut m.faces);
        copy_arena(&other.shells, &mut self.shells, &mut m.shells);
        copy_arena(&other.solids, &mut self.solids, &mut m.regions);

        for &new in m.vertices.values() {
            let v = &mut self.vertices[new];
            v.outgoing = remap(&m.half_edges, v.outgoing);
        }
        for &new in m.edges.values() {
            let e = &mut self.edges[new];
            e.half_edge = remap(&m.half_edges, e.half_edge);
        }
        for &new in m.half_edges.values() {
            let f = &mut self.half_edges[new];
            f.origin = remap(&m.vertices, f.origin);
            f.ring.next = remap(&m.half_edges, f.ring.next);
            f.ring.prev = remap(&m.half_edges, f.ring.prev);
            f.next = remap(&m.half_edges, f.next);
            f.prev = remap(&m.half_edges, f.prev);
            f.edge = remap(&m.edges, f.edge);
            f.loop_ = remap(&m.loops, f.loop_);
        }
        for &new in m.loops.values() {
            let l = &mut self.loops[new];
            l.half_edge = remap(&m.half_edges, l.half_edge);
            l.face = remap(&m.faces, l.face);
        }
        for &new in m.faces.values() {
            let f = &mut self.faces[new];
            f.outer_loop = remap(&m.loops, f.outer_loop);
            for l in &mut f.inner_loops {
                *l = remap(&m.loops, *l);
            }
            f.shell = remap(&m.shells, f.shell);
        }
        for &new in m.shells.values() {
            let s = &mut self.shells[new];
            s.face = remap(&m.faces, s.face);
            s.solid = remap(&m.regions, s.solid);
        }
        for &new in m.regions.values() {
            let r = &mut self.solids[new];
            r.shell = remap(&m.shells, r.shell);
        }
        for contact in &other.topology_contacts {
            let mut moved = contact.clone();
            moved.entities = contact
                .entities
                .iter()
                .filter_map(|e| match *e {
                    ContactEntity::Vertex(v) => {
                        m.vertices.get(&v).map(|&v| ContactEntity::Vertex(v))
                    }
                    ContactEntity::Edge(x) => m.edges.get(&x).map(|&x| ContactEntity::Edge(x)),
                    ContactEntity::Face(f) => m.faces.get(&f).map(|&f| ContactEntity::Face(f)),
                })
                .collect();
            self.topology_contacts.push(moved);
        }
        m
    }
}

#[cfg(test)]
#[path = "brep_append_tests.rs"]
mod tests;
