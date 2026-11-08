//! A body's pieces taken out whole (PK-1115): each kept piece with the
//! faces, edges and vertices bounding it, their handles, names and
//! geometry unchanged.
//!
//! A body whose pieces share no entity (disjoint lumps, the pieces
//! [`BRep::separate_at_contacts`] leaves, a shell imported as pieces)
//! comes apart into a body per piece without anything being cut;
//! pieces sharing a vertex or edge are one piece, a general body.

use std::collections::{BTreeMap, HashSet};

use super::EulerError;
use crate::topology::{BRep, EdgeHandle, FaceHandle, RegionHandle, VertexHandle};

impl BRep {
    /// This body with only the solids in `keep`: every entity bounding a
    /// dropped solid is removed, every kept one keeps its handle, name
    /// and geometry.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotSplittable`] when a kept and a dropped solid
    /// share a half-edge, edge or vertex.
    pub fn retain_pieces(&self, keep: &[RegionHandle]) -> Result<BRep, EulerError> {
        let solids: HashSet<_> = keep
            .iter()
            .copied()
            .filter(|s| self.solid_regions().contains(s))
            .collect();
        let shells: HashSet<_> = solids.iter().flat_map(|&s| self.region_shells(s)).collect();
        let faces: HashSet<_> = shells.iter().flat_map(|&s| self.shell_faces(s)).collect();
        let loops: HashSet<_> = self
            .loops
            .iter()
            .filter(|(_, l)| faces.contains(&l.face))
            .map(|(h, _)| h)
            .collect();
        let half_edges: HashSet<_> = self
            .half_edges
            .iter()
            .filter(|(_, he)| loops.contains(&he.loop_))
            .map(|(h, _)| h)
            .collect();
        let edges: HashSet<_> = half_edges
            .iter()
            .filter_map(|h| self.half_edges.get(*h).map(|he| he.edge))
            .collect();
        let lone = loops.iter().filter_map(|&l| self.loops[l].vertex);
        let vertices: HashSet<_> = half_edges
            .iter()
            .filter_map(|h| self.half_edges.get(*h).map(|he| he.origin))
            .chain(lone)
            .collect();
        let shared = |reason| EulerError::NotSplittable { reason };
        for h in &half_edges {
            let he = &self.half_edges[*h];
            let links = [he.ring.next, he.ring.prev, he.next, he.prev];
            if links
                .iter()
                .any(|l| self.half_edges.contains_key(*l) && !half_edges.contains(l))
            {
                return Err(shared("the body's solids share a half-edge"));
            }
        }
        for (h, he) in &self.half_edges {
            if half_edges.contains(&h) {
                continue;
            }
            if edges.contains(&he.edge) {
                return Err(shared("the body's solids share an edge"));
            }
            if vertices.contains(&he.origin) {
                return Err(shared("the body's solids share a vertex"));
            }
        }
        let mut out = self.clone();
        out.retain_solid_regions(&solids.iter().copied().collect::<Vec<_>>());
        out.loops.retain(|h, _| loops.contains(&h));
        out.half_edges.retain(|h, _| half_edges.contains(&h));
        out.edges.retain(|h, _| edges.contains(&h));
        out.vertices.retain(|h, _| vertices.contains(&h));
        out.surface_boundary_hints
            .retain(|h, _| half_edges.contains(h));
        out.surface_boundary_declines
            .retain(|h, _| half_edges.contains(h));
        out.surface_boundary_audit_events
            .retain(|e| half_edges.contains(&e.half_edge));
        out.face_use_pcurves
            .retain(|k, _| faces.contains(&k.face) && half_edges.contains(&k.half_edge));
        Ok(out)
    }

    /// This body's pieces: its face sets joined through a shared edge or
    /// vertex, or through a solid region they bound together (a solid's
    /// cavity shell goes with it), in face slot order.
    pub fn piece_face_sets(&self) -> Vec<Vec<FaceHandle>> {
        let faces: Vec<FaceHandle> = self.faces.keys().collect();
        let index: BTreeMap<FaceHandle, usize> =
            faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let join = |parent: &mut Vec<usize>, a: usize, b: usize| {
            let (a, b) = (root(parent, a), root(parent, b));
            parent[a] = b;
        };
        let mut by_vertex: BTreeMap<VertexHandle, usize> = BTreeMap::new();
        let mut by_edge: BTreeMap<EdgeHandle, usize> = BTreeMap::new();
        let mut by_region: BTreeMap<RegionHandle, usize> = BTreeMap::new();
        for (l, data) in &self.loops {
            let Some(&i) = index.get(&data.face) else {
                continue;
            };
            let fins = if self.half_edges.contains_key(data.half_edge) {
                self.loop_cycle(l)
            } else {
                Vec::new()
            };
            let vertices = fins
                .iter()
                .map(|&h| self.half_edges[h].origin)
                .chain(data.vertex);
            for v in vertices {
                match by_vertex.get(&v) {
                    Some(&j) => join(&mut parent, i, j),
                    None => {
                        by_vertex.insert(v, i);
                    }
                }
            }
            for &h in &fins {
                let e = self.half_edges[h].edge;
                match by_edge.get(&e) {
                    Some(&j) => join(&mut parent, i, j),
                    None => {
                        by_edge.insert(e, i);
                    }
                }
            }
        }
        for (&f, &i) in &index {
            let face = &self.faces[f];
            for side in [face.back, face.front] {
                let Some(r) = self.shell_region(side) else {
                    continue;
                };
                if self.region_kind(r) != Some(crate::topology::RegionKind::Solid) {
                    continue;
                }
                match by_region.get(&r) {
                    Some(&j) => join(&mut parent, i, j),
                    None => {
                        by_region.insert(r, i);
                    }
                }
            }
        }
        let mut sets: BTreeMap<usize, Vec<FaceHandle>> = BTreeMap::new();
        let mut order: Vec<usize> = Vec::new();
        for (i, &f) in faces.iter().enumerate() {
            let r = root(&mut parent, i);
            if !sets.contains_key(&r) {
                order.push(r);
            }
            sets.entry(r).or_default().push(f);
        }
        order.into_iter().filter_map(|r| sets.remove(&r)).collect()
    }

    /// Each piece of this body ([`BRep::piece_face_sets`]) as a body of
    /// its own: its faces with the loops, fins, edges and vertices that
    /// bound them, their handles, names and geometry unchanged, and the
    /// shells and regions they lie in. A body of one piece comes back
    /// whole.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotSplittable`] for a body holding wire edges or
    /// acorns (no face set holds them).
    pub fn piece_bodies(&self) -> Result<Vec<BRep>, EulerError> {
        if self
            .shells
            .values()
            .any(|s| !s.wire_edges.is_empty() || s.acorn.is_some())
        {
            return Err(EulerError::NotSplittable {
                reason: "the body holds wire edges or acorns",
            });
        }
        let sets = self.piece_face_sets();
        if sets.len() <= 1 {
            return Ok(vec![self.clone()]);
        }
        Ok(sets.iter().map(|set| self.retain_faces(set)).collect())
    }

    /// This body with only `keep`, a union of pieces: every other face
    /// goes with its loops and fins, then every edge and vertex no kept
    /// fin or vertex loop uses; shells and regions left empty go with
    /// their last face.
    fn retain_faces(&self, keep: &[FaceHandle]) -> BRep {
        let keep: HashSet<FaceHandle> = keep.iter().copied().collect();
        let mut out = self.clone();
        let dropped: Vec<FaceHandle> = out.faces.keys().filter(|f| !keep.contains(f)).collect();
        let loops: Vec<_> = out
            .loops
            .iter()
            .filter(|(_, l)| !keep.contains(&l.face))
            .map(|(h, _)| h)
            .collect();
        let fins: HashSet<_> = out
            .half_edges
            .iter()
            .filter(|(_, h)| loops.contains(&h.loop_))
            .map(|(h, _)| h)
            .collect();
        out.half_edges.retain(|h, _| !fins.contains(&h));
        for l in loops {
            out.loops.remove(l);
        }
        for f in dropped {
            out.remove_face_reanchored(f);
        }
        let live_edges: HashSet<_> = out.half_edges.values().map(|h| h.edge).collect();
        out.edges.retain(|e, _| live_edges.contains(&e));
        let live_vertices: HashSet<_> = out
            .half_edges
            .values()
            .map(|h| h.origin)
            .chain(out.loops.values().filter_map(|l| l.vertex))
            .collect();
        out.vertices.retain(|v, _| live_vertices.contains(&v));
        let live_fins: HashSet<_> = out.half_edges.keys().collect();
        out.surface_boundary_hints
            .retain(|h, _| live_fins.contains(h));
        out.surface_boundary_declines
            .retain(|h, _| live_fins.contains(h));
        out.surface_boundary_audit_events
            .retain(|e| live_fins.contains(&e.half_edge));
        out.face_use_pcurves
            .retain(|k, _| keep.contains(&k.face) && live_fins.contains(&k.half_edge));
        out
    }
}
