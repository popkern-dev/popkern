//! A body's recorded contacts glued into one general body, and a
//! general body split back into manifold pieces meeting at recorded
//! contacts.
//!
//! Manifold pieces that touch (the split policy's boolean results,
//! PK-964) record each contact with the coincident entities from every
//! piece: copies of a vertex, pairs of coincident edges. Gluing them
//! (vertices first, then edges, in radial order) makes the general body
//! those pieces are; a contact of a vertex on a face's or an edge's
//! interior has no loop form here and stays recorded. Splitting is the
//! inverse: every ring of more than two fins comes apart into the
//! manifold pairs bounding its solid wedges, every vertex with more than
//! one use into one copy per use, each contact recorded as the split
//! policy records it.

use std::collections::BTreeMap;

use super::EulerError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::topology::{
    BRep, ContactEntity, EdgeHandle, FaceHandle, RegionKind, TopologyContact, TopologyContactKind,
    VertexHandle,
};

/// What [`BRep::glue_contacts`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactGlue {
    /// Contacts glued (their records gone).
    pub glued: usize,
    /// Contacts kept recorded: a vertex on a face's or an edge's
    /// interior.
    pub kept: Vec<TopologyContact>,
}

/// Why [`BRep::glue_contacts`] declined.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[error("contact at {:?} could not be glued: {error}", contact.location)]
pub struct ContactGlueError {
    /// The contact.
    pub contact: TopologyContact,
    /// The operator's decline.
    pub error: EulerError,
}

impl BRep {
    /// Glue every recorded vertex and edge contact: copies of a vertex
    /// become one vertex, coincident edges one ring in radial order. A
    /// contact naming a face (a vertex on a face's interior) stays
    /// recorded. All or nothing: on an error the body is unchanged.
    ///
    /// # Errors
    ///
    /// [`ContactGlueError`] naming the contact and the operator's
    /// decline (curves not coincident, a ring whose order is undecided,
    /// a wedge between two regions).
    pub fn glue_contacts(&mut self) -> Result<ContactGlue, ContactGlueError> {
        let mut work = self.clone();
        let contacts = std::mem::take(&mut work.topology_contacts);
        let mut vertex_alias: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
        let mut edge_alias: BTreeMap<EdgeHandle, EdgeHandle> = BTreeMap::new();
        let mut kept = Vec::new();
        let mut glued = 0;
        let fail = |contact: &TopologyContact, error| ContactGlueError {
            contact: contact.clone(),
            error,
        };
        let on_face = |c: &TopologyContact| {
            c.entities
                .iter()
                .any(|e| matches!(e, ContactEntity::Face(_)))
        };
        for contact in contacts
            .iter()
            .filter(|c| c.kind == TopologyContactKind::Vertex)
        {
            let vertices: Vec<VertexHandle> = contact
                .entities
                .iter()
                .filter_map(|e| match e {
                    ContactEntity::Vertex(v) => Some(resolve(&vertex_alias, *v)),
                    _ => None,
                })
                .collect();
            if on_face(contact) || vertices.len() < 2 {
                kept.push(contact.clone());
                continue;
            }
            for &v in &vertices[1..] {
                let keep = resolve(&vertex_alias, vertices[0]);
                if v != keep {
                    work.glue_vertices(keep, v).map_err(|e| fail(contact, e))?;
                    vertex_alias.insert(v, keep);
                }
            }
            glued += 1;
        }
        for contact in contacts
            .iter()
            .filter(|c| c.kind == TopologyContactKind::Edge)
        {
            let edges: Vec<EdgeHandle> = contact
                .entities
                .iter()
                .filter_map(|e| match e {
                    ContactEntity::Edge(x) => Some(resolve(&edge_alias, *x)),
                    _ => None,
                })
                .collect();
            if on_face(contact) || edges.len() < 2 {
                kept.push(contact.clone());
                continue;
            }
            let keep = edges[0];
            for &gone in &edges[1..] {
                if gone == keep {
                    continue;
                }
                work.glue_edge_ends(keep, gone, &mut vertex_alias)
                    .map_err(|e| fail(contact, e))?;
                work.glue_edges(keep, gone).map_err(|e| fail(contact, e))?;
                edge_alias.insert(gone, keep);
            }
            glued += 1;
        }
        // Contacts named by the faces that touch (the split policy's
        // shell findings): glue the faces' coincident boundary vertices,
        // then the edges sharing both ends; resolved when the faces now
        // share a vertex at the contact.
        let mut still = Vec::new();
        for contact in kept {
            let faces: Vec<FaceHandle> = contact
                .entities
                .iter()
                .filter_map(|e| match e {
                    ContactEntity::Face(f) => Some(*f),
                    _ => None,
                })
                .filter(|&f| work.faces.contains_key(f))
                .collect();
            if faces.len() < 2 {
                still.push(contact);
                continue;
            }
            work.glue_face_boundaries(&faces, &mut vertex_alias)
                .map_err(|e| fail(&contact, e))?;
            if work.faces_share_vertex_at(&faces, contact.location) {
                glued += 1;
            } else {
                still.push(contact);
            }
        }
        work.topology_contacts = still.clone();
        *self = work;
        Ok(ContactGlue { glued, kept: still })
    }

    /// Glue the coincident boundary vertices of `faces` (within the
    /// linear tolerance), then every two of their edges sharing both
    /// ends.
    fn glue_face_boundaries(
        &mut self,
        faces: &[FaceHandle],
        alias: &mut BTreeMap<VertexHandle, VertexHandle>,
    ) -> Result<(), EulerError> {
        let tolerance = self.tolerance.linear();
        let vertices = |b: &BRep| -> Vec<VertexHandle> {
            let mut out: Vec<VertexHandle> = faces
                .iter()
                .flat_map(|&f| b.face_fins(f))
                .map(|h| b.half_edges[h].origin)
                .collect();
            out.sort_unstable();
            out.dedup();
            out
        };
        let vs = vertices(self);
        for (i, &a) in vs.iter().enumerate() {
            for &b in &vs[i + 1..] {
                let (a, b) = (resolve(alias, a), resolve(alias, b));
                if a == b || !self.vertices.contains_key(a) || !self.vertices.contains_key(b) {
                    continue;
                }
                if self.vertices[a]
                    .position
                    .distance_to(&self.vertices[b].position)
                    <= tolerance
                {
                    self.glue_vertices(a, b)?;
                    alias.insert(b, a);
                }
            }
        }
        loop {
            let edges: Vec<EdgeHandle> = {
                let mut e: Vec<EdgeHandle> = faces
                    .iter()
                    .flat_map(|&f| self.face_fins(f))
                    .map(|h| self.half_edges[h].edge)
                    .collect();
                e.sort_unstable();
                e.dedup();
                e
            };
            let ends = |e: EdgeHandle| {
                self.edge_vertices(e)
                    .map(|(a, b)| if a < b { (a, b) } else { (b, a) })
            };
            let pair = edges.iter().enumerate().find_map(|(i, &a)| {
                edges[i + 1..]
                    .iter()
                    .copied()
                    .find(|&b| ends(a).is_some() && ends(a) == ends(b))
                    .map(|b| (a, b))
            });
            let Some((a, b)) = pair else {
                return Ok(());
            };
            self.glue_edges(a, b)?;
        }
    }

    /// Whether every one of `faces` has a vertex at `p` that all of them
    /// share.
    fn faces_share_vertex_at(&self, faces: &[FaceHandle], p: Point3) -> bool {
        let tolerance = self.tolerance.linear();
        let at = |f: FaceHandle| -> std::collections::BTreeSet<VertexHandle> {
            self.face_fins(f)
                .into_iter()
                .map(|h| self.half_edges[h].origin)
                .filter(|&v| self.vertices[v].position.distance_to(&p) <= tolerance)
                .collect()
        };
        let mut sets = faces.iter().map(|&f| at(f));
        let Some(first) = sets.next() else {
            return false;
        };
        let common = sets.fold(first, |acc, s| acc.intersection(&s).copied().collect());
        !common.is_empty()
    }

    /// Every fin of `face`'s loops.
    fn face_fins(&self, face: FaceHandle) -> Vec<crate::topology::HalfEdgeHandle> {
        let Some(f) = self.faces.get(face) else {
            return Vec::new();
        };
        std::iter::once(f.outer_loop)
            .chain(f.inner_loops.iter().copied())
            .flat_map(|l| self.loop_cycle(l))
            .collect()
    }

    /// Glue `gone`'s end vertices onto `keep`'s, each to the nearer.
    pub(crate) fn glue_edge_ends(
        &mut self,
        keep: EdgeHandle,
        gone: EdgeHandle,
        alias: &mut BTreeMap<VertexHandle, VertexHandle>,
    ) -> Result<(), EulerError> {
        let (a, b) = self
            .edge_vertices(keep)
            .ok_or(EulerError::UnknownEdge(keep))?;
        let (c, d) = self
            .edge_vertices(gone)
            .ok_or(EulerError::UnknownEdge(gone))?;
        let at = |v: VertexHandle| self.vertices[v].position;
        let straight = at(a).distance_to(&at(c)) + at(b).distance_to(&at(d))
            <= at(a).distance_to(&at(d)) + at(b).distance_to(&at(c));
        let pairs = if straight {
            [(a, c), (b, d)]
        } else {
            [(a, d), (b, c)]
        };
        for (k, g) in pairs {
            if k != g && self.vertices.contains_key(g) {
                self.glue_vertices(k, g)?;
                alias.insert(g, k);
            }
        }
        Ok(())
    }

    /// Split this general body into manifold pieces meeting at recorded
    /// contacts, the inverse of [`BRep::glue_contacts`]: each ring of
    /// more than two fins into the pairs bounding its solid wedges,
    /// then each vertex with more than one use into a copy per use.
    /// Returns the contacts recorded. On an error the body is unchanged.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotSplittable`] for a body with wire edges, acorns
    /// or a face between two solids (no manifold pieces hold those), or
    /// a ring whose fins do not pair off across solid wedges; an
    /// unglue's own decline.
    pub fn split_to_manifold(&mut self) -> Result<Vec<TopologyContact>, EulerError> {
        let unsplittable = |reason: &'static str| EulerError::NotSplittable { reason };
        if self
            .shells
            .values()
            .any(|s| !s.wire_edges.is_empty() || s.acorn.is_some())
        {
            return Err(unsplittable("the body holds wire edges or acorns"));
        }
        let solid = |b: &BRep, s| {
            b.shell_region(s).and_then(|r| b.region_kind(r)) == Some(RegionKind::Solid)
        };
        if self
            .faces
            .values()
            .any(|f| f.back != f.front && solid(self, f.back) && solid(self, f.front))
        {
            return Err(unsplittable("a face lies between two solids"));
        }
        let mut work = self.clone();
        let mut contacts = Vec::new();
        let long: Vec<EdgeHandle> = work
            .edges
            .keys()
            .filter(|&e| work.ring_summary(e).count > 2)
            .collect();
        for edge in long {
            let wedges = work
                .edge_wedges(edge)
                .ok_or(unsplittable("a fin's sense along its edge is undecidable"))?;
            let pairs: Vec<_> = wedges
                .iter()
                .filter(|w| {
                    let face = work.fin_face(w.from);
                    face.and_then(|f| work.face_side_shell(f, w.from_side))
                        .is_some_and(|s| solid(&work, s))
                })
                .map(|w| [w.from, w.to])
                .collect();
            let paired: usize = pairs.len() * 2;
            if paired != wedges.len() || pairs.iter().any(|p| p[0] == p[1]) {
                return Err(unsplittable(
                    "a ring's fins do not pair off across solid wedges",
                ));
            }
            let first = work.edges[edge].half_edge;
            let mut copies = vec![edge];
            for pair in pairs.iter().filter(|p| !p.contains(&first)) {
                copies.push(work.unglue_edge(edge, pair)?);
            }
            let location = work.edge_midpoint(edge);
            for (i, &a) in copies.iter().enumerate() {
                for &b in &copies[i + 1..] {
                    contacts.push(TopologyContact {
                        kind: TopologyContactKind::Edge,
                        location,
                        entities: vec![ContactEntity::Edge(a), ContactEntity::Edge(b)],
                    });
                }
            }
        }
        let shared: Vec<VertexHandle> = work
            .vertices
            .keys()
            .filter(|&v| work.vertex_uses(v).len() > 1)
            .collect();
        for v in shared {
            let mut copies = vec![v];
            copies.extend(work.unglue_vertex(v)?);
            contacts.push(TopologyContact {
                kind: TopologyContactKind::Vertex,
                location: work.vertices[v].position,
                entities: copies.into_iter().map(ContactEntity::Vertex).collect(),
            });
        }
        work.topology_contacts.extend(contacts.iter().cloned());
        *self = work;
        Ok(contacts)
    }

    /// A point of `edge`: its curve's midpoint, else its chord's.
    fn edge_midpoint(&self, edge: EdgeHandle) -> Point3 {
        if let Some(c) = self.edges.get(edge).and_then(|e| e.curve.as_ref()) {
            return c.geometry.point_at(0.5 * (c.t_start + c.t_end));
        }
        self.edge_vertices(edge)
            .map_or(Point3::new(0.0, 0.0, 0.0), |(a, b)| {
                let (a, b) = (self.vertices[a].position, self.vertices[b].position);
                a + (b - a) * 0.5
            })
    }
}

/// Follow `alias` from `h` to the handle it was glued onto.
fn resolve<H: Copy + Ord>(alias: &BTreeMap<H, H>, mut h: H) -> H {
    while let Some(&next) = alias.get(&h) {
        h = next;
    }
    h
}
