//! A general body taken apart into its manifold pieces (PK-1106,
//! PK-1115).
//!
//! Gluing touching pieces into one general body is
//! [`BRep::glue_touching`], which reads their geometry; this is the
//! inverse: every ring of more than two fins comes apart into the
//! manifold pairs bounding its solid wedges, every vertex with more
//! than one use into one copy per use. The pieces then touch without
//! sharing an entity, which no single body may hold, so
//! [`BRep::split_to_manifold`] returns them as separate bodies;
//! [`BRep::separate_at_contacts`] leaves them in one body for operators
//! that take pieces as an intermediate.

use std::collections::BTreeMap;

use super::EulerError;
use crate::topology::{BRep, EdgeHandle, RegionKind, VertexHandle};

impl BRep {
    /// This general body's manifold pieces as separate bodies, one per
    /// solid, each valid on its own; touching pieces meet where their
    /// own copies of the contact edges and vertices coincide. A body
    /// with nothing to take apart comes back as its solids (one body
    /// for a solid).
    ///
    /// # Errors
    ///
    /// As [`BRep::separate_at_contacts`].
    pub fn split_to_manifold(&self) -> Result<Vec<BRep>, EulerError> {
        let mut work = self.clone();
        work.separate_at_contacts()?;
        work.piece_bodies()
    }

    /// Take this general body apart into its manifold pieces, left in
    /// this body touching without sharing an entity: each ring of more
    /// than two fins into the pairs bounding its solid wedges, then each
    /// vertex with more than one use into a copy per use (a copy
    /// standing in a vertex loop dropped with its loop). The result is
    /// an intermediate for operators that take pieces (a boolean's
    /// general operand), not a valid body: touching lumps are one
    /// general body or separate bodies ([`BRep::split_to_manifold`]).
    /// On an error the body is unchanged.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotSplittable`] for a body with wire edges, acorns
    /// or a face between two solids (no manifold pieces hold those), or
    /// a ring whose fins do not pair off across solid wedges; an
    /// unglue's own decline.
    pub fn separate_at_contacts(&mut self) -> Result<(), EulerError> {
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
            if pairs.len() * 2 != wedges.len() || pairs.iter().any(|p| p[0] == p[1]) {
                return Err(unsplittable(
                    "a ring's fins do not pair off across solid wedges",
                ));
            }
            let first = work.edges[edge].half_edge;
            for pair in pairs.iter().filter(|p| !p.contains(&first)) {
                work.unglue_edge(edge, pair)?;
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
            for c in copies {
                if let Some((l, _)) = work.vertex_loops().into_iter().find(|&(_, w)| w == c) {
                    work.kill_vertex_loop(l)?;
                }
            }
        }
        *self = work;
        Ok(())
    }

    /// Glue `gone`'s end vertices onto `keep`'s, each to the nearer,
    /// noting each glued vertex in `alias` (sewing's edge glue).
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
}
