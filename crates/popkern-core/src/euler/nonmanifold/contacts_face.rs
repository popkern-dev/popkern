//! Contacts named by the faces that touch (PK-1108, PK-1116): the
//! faces' coincident boundary vertices and edges glued, a vertex on
//! another face's interior held by a vertex loop of that face, one on
//! another face's edge splitting it.

use std::collections::BTreeMap;

use super::EulerError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::contacts::resolve;

impl BRep {
    /// Glue the coincident boundary vertices of `faces` (within the
    /// linear tolerance), then every two of their edges sharing both
    /// ends.
    pub(super) fn glue_face_boundaries(
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

    /// A boundary vertex of one of `faces` at `p` lying on another's
    /// interior (a cone's apex on a block's face): that face gets a
    /// vertex loop there ([`BRep::make_vertex_loop`]) and the vertex is
    /// glued onto it; lying on another's edge away from its ends, the
    /// edge is split there and the vertex glued onto the split's.
    /// Nothing when no such pair is found.
    pub(super) fn glue_vertex_into_face(
        &mut self,
        faces: &[FaceHandle],
        p: Point3,
    ) -> Result<(), EulerError> {
        let tolerance = self.tolerance.linear();
        for &f in faces {
            let near: Vec<VertexHandle> = self
                .face_fins(f)
                .into_iter()
                .map(|h| self.half_edges[h].origin)
                .filter(|&v| self.vertices[v].position.distance_to(&p) <= tolerance)
                .collect();
            for v in near {
                let at = self.vertices[v].position;
                let host = faces.iter().copied().find(|&g| {
                    g != f
                        && self.face_point_trim_containment(g, &at).is_some_and(|c| {
                            c.distance <= tolerance
                                && c.containment == crate::validation::FaceUvContainment::Interior
                        })
                });
                if let Some(g) = host {
                    let made = self.make_vertex_loop(g, at)?;
                    self.glue_vertices(v, made.vertex)?;
                    return Ok(());
                }
                // On another face's edge, away from its ends: the edge is
                // split there and the vertex glued onto the split's.
                if let Some(edge) = self.edge_touched_inside(faces, f, at) {
                    let split = self
                        .split_edge_at_point(edge, at, &self.tolerance.clone())
                        .map_err(|e| EulerError::EdgeSplit {
                            edge,
                            reason: e.to_string(),
                        })?;
                    self.glue_vertices(v, split.vertex)?;
                    return Ok(());
                }
            }
        }
        Ok(())
    }

    /// An edge of one of `faces` other than `f` whose interior passes
    /// within the linear tolerance of `at`, its ends farther.
    fn edge_touched_inside(
        &self,
        faces: &[FaceHandle],
        f: FaceHandle,
        at: Point3,
    ) -> Option<EdgeHandle> {
        let tolerance = self.tolerance.linear();
        faces
            .iter()
            .filter(|&&g| g != f)
            .flat_map(|&g| self.face_fins(g))
            .map(|h| self.half_edges[h].edge)
            .find(|&e| {
                let Some((a, b)) = self.edge_vertices(e) else {
                    return false;
                };
                let far = |v: VertexHandle| self.vertices[v].position.distance_to(&at) > tolerance;
                let Some(c) = self.edges.get(e).and_then(|d| d.curve.as_ref()) else {
                    return false;
                };
                let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                far(a)
                    && far(b)
                    && c.geometry.closest_point(&at).is_some_and(|t| {
                        t > lo && t < hi && c.geometry.point_at(t).distance_to(&at) <= tolerance
                    })
            })
    }

    /// Whether every one of `faces` has a vertex at `p` that all of them
    /// share (a vertex loop of the face counting).
    pub(super) fn faces_share_vertex_at(&self, faces: &[FaceHandle], p: Point3) -> bool {
        let tolerance = self.tolerance.linear();
        let at = |f: FaceHandle| -> std::collections::BTreeSet<VertexHandle> {
            let lone = self
                .face_vertex_loops(f)
                .into_iter()
                .filter_map(|l| self.loops.get(l).and_then(|l| l.vertex));
            self.face_fins(f)
                .into_iter()
                .map(|h| self.half_edges[h].origin)
                .chain(lone)
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
    pub(super) fn face_fins(&self, face: FaceHandle) -> Vec<crate::topology::HalfEdgeHandle> {
        let Some(f) = self.faces.get(face) else {
            return Vec::new();
        };
        std::iter::once(f.outer_loop)
            .chain(f.inner_loops.iter().copied())
            .flat_map(|l| self.loop_cycle(l))
            .collect()
    }
}
