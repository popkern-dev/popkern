//! imprint_edge / kill_imprinted_edge: an edge made across a face's
//! interior between two of its vertices, the face's surface unchanged
//! (PK-1117).
//!
//! A contact curve where one piece's edge, or a tangent ruling, lies on
//! another piece's face is held as topology by an edge of that face
//! along the curve, which the glue then joins with the other piece's
//! (`BRep::glue_touching`). The ends are vertices of the face already:
//! on its boundary loops (split there first) or vertex loops. By where
//! they stand:
//!
//! - both on one boundary loop (the outer loop): the face splits in two
//!   (`mef`), holes and vertex loops going to the side holding them;
//! - on the outer loop and a hole: the hole joins the outer loop
//!   (`mekr`);
//! - on a boundary loop and a vertex loop: the edge is a spur of the
//!   loop, the vertex loop absorbed;
//! - on two vertex loops: a slit, one inner loop of the edge's two fins.
//!
//! Counts move as those Euler steps do (a spur and a slit: edges +1,
//! fins +2, loops −1). `kill_imprinted_edge` is the inverse of each.

use super::EulerError;
use crate::geometry::curve::CurveOnEdge;
use crate::topology::{
    BRep, Edge, EdgeHandle, FaceHandle, HalfEdge, HalfEdgeHandle, LoopHandle, RadialLinks,
    VertexHandle,
};
use crate::validation::FaceUvContainment;

/// What [`BRep::imprint_edge`] made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImprintedEdge {
    /// The new edge, its canonical fin running from the first vertex.
    pub edge: EdgeHandle,
    /// The face split off when the edge crossed the face between two
    /// points of its outer loop.
    pub new_face: Option<FaceHandle>,
}

/// Where a vertex stands on a face.
#[derive(Debug, Clone, Copy)]
enum Stand {
    /// On a boundary loop, leaving by this fin.
    Fin(HalfEdgeHandle, LoopHandle),
    /// At this vertex loop.
    Lone(LoopHandle),
}

impl BRep {
    /// Imprint an edge along `curve` across `face` from `v0` to `v1`,
    /// both vertices of the face (on a boundary loop or a vertex loop).
    /// The curve runs from `v0` to `v1` and must lie on the face's
    /// surface; the caller decides that. On an error the body is
    /// unchanged.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotImprintable`] naming why.
    pub fn imprint_edge(
        &mut self,
        face: FaceHandle,
        v0: VertexHandle,
        v1: VertexHandle,
        curve: CurveOnEdge,
    ) -> Result<ImprintedEdge, EulerError> {
        let mut work = self.clone();
        let made = work.apply_imprint(face, v0, v1, curve)?;
        *self = work;
        Ok(made)
    }

    fn apply_imprint(
        &mut self,
        face: FaceHandle,
        v0: VertexHandle,
        v1: VertexHandle,
        curve: CurveOnEdge,
    ) -> Result<ImprintedEdge, EulerError> {
        let decline = |reason: &str| EulerError::NotImprintable {
            face,
            reason: reason.into(),
        };
        let core = |e: crate::error::CoreError| EulerError::NotImprintable {
            face,
            reason: e.to_string(),
        };
        let Some(outer) = self.faces.get(face).map(|f| f.outer_loop) else {
            return Err(decline("unknown face"));
        };
        if v0 == v1 {
            return Err(decline("the two ends are one vertex"));
        }
        let (s0, s1) = (self.stand(face, v0)?, self.stand(face, v1)?);
        let (edge, from_v0, new_face) = match (s0, s1) {
            (Stand::Fin(h0, l0), Stand::Fin(h1, l1)) if l0 == l1 => {
                if l0 != outer {
                    return Err(decline("both ends lie on one hole"));
                }
                let m = self.mef(h0, h1).map_err(core)?;
                (m.edge, m.half_edge, Some(m.face))
            }
            (Stand::Fin(h0, l0), Stand::Fin(h1, _)) => {
                let (ho, hi) = if l0 == outer { (h0, h1) } else { (h1, h0) };
                if self.half_edges[ho].loop_ != outer {
                    return Err(decline("the ends lie on two holes"));
                }
                self.mekr(ho, hi).map_err(core)?;
                let e = self
                    .edge_between(face, v0, v1)
                    .ok_or(decline("mekr made no edge"))?;
                let h = self.fin_from(e, v0).ok_or(decline("mekr made no edge"))?;
                (e, h, None)
            }
            (Stand::Fin(h, _), Stand::Lone(l)) => {
                let (out, _) = self.spur(h, l);
                (self.half_edges[out].edge, out, None)
            }
            (Stand::Lone(l), Stand::Fin(h, _)) => {
                let (out, back) = self.spur(h, l);
                (self.half_edges[out].edge, back, None)
            }
            (Stand::Lone(l0), Stand::Lone(l1)) => {
                let h = self.slit(face, l0, l1);
                (self.half_edges[h].edge, h, None)
            }
        };
        self.edges[edge].half_edge = from_v0;
        self.edges[edge].curve = Some(curve);
        self.recompute_face_pcurves(face).map_err(core)?;
        if let Some(new) = new_face {
            self.recompute_face_pcurves(new).map_err(core)?;
            self.rehome_loops(face, new);
        }
        Ok(ImprintedEdge { edge, new_face })
    }

    /// Where `v` stands on `face`: the one fin of a boundary loop
    /// leaving it, or the vertex loop at it.
    fn stand(&self, face: FaceHandle, v: VertexHandle) -> Result<Stand, EulerError> {
        let decline = |reason: &str| EulerError::NotImprintable {
            face,
            reason: reason.into(),
        };
        let leaving: Vec<HalfEdgeHandle> = self
            .face_fins(face)
            .into_iter()
            .filter(|&h| self.half_edges[h].origin == v)
            .collect();
        match leaving.as_slice() {
            [h] => Ok(Stand::Fin(*h, self.half_edges[*h].loop_)),
            [] => self
                .face_vertex_loops(face)
                .into_iter()
                .find(|&l| self.loops[l].vertex == Some(v))
                .map(Stand::Lone)
                .ok_or(decline("a vertex is not on the face")),
            _ => Err(decline("a vertex is visited twice by the face's loops")),
        }
    }

    /// A new edge from vertex loop `lone`'s vertex into the boundary
    /// loop at `h`'s origin, spliced before `h`; the vertex loop goes.
    /// Returns the fin leaving `h`'s origin and the one coming back.
    fn spur(&mut self, h: HalfEdgeHandle, lone: LoopHandle) -> (HalfEdgeHandle, HalfEdgeHandle) {
        let (at, the_loop, prev) = {
            let d = &self.half_edges[h];
            (d.origin, d.loop_, d.prev)
        };
        let Some(far) = self.loops[lone].vertex else {
            return (h, h);
        };
        let edge = self.new_edge();
        let out = self.half_edges.insert(HalfEdge {
            origin: at,
            ring: RadialLinks::alone(),
            next: Default::default(),
            prev,
            edge,
            loop_: the_loop,
            pcurve: None,
        });
        let back = self.half_edges.insert(HalfEdge {
            origin: far,
            ring: RadialLinks::pair(out),
            next: h,
            prev: out,
            edge,
            loop_: the_loop,
            pcurve: None,
        });
        self.half_edges[out].set_mate(back);
        self.half_edges[out].next = back;
        self.half_edges[prev].next = out;
        self.half_edges[h].prev = back;
        self.edges[edge].half_edge = out;
        if !self.half_edges.contains_key(self.vertices[far].outgoing) {
            self.vertices[far].outgoing = back;
        }
        self.loops.remove(lone);
        (out, back)
    }

    /// A new edge between vertex loops `l0` and `l1` of `face`: `l0`
    /// becomes an inner loop of the edge's two fins, `l1` goes. Returns
    /// the fin leaving `l0`'s vertex.
    fn slit(&mut self, face: FaceHandle, l0: LoopHandle, l1: LoopHandle) -> HalfEdgeHandle {
        let (Some(a), Some(b)) = (self.loops[l0].vertex, self.loops[l1].vertex) else {
            return Default::default();
        };
        let edge = self.new_edge();
        let out = self.half_edges.insert(HalfEdge {
            origin: a,
            ring: RadialLinks::alone(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_: l0,
            pcurve: None,
        });
        let back = self.half_edges.insert(HalfEdge {
            origin: b,
            ring: RadialLinks::pair(out),
            next: out,
            prev: out,
            edge,
            loop_: l0,
            pcurve: None,
        });
        let d = &mut self.half_edges[out];
        d.set_mate(back);
        d.next = back;
        d.prev = back;
        self.edges[edge].half_edge = out;
        for (v, h) in [(a, out), (b, back)] {
            if !self.half_edges.contains_key(self.vertices[v].outgoing) {
                self.vertices[v].outgoing = h;
            }
        }
        self.loops[l0].half_edge = out;
        self.loops[l0].vertex = None;
        self.loops.remove(l1);
        self.faces[face].inner_loops.push(l0);
        out
    }

    fn new_edge(&mut self) -> EdgeHandle {
        self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        })
    }

    /// The edge of `face` between `a` and `b`, if one.
    fn edge_between(
        &self,
        face: FaceHandle,
        a: VertexHandle,
        b: VertexHandle,
    ) -> Option<EdgeHandle> {
        self.face_fins(face).into_iter().find_map(|h| {
            let e = self.half_edges[h].edge;
            let (p, q) = self.edge_vertices(e)?;
            ((p, q) == (a, b) || (p, q) == (b, a)).then_some(e)
        })
    }

    /// The fin of `edge` leaving `v`.
    fn fin_from(&self, edge: EdgeHandle, v: VertexHandle) -> Option<HalfEdgeHandle> {
        self.radial_fins(edge)
            .into_iter()
            .find(|&h| self.half_edges[h].origin == v)
    }

    /// After `face` split off `new`: each hole and vertex loop of `face`
    /// standing inside `new` moves there.
    fn rehome_loops(&mut self, face: FaceHandle, new: FaceHandle) {
        let inside = |b: &BRep, p: &crate::geometry::Point3| {
            b.face_point_trim_containment(new, p)
                .is_some_and(|c| c.containment == FaceUvContainment::Interior)
        };
        let holes = self.faces[face].inner_loops.clone();
        let mut moved = Vec::new();
        for l in holes {
            let Some(p) = self
                .loop_cycle(l)
                .first()
                .map(|&h| self.vertices[self.half_edges[h].origin].position)
            else {
                continue;
            };
            if inside(self, &p) {
                moved.push(l);
            }
        }
        for l in moved {
            self.faces[face].inner_loops.retain(|&x| x != l);
            self.faces[new].inner_loops.push(l);
            self.loops[l].face = new;
        }
        for l in self.face_vertex_loops(face) {
            let Some(p) = self.loops[l].vertex.map(|v| self.vertices[v].position) else {
                continue;
            };
            if inside(self, &p) {
                self.loops[l].face = new;
            }
        }
    }

    /// Remove `edge`, an imprinted edge, the inverse of
    /// [`BRep::imprint_edge`]: two faces it separates merge (`kef`), a
    /// hole it joined to the outer loop parts again (`kemr`), a spur's
    /// far end and a slit's ends stand as vertex loops again. On an
    /// error the body is unchanged.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotAnImprintedEdge`] for an edge of other than two
    /// fins or one whose faces' surfaces differ.
    pub fn kill_imprinted_edge(&mut self, edge: EdgeHandle) -> Result<(), EulerError> {
        let decline = |reason: &str| EulerError::NotAnImprintedEdge {
            edge,
            reason: reason.into(),
        };
        let fins = self.radial_fins(edge);
        let [a, b] = fins.as_slice() else {
            return Err(decline("the edge does not have two fins"));
        };
        let (a, b) = (*a, *b);
        let (la, lb) = (self.half_edges[a].loop_, self.half_edges[b].loop_);
        let (fa, fb) = (self.loops[la].face, self.loops[lb].face);
        let core = |e: crate::error::CoreError| EulerError::NotAnImprintedEdge {
            edge,
            reason: e.to_string(),
        };
        let mut work = self.clone();
        let mut survivor = fa;
        if la != lb {
            let same = |f: FaceHandle| work.faces[f].surface.as_ref().map(|s| &s.geometry);
            if same(fa) != same(fb) {
                return Err(decline("its faces lie on different surfaces"));
            }
            let kept = work.kef(a).map_err(core)?.surviving_loop;
            survivor = work.loops[kept].face;
        } else {
            let (na, nb) = (work.half_edges[a].next, work.half_edges[b].next);
            match (na == b, nb == a) {
                (true, true) => work.unslit(la, a, b),
                (true, false) => work.unspur(a, b),
                (false, true) => work.unspur(b, a),
                (false, false) => {
                    work.kemr(a).map_err(core)?;
                }
            }
        }
        work.recompute_face_pcurves(survivor).map_err(core)?;
        *self = work;
        Ok(())
    }

    /// Remove spur `out` (followed by `back`, which leaves the far end):
    /// the far end stands as a vertex loop again.
    fn unspur(&mut self, out: HalfEdgeHandle, back: HalfEdgeHandle) {
        let (prev, the_loop, edge) = {
            let d = &self.half_edges[out];
            (d.prev, d.loop_, d.edge)
        };
        let (far, next) = (self.half_edges[back].origin, self.half_edges[back].next);
        self.half_edges[prev].next = next;
        self.half_edges[next].prev = prev;
        if self.loops[the_loop].half_edge == out || self.loops[the_loop].half_edge == back {
            self.loops[the_loop].half_edge = next;
        }
        let face = self.loops[the_loop].face;
        self.drop_edge(edge, out, back);
        self.loops.insert(crate::topology::Loop {
            half_edge: Default::default(),
            face,
            vertex: Some(far),
        });
    }

    /// Remove slit `a`/`b`, loop `l`'s only fins: both ends stand as
    /// vertex loops.
    fn unslit(&mut self, l: LoopHandle, a: HalfEdgeHandle, b: HalfEdgeHandle) {
        let (va, vb, edge) = (
            self.half_edges[a].origin,
            self.half_edges[b].origin,
            self.half_edges[a].edge,
        );
        let face = self.loops[l].face;
        self.faces[face].inner_loops.retain(|&x| x != l);
        self.drop_edge(edge, a, b);
        self.loops[l].half_edge = Default::default();
        self.loops[l].vertex = Some(va);
        self.loops.insert(crate::topology::Loop {
            half_edge: Default::default(),
            face,
            vertex: Some(vb),
        });
    }

    /// Remove `edge` and its fins `a`, `b`, re-anchoring a vertex that
    /// left by one of them.
    fn drop_edge(&mut self, edge: EdgeHandle, a: HalfEdgeHandle, b: HalfEdgeHandle) {
        let ends = [self.half_edges[a].origin, self.half_edges[b].origin];
        // Both at once, so neither end is re-anchored onto the other fin.
        self.reanchor_vertex_outgoing(&[a, b]);
        self.remove_half_edge_reanchored(a);
        self.remove_half_edge_reanchored(b);
        self.edges.remove(edge);
        // An end left with no fin stands at a vertex loop, anchored on
        // nothing (as kev leaves the last vertex).
        for v in ends {
            if !self.half_edges.contains_key(self.vertices[v].outgoing) {
                self.vertices[v].outgoing = Default::default();
            }
        }
    }
}

#[cfg(test)]
#[path = "imprint_edge_tests.rs"]
mod tests;
