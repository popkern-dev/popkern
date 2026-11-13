//! Gluing touching pieces into one general body from their geometry
//! alone (PK-1117): no record of where they touch is read.
//!
//! [`BRep::glue_touching`] finds where the pieces touch
//! (`touch_detect`), gives every piece a vertex at each contact point
//! (its vertex there, its edge split there, or a vertex loop of its
//! face), imprints each contact curve on the faces it crosses
//! ([`BRep::imprint_edge`]), and glues what now coincides
//! ([`BRep::glue_coincident`]): one vertex at each point, one ring in
//! radial order along each curve. Space the glued pieces close in all
//! round becomes a bounded void ([`BRep::fill_enclosed_cells`]).

use std::collections::BTreeMap;

use super::{CoincidentGlue, EulerError};
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// What [`BRep::glue_touching`] did.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TouchGlue {
    /// Contact points found.
    pub points: usize,
    /// Contact curves found.
    pub curves: usize,
    /// Edges split where a point touched their interior.
    pub edges_split: usize,
    /// Vertex loops made where a point touched a face's interior.
    pub vertex_loops: usize,
    /// Edges imprinted along contact curves.
    pub imprinted: usize,
    /// What the final glue joined.
    pub glued: CoincidentGlue,
    /// Bounded voids the glued pieces enclose, each made a region
    /// ([`BRep::fill_enclosed_cells`]): the space between pieces that
    /// touch all round it (the cell inside two of three crossing rods
    /// when exactly one is kept, PK-1118).
    pub enclosed: usize,
}

impl BRep {
    /// Make one general body of this body's touching pieces (face sets
    /// joined through shared edges), found from their geometry: each
    /// contact point a vertex every touching piece uses, each contact
    /// curve an edge every touching piece uses, glued. Wire edges and
    /// acorns take no part. A body whose pieces touch nowhere is
    /// unchanged. All or nothing: on an error the body is unchanged.
    ///
    /// # Errors
    ///
    /// The operators' declines: an edge split, a vertex loop, an imprint
    /// ([`EulerError::NotImprintable`]), or the glue's (a ring whose
    /// order is undecided, a wedge between two regions).
    pub fn glue_touching(&mut self) -> Result<TouchGlue, EulerError> {
        let mut work = self.clone();
        let plan = work.touch_plan();
        let tol = work.tolerance.linear_reach();
        let mut out = TouchGlue {
            points: plan.points.len(),
            curves: plan.curves.len(),
            ..TouchGlue::default()
        };
        let by_piece = |b: &BRep| -> Vec<Vec<FaceHandle>> {
            let mut groups: BTreeMap<usize, Vec<FaceHandle>> = BTreeMap::new();
            for (f, k) in b.face_pieces() {
                groups.entry(k).or_default().push(f);
            }
            groups.into_values().collect()
        };
        // Points: splits and vertex loops leave the pieces as they were.
        let pieces = by_piece(&work);
        for p in &plan.points {
            for faces in &pieces {
                work.vertex_for(faces, *p, tol, &mut out)?;
            }
        }
        // Curves: each imprint may split a face, so pieces are re-read.
        for c in &plan.curves {
            let (p0, p1) = (c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end));
            let m = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            for faces in by_piece(&work) {
                // Every face of the piece the curve lies across: one
                // for a piece touching another, two for a piece
                // touching itself.
                let across: Vec<FaceHandle> = faces
                    .iter()
                    .copied()
                    .filter(|&f| work.face_interior_holds(f, &m, tol))
                    .collect();
                for f in across {
                    let v0 = work.vertex_at_curve_end(f, p0, tol, &mut out)?;
                    let v1 = work.vertex_at_curve_end(f, p1, tol, &mut out)?;
                    work.imprint_edge(f, v0, v1, c.clone())?;
                    out.imprinted += 1;
                }
            }
        }
        work.contain_vertex_reach();
        out.glued = work.glue_coincident()?;
        if out.glued != CoincidentGlue::default() {
            out.enclosed = work.fill_enclosed_cells().len();
        }
        *self = work;
        Ok(out)
    }

    /// The vertex of face `f` at a contact curve's end `p`: one within
    /// `tol`, else the nearest boundary edge of `f` split at `p` when
    /// that edge lies within two tolerances of it. The curve lies on the
    /// other piece, `p` with it, and this face is the gap away, which
    /// the detection accepted within the tolerance: its boundary there
    /// can lie the gap and its own placement off `p`. A cap a tolerance
    /// short of a ring's extreme found its touch 1.00000000003e-6 off
    /// the ring, whose arc, built 1.7e-11 off its axis, lay 1.0000041e-6
    /// from the curve's end, and the glue declined (PK-1170).
    fn vertex_at_curve_end(
        &mut self,
        f: FaceHandle,
        p: Point3,
        tol: f64,
        out: &mut TouchGlue,
    ) -> Result<VertexHandle, EulerError> {
        let at = |b: &BRep| {
            b.face_vertices(f)
                .into_iter()
                .find(|&v| b.vertices[v].position.distance_to(&p) <= tol)
        };
        if let Some(v) = at(self) {
            return Ok(v);
        }
        let reach = 2.0 * self.tolerance.linear_reach();
        let vertex_near = self
            .face_vertices(f)
            .into_iter()
            .map(|v| (v, self.vertices[v].position.distance_to(&p)))
            .filter(|&(_, d)| d <= reach)
            .min_by(|a, b| a.1.total_cmp(&b.1));
        if let Some((v, d)) = vertex_near {
            self.record_vertex_reach(v, d);
            return Ok(v);
        }
        let edge_distance = |e: EdgeHandle| -> f64 {
            self.edges
                .get(e)
                .and_then(|d| d.curve.as_ref())
                .and_then(|c| {
                    let t = c.geometry.closest_point(&p)?;
                    Some(c.geometry.point_at(t).distance_to(&p))
                })
                .unwrap_or(f64::INFINITY)
        };
        let nearest = self
            .face_edges(f)
            .into_iter()
            .filter(|&e| self.edge_interior_holds(e, &p, reach))
            .min_by(|&a, &b| edge_distance(a).total_cmp(&edge_distance(b)));
        let missing = || EulerError::NotImprintable {
            face: f,
            reason: format!("no vertex of the face at the contact curve's end {p:?}"),
        };
        let edge = nearest.ok_or_else(missing)?;
        // Split at the edge's own point nearest `p` (the edge is the gap
        // away; `p` itself is off it).
        let on_edge = self
            .edges
            .get(edge)
            .and_then(|d| d.curve.as_ref())
            .and_then(|c| c.geometry.closest_point(&p).map(|t| c.geometry.point_at(t)))
            .ok_or_else(missing)?;
        let tolerance = self.tolerance;
        self.split_edge_at_point(edge, on_edge, &tolerance)
            .map_err(|e| EulerError::EdgeSplit {
                edge,
                reason: e.to_string(),
            })?;
        out.edges_split += 1;
        let v = self
            .face_vertices(f)
            .into_iter()
            .min_by(|&a, &b| {
                let d = |v: VertexHandle| self.vertices[v].position.distance_to(&p);
                d(a).total_cmp(&d(b))
            })
            .filter(|&v| self.vertices[v].position.distance_to(&p) <= reach)
            .ok_or_else(missing)?;
        let d = self.vertices[v].position.distance_to(&p);
        self.record_vertex_reach(v, d);
        Ok(v)
    }

    /// Every edge ending at a vertex whose recorded tolerance exceeds
    /// the body's takes at least that tolerance (an edge's tolerance
    /// contains its vertices').
    fn contain_vertex_reach(&mut self) {
        let linear = self.tolerance.linear();
        let raised: Vec<(EdgeHandle, f64)> = self
            .edges
            .keys()
            .filter_map(|e| {
                let (a, b) = self.edge_vertices(e)?;
                let t = self.vertex_tolerance(a).max(self.vertex_tolerance(b));
                (t > linear && t > self.edge_tolerance(e)).then_some((e, t))
            })
            .collect();
        for (e, t) in raised {
            if let Some(d) = self.edges.get_mut(e) {
                d.tolerance = Some(t);
            }
        }
    }

    /// Record on vertex `v` that a curve end it takes lies `d` from it:
    /// its tolerance becomes at least that, as an imported curve's end
    /// off its vertex records the file's precision (PK-639), so the
    /// contact curve imprinted there and glued closes on it.
    fn record_vertex_reach(&mut self, v: VertexHandle, d: f64) {
        if d > self.tolerance.linear()
            && let Some(vx) = self.vertices.get_mut(v)
        {
            let at = d * (1.0 + 1e-9);
            vx.tolerance = Some(vx.tolerance.map_or(at, |t| t.max(at)));
        }
    }

    /// Give the piece of faces `faces` a vertex at `p` wherever it
    /// touches `p`: each of its edges holding `p` inside is split there,
    /// and each face holding it inside, with no vertex there, gets a
    /// vertex loop. A piece touching another holds `p` once (a vertex,
    /// an edge or a face); a piece touching itself holds it on two of
    /// its edges or faces, and each gets the vertex the glue joins.
    fn vertex_for(
        &mut self,
        faces: &[FaceHandle],
        p: Point3,
        tol: f64,
        out: &mut TouchGlue,
    ) -> Result<(), EulerError> {
        let mut edges: Vec<_> = faces
            .iter()
            .flat_map(|&f| self.face_edges(f))
            .filter(|&e| self.edge_interior_holds(e, &p, tol))
            .collect();
        edges.sort_unstable();
        edges.dedup();
        for edge in edges {
            let tolerance = self.tolerance;
            self.split_edge_at_point(edge, p, &tolerance)
                .map_err(|e| EulerError::EdgeSplit {
                    edge,
                    reason: e.to_string(),
                })?;
            out.edges_split += 1;
        }
        let holding: Vec<FaceHandle> = faces
            .iter()
            .copied()
            .filter(|&f| self.faces.contains_key(f))
            .filter(|&f| {
                self.face_interior_holds(f, &p, tol)
                    && !self
                        .face_vertices(f)
                        .into_iter()
                        .any(|v| self.vertices[v].position.distance_to(&p) <= tol)
            })
            .collect();
        for f in holding {
            self.make_vertex_loop(f, p)?;
            out.vertex_loops += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "touch_tests.rs"]
mod tests;
