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
use crate::topology::{BRep, FaceHandle, VertexHandle};

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
        let tol = work.tolerance.linear();
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
                let Some(f) = faces
                    .iter()
                    .copied()
                    .find(|&f| work.face_interior_holds(f, &m, tol))
                else {
                    continue;
                };
                let end = |b: &BRep, p: Point3| -> Result<VertexHandle, EulerError> {
                    b.face_vertices(f)
                        .into_iter()
                        .find(|&v| b.vertices[v].position.distance_to(&p) <= tol)
                        .ok_or_else(|| EulerError::NotImprintable {
                            face: f,
                            reason: format!(
                                "no vertex of the face at the contact curve's end {p:?}"
                            ),
                        })
                };
                let (v0, v1) = (end(&work, p0)?, end(&work, p1)?);
                work.imprint_edge(f, v0, v1, c.clone())?;
                out.imprinted += 1;
            }
        }
        out.glued = work.glue_coincident()?;
        if out.glued != CoincidentGlue::default() {
            out.enclosed = work.fill_enclosed_cells().len();
        }
        *self = work;
        Ok(out)
    }

    /// Give the piece of faces `faces` a vertex at `p` when it touches
    /// `p`: its vertex there, else its edge split there, else a vertex
    /// loop of its face.
    fn vertex_for(
        &mut self,
        faces: &[FaceHandle],
        p: Point3,
        tol: f64,
        out: &mut TouchGlue,
    ) -> Result<(), EulerError> {
        let has_vertex = faces.iter().any(|&f| {
            self.face_vertices(f)
                .into_iter()
                .any(|v| self.vertices[v].position.distance_to(&p) <= tol)
        });
        if has_vertex {
            return Ok(());
        }
        let edge = faces
            .iter()
            .flat_map(|&f| self.face_edges(f))
            .find(|&e| self.edge_interior_holds(e, &p, tol));
        if let Some(edge) = edge {
            let tolerance = self.tolerance;
            self.split_edge_at_point(edge, p, &tolerance)
                .map_err(|e| EulerError::EdgeSplit {
                    edge,
                    reason: e.to_string(),
                })?;
            out.edges_split += 1;
            return Ok(());
        }
        if let Some(&f) = faces
            .iter()
            .find(|&&f| self.face_interior_holds(f, &p, tol))
        {
            self.make_vertex_loop(f, p)?;
            out.vertex_loops += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "touch_tests.rs"]
mod tests;
