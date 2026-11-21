//! make_vertex_loop / kill_vertex_loop: a vertex on a face's interior,
//! held by the face as a loop of no fins (PK-1116).
//!
//! A vertex of one piece touching another piece's face away from its
//! edges (a cone's apex on a block's face) is a contact the face must
//! hold as topology for the two to be one general body: the face gets a
//! vertex loop at the point, and the apex is glued onto that vertex
//! ([`BRep::glue_vertices`]). The loop bounds nothing, so it is none of
//! the face's boundary loops (`outer_loop`, `inner_loops`: its trim is
//! unchanged and every boundary walk is too); it names its face, which
//! [`BRep::face_vertex_loops`] reads back, and its vertex is one use of
//! its own ([`crate::topology::VertexUse::Loop`]).

use super::EulerError;
use crate::geometry::Point3;
use crate::topology::{BRep, FaceHandle, Loop, LoopHandle, Vertex, VertexHandle};
use crate::validation::FaceUvContainment;

/// What [`BRep::make_vertex_loop`] made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VertexLoopResult {
    /// The vertex.
    pub vertex: VertexHandle,
    /// The face's new loop standing at it.
    pub loop_: LoopHandle,
}

impl BRep {
    /// Whether a vertex loop of `face` may stand at `position`, the one
    /// rule [`Self::make_vertex_loop`] and the validator read: on the
    /// face's surface within `tolerance`, inside its trim by even-odd
    /// parity, and farther than `tolerance` from its edges in space (the
    /// trim classifier's parametric band of ten tolerances is no
    /// boundary: a touch beside a cylinder's seam fell in it). The
    /// reason with its measurement when not.
    pub(crate) fn vertex_loop_site(
        &self,
        face: FaceHandle,
        position: &Point3,
        tolerance: f64,
    ) -> Result<(), String> {
        let Some(place) = self.face_point_trim_containment_within(face, position, 1e-12) else {
            return Err("the face's surface or trim cannot be read there".into());
        };
        // Within the tolerance is inclusive with its rounding allowance
        // (`Tolerance::within_linear`): a contact placed a tolerance off
        // read 1.000e-6 over 1.000e-6.
        if place.distance > tolerance * (1.0 + crate::tolerance::LINEAR_ROUNDING) {
            return Err(format!(
                "the point lies {:.9e} off the surface (tolerance {tolerance:.3e})",
                place.distance
            ));
        }
        if place.containment != FaceUvContainment::Interior {
            return Err(format!(
                "the point is {:?} the face's trim, not inside it",
                place.containment
            ));
        }
        if self.point_near_face_boundary(face, position, tolerance) {
            return Err("the point lies on the face's boundary within the tolerance".into());
        }
        Ok(())
    }

    /// Make a vertex at `position` on `face`'s interior, held by a new
    /// vertex loop of the face (a loop of no fins). The point must lie on the face's surface
    /// within the linear tolerance, inside its trim by even-odd parity
    /// and farther than the tolerance from its edges in space (the trim
    /// classifier's parametric band of ten tolerances is no boundary: a
    /// touch beside a cylinder's seam fell in it), and off the face's
    /// other vertex loops. Counts: vertices +1, loops +1.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotOnFaceInterior`] naming what fails, with the
    /// measured distance where there is one.
    pub fn make_vertex_loop(
        &mut self,
        face: FaceHandle,
        position: Point3,
    ) -> Result<VertexLoopResult, EulerError> {
        let decline = |reason: String| Err(EulerError::NotOnFaceInterior { face, reason });
        let Some(data) = self.faces.get(face) else {
            return decline("unknown face".into());
        };
        let tolerance = self.tolerance.linear().max(data.tolerance.unwrap_or(0.0));
        if let Err(reason) = self.vertex_loop_site(face, &position, tolerance) {
            return decline(reason);
        }
        for l in self.face_vertex_loops(face) {
            if let Some(v) = self.loops.get(l).and_then(|l| l.vertex)
                && let Some(at) = self.vertices.get(v)
                && at.position.distance_to(&position) <= tolerance
            {
                return decline(format!("vertex loop {l:?} already stands there"));
            }
        }
        let vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let loop_ = self.loops.insert(Loop {
            half_edge: Default::default(),
            face,
            vertex: Some(vertex),
        });
        Ok(VertexLoopResult { vertex, loop_ })
    }

    /// Remove vertex loop `loop_` and its vertex, the inverse of
    /// [`BRep::make_vertex_loop`]; the vertex must have no other use
    /// (unglue it first). Counts: vertices −1, loops −1.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotALoneVertexLoop`] when the loop has fins, names
    /// no vertex, is a face's outer loop, or its vertex has another use.
    pub fn kill_vertex_loop(&mut self, loop_: LoopHandle) -> Result<(), EulerError> {
        let decline = Err(EulerError::NotALoneVertexLoop { loop_ });
        let Some(data) = self.loops.get(loop_) else {
            return decline;
        };
        let (Some(vertex), false) = (data.vertex, self.half_edges.contains_key(data.half_edge))
        else {
            return decline;
        };
        if !self.face_vertex_loops(data.face).contains(&loop_)
            || self.vertex_uses(vertex).len() != 1
        {
            return decline;
        }
        self.loops.remove(loop_);
        self.vertices.remove(vertex);
        Ok(())
    }
}

#[cfg(test)]
#[path = "vertex_loop_tests.rs"]
mod tests;
