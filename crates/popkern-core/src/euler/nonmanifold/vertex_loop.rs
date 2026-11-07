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
    /// Make a vertex at `position` on `face`'s interior, held by a new
    /// vertex loop of the face (a loop of no fins). The point must lie on the face's surface
    /// within the linear tolerance, strictly inside its trim (off its
    /// edges by the trim classifier's band), and off the face's other
    /// vertex loops. Counts: vertices +1, loops +1.
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
        let Some(place) = self.face_point_trim_containment(face, &position) else {
            return decline("the face's surface or trim cannot be read there".into());
        };
        let tolerance = self.tolerance.linear().max(data.tolerance.unwrap_or(0.0));
        if place.distance > tolerance {
            return decline(format!(
                "the point lies {:.3e} off the surface (tolerance {tolerance:.3e})",
                place.distance
            ));
        }
        if place.containment != FaceUvContainment::Interior {
            return decline(format!(
                "the point is {:?} the face's trim, not inside it",
                place.containment
            ));
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
