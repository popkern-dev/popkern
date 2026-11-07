//! Vertex loops stand on their face's interior (PK-1116).
//!
//! A loop of no fins names the lone vertex it stands at. On a face with
//! edges (a vertex loop, not the first Euler step's empty face) the
//! vertex must lie on the face's surface within tolerance and strictly
//! inside its trim: a vertex loop off the face, on an edge, or in a hole
//! claims a contact that is not there.

use crate::error::ValidationError;
use crate::topology::BRep;
use crate::validation::FaceUvContainment;

impl BRep {
    /// Every loop of no fins that names no live vertex, and every vertex
    /// loop whose vertex lies off its face or not inside its trim.
    pub fn validate_vertex_loops_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for (l, data) in &self.loops {
            if self.half_edges.contains_key(data.half_edge) {
                continue;
            }
            let finding = |reason: String| ValidationError::VertexLoop {
                loop_: l,
                face: data.face,
                reason,
            };
            let Some(vertex) = data.vertex.and_then(|v| self.vertices.get(v)) else {
                errors.push(finding("names no live vertex".into()));
                continue;
            };
            let Some(face) = self.faces.get(data.face) else {
                errors.push(finding("its face is gone".into()));
                continue;
            };
            // The first Euler step's empty face has no trim to stand in.
            if face.outer_loop == l || face.surface.is_none() {
                continue;
            }
            let Some(place) = self.face_point_trim_containment(data.face, &vertex.position) else {
                errors.push(finding(
                    "the face's surface or trim cannot be read there".into(),
                ));
                continue;
            };
            let tolerance = self
                .tolerance
                .linear()
                .max(face.tolerance.unwrap_or(0.0))
                .max(vertex.tolerance.unwrap_or(0.0));
            if place.distance > tolerance {
                errors.push(finding(format!(
                    "its vertex lies {:.3e} off the surface (tolerance {tolerance:.3e})",
                    place.distance
                )));
            } else if place.containment != FaceUvContainment::Interior {
                errors.push(finding(format!(
                    "its vertex is {:?} the trim, not inside it",
                    place.containment
                )));
            }
        }
        errors
    }
}
