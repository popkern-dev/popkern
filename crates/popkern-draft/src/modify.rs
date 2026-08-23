//! In-place draft-angle modification (PK-280 stage 2).
//!
//! Draft is already topology-preserving by construction — the
//! operation tilts surfaces about the neutral line and recomputes the
//! touched vertices and edges. Rotations about the SAME neutral line
//! compose additively, and the neutral line itself is stationary
//! under the rotation (it lies on the face before and after), so
//! changing a drafted face's angle from `old` to `new` is EXACTLY a
//! fresh draft of the already-drafted face by `new − old` with the
//! same pull direction and neutral reference. No unfillet-style
//! restore is needed; the whole cost of the parameter tweak is one
//! planar rotation plus the neighbour recomputation — the mechanism
//! the feature-tree replay routes parameter edits through.

use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::DraftError;
use crate::types::{AngleVariation, DraftOptions, DraftResult, DraftSpec, NeutralReference};

/// Change an already-drafted planar face's draft angle in place.
///
/// `old_angle` and `new_angle` are the absolute angles (radians) of
/// the existing and the desired draft with respect to the SAME pull
/// direction and neutral reference the face was originally drafted
/// with; the modification applies the delta. The face must be planar
/// with a constant angle (the additive-composition argument is exact
/// for the planar Rodrigues rotation; variable-angle and non-planar
/// drafts re-execute through [`crate::draft_faces`] instead).
///
/// # Errors
///
/// - [`DraftError::InvalidFace`] when the face is missing or not
///   planar.
/// - Any error from the underlying [`crate::draft_faces`] delta
///   application.
pub fn modify_draft_angle(
    brep: &BRep,
    face: FaceHandle,
    old_angle: f64,
    new_angle: f64,
    pull_direction: Vec3,
    neutral: NeutralReference,
    tolerance: Tolerance,
) -> Result<DraftResult, DraftError> {
    let is_planar = brep.faces.get(face).is_some_and(|fd| {
        fd.surface.as_ref().is_some_and(|sof| {
            matches!(
                sof.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
            )
        })
    });
    if !is_planar {
        return Err(DraftError::InvalidFace(format!(
            "in-place draft-angle modification supports planar faces only \
             (face {face:?}); re-execute the draft via draft_faces for other classes"
        )));
    }
    if matches!(neutral, NeutralReference::Auto) {
        return Err(DraftError::InvalidNeutral(
            "in-place angle modification needs an anchored neutral reference (Plane or              Edge): Auto re-derives the line on the tilted face and is not guaranteed to              reproduce the original hinge, breaking the additive composition"
                .to_string(),
        ));
    }
    let delta = new_angle - old_angle;
    if delta.abs() < 1.0e-15 {
        return Ok(DraftResult {
            brep: brep.clone(),
            warnings: Vec::new(),
        });
    }
    crate::draft_faces(
        brep,
        &[DraftSpec {
            face,
            angle: delta,
            pull_direction,
            neutral,
            angle_variation: AngleVariation::Constant,
        }],
        &DraftOptions::default(),
        tolerance,
        None,
    )
}

#[cfg(test)]
#[path = "modify_tests.rs"]
mod tests;
