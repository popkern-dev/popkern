#![warn(missing_docs)]
//! popkern-draft: Draft/taper operations for the popkern BRep kernel.
//!
//! Provides [`draft_faces`] and [`draft_face`] for applying draft angles
//! to BRep faces, used in mold design to ensure parts can be withdrawn
//! from the mold. Also provides [`draft_parting_line`] for split-direction
//! drafts.

pub mod adjacent;
pub mod error;
pub mod neutral;
pub mod parting;
pub mod recompute;
pub mod surface;
pub mod topology;
pub mod types;

pub use error::{DraftError, DraftWarning};
pub use parting::draft_parting_line;
pub use types::*;

use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Apply draft angles to the specified faces of a BRep solid.
///
/// Each [`DraftSpec`] identifies a face, a draft angle, and a pull direction.
/// The face surface is rotated about the neutral line by the draft angle,
/// and adjacent edges and vertices are recomputed.
///
/// # Errors
///
/// Returns [`DraftError`] if no specs are provided, a face is invalid,
/// or the draft angle exceeds the maximum.
pub fn draft_faces(
    brep: &BRep,
    specs: &[DraftSpec],
    options: &DraftOptions,
    tolerance: Tolerance,
) -> Result<DraftResult, DraftError> {
    if specs.is_empty() {
        return Err(DraftError::EmptySpecs);
    }

    // Validate angles and faces.
    for spec in specs {
        if spec.angle.abs() > options.max_angle {
            return Err(DraftError::AngleTooLarge {
                angle: spec.angle,
                max: options.max_angle,
            });
        }
        if !brep.faces.contains_key(spec.face) {
            return Err(DraftError::InvalidFace(format!(
                "face {:?} not in BRep",
                spec.face
            )));
        }
    }

    // Stage 1: Compute neutral lines for each face.
    let mut neutral_lines = Vec::new();
    for spec in specs {
        let nl = neutral::compute_neutral_line(
            brep,
            spec.face,
            &spec.neutral,
            &spec.pull_direction,
            &tolerance,
        )?;
        neutral_lines.push(nl);
    }

    // Stage 2: Compute drafted surfaces.
    let mut drafted_surfaces = std::collections::HashMap::new();
    for (i, spec) in specs.iter().enumerate() {
        let new_surface = surface::compute_drafted_surface(
            brep,
            spec.face,
            &neutral_lines[i],
            spec.angle,
            &spec.angle_variation,
            &spec.pull_direction,
            options,
            &tolerance,
        )?;
        drafted_surfaces.insert(spec.face, new_surface);
    }

    // Stages 3-4: Apply topology modifications (clone, surface replacement,
    // vertex/edge recomputation, adjacent face handling).
    let (result_brep, warnings) =
        topology::apply_draft_topology(brep, &drafted_surfaces, options, &tolerance)?;

    Ok(DraftResult {
        brep: result_brep,
        warnings,
    })
}

/// Convenience wrapper to draft a single face.
///
/// Equivalent to calling [`draft_faces`] with one [`DraftSpec`] and
/// default [`DraftOptions`].
///
/// # Errors
///
/// Returns [`DraftError`] if the face is invalid or the angle is too large.
pub fn draft_face(
    brep: &BRep,
    face: FaceHandle,
    angle: f64,
    pull_direction: Vec3,
    neutral: NeutralReference,
    tolerance: Tolerance,
) -> Result<DraftResult, DraftError> {
    draft_faces(
        brep,
        &[DraftSpec {
            face,
            angle,
            pull_direction,
            neutral,
            angle_variation: AngleVariation::Constant,
        }],
        &DraftOptions::default(),
        tolerance,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn draft_empty_specs_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = draft_faces(&brep, &[], &DraftOptions::default(), tol());
        assert!(matches!(result, Err(DraftError::EmptySpecs)));
    }

    #[test]
    fn draft_angle_too_large_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let specs = vec![DraftSpec {
            face,
            angle: 1.0, // > pi/4
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        }];
        let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol());
        assert!(matches!(result, Err(DraftError::AngleTooLarge { .. })));
    }

    #[test]
    fn draft_invalid_face_returns_error() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let specs = vec![DraftSpec {
            face: Default::default(),
            angle: 0.1,
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        }];
        let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol());
        assert!(matches!(result, Err(DraftError::InvalidFace(_))));
    }

    #[test]
    fn draft_cube_face_produces_result() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let specs = vec![DraftSpec {
            face,
            angle: 0.087, // ~5 degrees
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        }];
        let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol());
        assert!(result.is_ok(), "draft failed: {:?}", result.err());
        assert_eq!(result.unwrap().brep.faces.len(), 6);
    }

    #[test]
    fn draft_face_convenience_matches_draft_faces() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let pull = Vec3::new(0.0, 0.0, 1.0);
        let angle = 0.087;

        let result1 = draft_face(&brep, face, angle, pull, NeutralReference::Auto, tol());
        let result2 = draft_faces(
            &brep,
            &[DraftSpec {
                face,
                angle,
                pull_direction: pull,
                neutral: NeutralReference::Auto,
                angle_variation: AngleVariation::Constant,
            }],
            &DraftOptions::default(),
            tol(),
        );

        assert!(result1.is_ok());
        assert!(result2.is_ok());
        let r1 = result1.unwrap();
        let r2 = result2.unwrap();
        assert_eq!(r1.brep.faces.len(), r2.brep.faces.len());
        assert_eq!(r1.brep.vertices.len(), r2.brep.vertices.len());
    }
}
