#![warn(missing_docs)]
//! popkern-shell: Shell and offset solid operations for the popkern BRep kernel.
//!
//! Provides the [`shell`] function for creating thin-walled solids by offsetting
//! faces of a BRep body.

pub mod classify;
pub mod edge_recompute;
pub mod error;
pub mod inner_skin;
pub mod skin_map;
pub mod surface_offset;
pub mod types;
pub mod vertex_recompute;

pub use error::{ShellError, ShellWarning};
pub use types::*;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Create a thin-walled solid by offsetting faces of a BRep body.
///
/// Removes the specified faces and offsets the remaining faces inward (or
/// outward, depending on `direction`) by `thickness`, producing a hollow
/// shell solid.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, all faces fail to
/// offset, or the offset distance is too large.
pub fn shell(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    // 1. Classify faces.
    let classification =
        classify::classify_faces(brep, thickness, face_thicknesses, direction, removed_faces)?;

    // 2. Compute offset surfaces.
    let offset_result = surface_offset::compute_offset_surfaces(
        brep,
        &classification,
        options.allow_partial,
        &tolerance,
    )?;

    // Merge auto-removed faces into the removed list.
    let mut all_removed = classification.removed_faces.clone();
    all_removed.extend_from_slice(&offset_result.auto_removed);

    // 3. Clone outer BRep into inner skin with handle mapping.
    let (mut inner, skin_map) = inner_skin::build_inner_skin(brep, &all_removed)?;

    // 4. Recompute inner skin vertex positions.
    vertex_recompute::recompute_vertices(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &tolerance,
    )?;

    // 5. Recompute inner skin edge curves.
    edge_recompute::recompute_edges(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &tolerance,
    )?;

    // 6. Replace face geometry and flip normals.
    inner_skin::apply_offset_geometry(&mut inner, &offset_result.surfaces, &skin_map)?;

    // Collect warnings.
    let mut warnings = offset_result.warnings;
    for face in &offset_result.auto_removed {
        warnings.push(ShellWarning::FaceAutoRemoved {
            face: *face,
            reason: "offset failed".into(),
        });
    }

    Ok(ShellResult {
        brep: inner,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn shell_empty_brep_returns_no_faces() {
        let brep = BRep::new(Tolerance::default());
        let result = shell(
            &brep,
            0.1,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            Tolerance::default(),
        );
        assert!(matches!(result, Err(ShellError::NoFaces)));
    }
}
