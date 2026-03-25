//! Face classification for shell operations.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::types::ShellDirection;
use popkern_core::topology::{BRep, FaceHandle};

/// Result of classifying faces for a shell operation.
#[derive(Debug, Clone)]
pub struct FaceClassification {
    /// Faces to be offset: (face_handle, signed_distance, inward_flag).
    pub offset_faces: Vec<(FaceHandle, f64, bool)>,
    /// Faces to be removed (open faces of the shell).
    pub removed_faces: Vec<FaceHandle>,
}

/// Classify faces for a shell operation.
///
/// Each face is either marked for removal or assigned an offset distance
/// and direction. The `face_thicknesses` map overrides `default_thickness`
/// for specific faces.
///
/// For `Inward` direction, the inner skin is offset inward (against the
/// face normal). For `Outward`, it is offset outward (along the face
/// normal). For `Midplane`, the thickness is halved and the inner skin
/// is offset inward.
///
/// The `SurfaceOnFace.orientation` flag is respected: when false (flipped
/// normal), the inward/outward sense is reversed for that face.
///
/// # Errors
///
/// Returns [`ShellError::NoFaces`] if the BRep has no faces.
pub fn classify_faces(
    brep: &BRep,
    default_thickness: f64,
    face_thicknesses: &[(FaceHandle, f64)],
    direction: ShellDirection,
    removed: &[FaceHandle],
) -> Result<FaceClassification, ShellError> {
    if brep.faces.is_empty() {
        return Err(ShellError::NoFaces);
    }

    let thickness_map: HashMap<FaceHandle, f64> = face_thicknesses.iter().copied().collect();

    let mut offset_faces = Vec::new();
    let mut removed_faces = Vec::new();

    for (fh, face) in &brep.faces {
        if removed.contains(&fh) {
            removed_faces.push(fh);
            continue;
        }

        let thickness = thickness_map.get(&fh).copied().unwrap_or(default_thickness);

        // Determine the base inward flag and distance from the direction.
        let (distance, base_inward) = match direction {
            ShellDirection::Inward => (thickness, true),
            ShellDirection::Outward => (thickness, false),
            ShellDirection::Midplane => (thickness / 2.0, true),
        };

        // If the face surface orientation is flipped, reverse the inward sense.
        let surface_flipped = face
            .surface
            .as_ref()
            .map(|sof| !sof.orientation)
            .unwrap_or(false);

        let inward = if surface_flipped {
            !base_inward
        } else {
            base_inward
        };

        offset_faces.push((fh, distance, inward));
    }

    Ok(FaceClassification {
        offset_faces,
        removed_faces,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn empty_brep_returns_no_faces() {
        let brep = BRep::new(tol());
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]);
        assert!(matches!(result, Err(ShellError::NoFaces)));
    }

    #[test]
    fn cube_all_faces_offset() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]).unwrap();
        assert_eq!(result.offset_faces.len(), 6);
        assert_eq!(result.removed_faces.len(), 0);
        // All faces should have distance 0.1 and inward=true.
        for &(_, dist, inward) in &result.offset_faces {
            assert!((dist - 0.1).abs() < 1e-12);
            assert!(inward);
        }
    }

    #[test]
    fn cube_one_face_removed() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result =
            classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[first_face]).unwrap();
        assert_eq!(result.offset_faces.len(), 5);
        assert_eq!(result.removed_faces.len(), 1);
        assert_eq!(result.removed_faces[0], first_face);
    }

    #[test]
    fn midplane_halves_thickness() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = classify_faces(&brep, 0.2, &[], ShellDirection::Midplane, &[]).unwrap();
        for &(_, dist, inward) in &result.offset_faces {
            assert!((dist - 0.1).abs() < 1e-12);
            assert!(inward);
        }
    }

    #[test]
    fn outward_direction_sets_inward_false() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Outward, &[]).unwrap();
        for &(_, _, inward) in &result.offset_faces {
            assert!(!inward);
        }
    }

    #[test]
    fn per_face_thickness_override() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = classify_faces(
            &brep,
            0.1,
            &[(first_face, 0.5)],
            ShellDirection::Inward,
            &[],
        )
        .unwrap();
        // Find the overridden face.
        let overridden = result
            .offset_faces
            .iter()
            .find(|(fh, _, _)| *fh == first_face)
            .unwrap();
        assert!((overridden.1 - 0.5).abs() < 1e-12);
        // Other faces should have default thickness.
        for &(fh, dist, _) in &result.offset_faces {
            if fh != first_face {
                assert!((dist - 0.1).abs() < 1e-12);
            }
        }
    }
}
