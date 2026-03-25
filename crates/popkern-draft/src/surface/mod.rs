//! Surface drafting algorithms: dispatch by surface type.

pub mod cylindrical;
pub mod general;
pub mod planar;
pub mod skinning;

use crate::error::DraftError;
use crate::neutral::NeutralLine;
use crate::types::{AngleVariation, DraftOptions};
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Compute the drafted surface for a face.
///
/// Dispatches to planar, cylindrical, or general NURBS drafting
/// based on the face's surface type.
///
/// # Errors
///
/// Returns [`DraftError`] if the surface cannot be drafted (e.g., angle
/// too large for the geometry, degenerate surface).
#[allow(clippy::too_many_arguments)]
pub fn compute_drafted_surface(
    brep: &BRep,
    face: FaceHandle,
    neutral: &NeutralLine,
    angle: f64,
    variation: &AngleVariation,
    pull_direction: &Vec3,
    options: &DraftOptions,
    tolerance: &Tolerance,
) -> Result<SurfaceGeometry, DraftError> {
    let face_data = &brep.faces[face];
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| DraftError::InvalidFace("face has no surface geometry".into()))?;

    match &surface.geometry {
        SurfaceGeometry::Plane { .. } => planar::draft_planar(&surface.geometry, neutral, angle),
        SurfaceGeometry::Cylinder { .. } => {
            cylindrical::draft_cylindrical(&surface.geometry, neutral, angle)
        }
        _ => general::draft_general(
            &surface.geometry,
            neutral,
            angle,
            variation,
            pull_direction,
            options,
            tolerance,
        ),
    }
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
    fn dispatch_planar_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        // Build a neutral line manually.
        let nl = crate::neutral::compute_neutral_line(
            &brep,
            face,
            &crate::types::NeutralReference::Auto,
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        if let Ok(neutral) = nl {
            let result = compute_drafted_surface(
                &brep,
                face,
                &neutral,
                0.087,
                &AngleVariation::Constant,
                &Vec3::new(0.0, 0.0, 1.0),
                &DraftOptions::default(),
                &tol(),
            );
            assert!(result.is_ok(), "draft surface failed: {:?}", result.err());
        }
    }
}
