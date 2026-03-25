//! Offset solid operations.
//!
//! Offsets all faces of a solid by a uniform or per-face distance,
//! producing a new solid that is larger (positive offset) or smaller
//! (negative offset) than the original.

use crate::ShellWarning;
use crate::classify::classify_faces;
use crate::edge_recompute::recompute_edges;
use crate::error::ShellError;
use crate::inner_skin::{apply_offset_geometry, build_inner_skin};
use crate::surface_offset::compute_offset_surfaces;
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use crate::vertex_recompute::recompute_vertices;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Offset all faces of a solid by a uniform distance.
///
/// A positive `distance` expands the solid outward; a negative distance
/// shrinks it inward.  No faces are removed and no side walls are built.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces or all faces fail.
pub fn offset_solid(
    brep: &BRep,
    distance: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    // ShellDirection controls offset relative to the face normal.
    // build_inner_skin() clones the BRep with flipped normals, so after
    // that flip: Inward (push against the *original* normal) actually
    // moves faces outward in world space, and Outward does the reverse.
    // Therefore positive distance (grow outward) maps to Inward, and
    // negative distance (shrink inward) maps to Outward.
    let direction = if distance >= 0.0 {
        ShellDirection::Inward
    } else {
        ShellDirection::Outward
    };
    let abs_dist = distance.abs();

    let classification = classify_faces(brep, abs_dist, &[], direction, &[])?;
    let stage_tol = Tolerance::default().with_linear(tolerance.linear / 4.0);
    let offset_result =
        compute_offset_surfaces(brep, &classification, options.allow_partial, &stage_tol)?;

    let (mut result, skin_map) = build_inner_skin(brep, &[])?;

    recompute_vertices(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    recompute_edges(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    apply_offset_geometry(&mut result, &offset_result.surfaces, &skin_map)?;

    // Restore face orientations (offset_solid doesn't flip normals).
    for (_, face) in &mut result.faces {
        if let Some(sof) = &mut face.surface {
            sof.orientation = !sof.orientation;
        }
    }

    let mut warnings = offset_result.warnings;
    for face in &offset_result.auto_removed {
        warnings.push(ShellWarning::FaceAutoRemoved {
            face: *face,
            reason: "offset failed".into(),
        });
    }

    Ok(ShellResult {
        brep: result,
        warnings,
    })
}

/// Offset all faces with per-face variable distances.
///
/// Faces not listed in `face_distances` use `default_distance`.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces or all faces fail.
pub fn offset_solid_variable(
    brep: &BRep,
    default_distance: f64,
    face_distances: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    let direction = if default_distance >= 0.0 {
        ShellDirection::Inward
    } else {
        ShellDirection::Outward
    };
    let abs_default = default_distance.abs();
    let abs_faces: Vec<(FaceHandle, f64)> = face_distances
        .iter()
        .map(|&(fh, d)| (fh, d.abs()))
        .collect();

    let classification = classify_faces(brep, abs_default, &abs_faces, direction, &[])?;
    let stage_tol = Tolerance::default().with_linear(tolerance.linear / 4.0);
    let offset_result =
        compute_offset_surfaces(brep, &classification, options.allow_partial, &stage_tol)?;

    let (mut result, skin_map) = build_inner_skin(brep, &[])?;

    recompute_vertices(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    recompute_edges(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    apply_offset_geometry(&mut result, &offset_result.surfaces, &skin_map)?;

    // Restore face orientations.
    for (_, face) in &mut result.faces {
        if let Some(sof) = &mut face.surface {
            sof.orientation = !sof.orientation;
        }
    }

    let mut warnings = offset_result.warnings;
    for face in &offset_result.auto_removed {
        warnings.push(ShellWarning::FaceAutoRemoved {
            face: *face,
            reason: "offset failed".into(),
        });
    }

    Ok(ShellResult {
        brep: result,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn offset_cube_outward() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
        assert_eq!(result.brep.faces.len(), 6);
        // Positive distance grows the solid outward: corners move from
        // {0,1} to {-0.1, 1.1}.
        for (_, v) in &result.brep.vertices {
            let coords = [v.position.x, v.position.y, v.position.z];
            for &c in &coords {
                assert!(
                    (c - (-0.1)).abs() < 0.05 || (c - 1.1).abs() < 0.05,
                    "vertex coord {c} unexpected after outward offset"
                );
            }
        }
    }

    #[test]
    fn offset_cube_negative() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = offset_solid(&brep, -0.1, &ShellOptions::default(), tol()).unwrap();
        assert_eq!(result.brep.faces.len(), 6);
        // Negative distance shrinks the solid inward: corners move from
        // {0,1} to {0.1, 0.9}.
        for (_, v) in &result.brep.vertices {
            let coords = [v.position.x, v.position.y, v.position.z];
            for &c in &coords {
                assert!(
                    (c - 0.1).abs() < 0.05 || (c - 0.9).abs() < 0.05,
                    "vertex coord {c} unexpected after inward offset"
                );
            }
        }
    }

    #[test]
    fn offset_empty_brep_errors() {
        let brep = BRep::new(tol());
        let result = offset_solid(&brep, 0.1, &ShellOptions::default(), tol());
        assert!(matches!(result, Err(ShellError::NoFaces)));
    }

    #[test]
    fn offset_variable_produces_result() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = offset_solid_variable(
            &brep,
            0.1,
            &[(first_face, 0.2)],
            &ShellOptions::default(),
            tol(),
        )
        .unwrap();
        assert_eq!(result.brep.faces.len(), 6);
    }
}
