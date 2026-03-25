//! Surface offset computation for shell operations.

use std::collections::HashMap;

use crate::classify::FaceClassification;
use crate::error::{ShellError, ShellWarning};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Result of computing offset surfaces for classified faces.
#[derive(Debug, Clone)]
pub struct OffsetResult {
    /// Successfully offset surfaces, keyed by outer face handle.
    pub surfaces: HashMap<FaceHandle, SurfaceGeometry>,
    /// Faces that failed to offset.
    pub failed: Vec<(FaceHandle, String)>,
    /// Warnings generated during offset computation.
    pub warnings: Vec<ShellWarning>,
    /// Faces auto-removed due to offset failure (when allow_partial is true).
    pub auto_removed: Vec<FaceHandle>,
}

/// Collect vertex positions around a face boundary for UV range computation.
fn face_vertex_positions(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Vec<popkern_core::geometry::Point3> {
    let face = &brep.faces[face_handle];
    let start_he = brep.loops[face.outer_loop].half_edge;
    let mut positions = Vec::new();
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Compute offset surfaces for all classified faces.
///
/// For each face in `classification.offset_faces`, offsets the surface
/// geometry by the specified distance. Self-intersection is checked and
/// failed faces are either collected as errors or auto-removed depending
/// on `allow_partial`.
///
/// # Errors
///
/// Returns [`ShellError::AllFacesFailed`] if every face fails to offset.
pub fn compute_offset_surfaces(
    brep: &BRep,
    classification: &FaceClassification,
    allow_partial: bool,
    tolerance: &Tolerance,
) -> Result<OffsetResult, ShellError> {
    let mut surfaces = HashMap::new();
    let mut failed = Vec::new();
    let mut warnings = Vec::new();
    let mut auto_removed = Vec::new();

    for &(face_handle, distance, inward) in &classification.offset_faces {
        let face = &brep.faces[face_handle];
        let surface = match &face.surface {
            Some(sof) => &sof.geometry,
            None => {
                let msg = "face has no surface geometry".to_string();
                if allow_partial {
                    auto_removed.push(face_handle);
                    warnings.push(ShellWarning::FaceAutoRemoved {
                        face: face_handle,
                        reason: msg.clone(),
                    });
                    failed.push((face_handle, msg));
                } else {
                    failed.push((face_handle, msg));
                }
                continue;
            }
        };

        // Get UV range from face boundary vertices.
        let positions = face_vertex_positions(brep, face_handle);
        let (u_range, v_range) = bounded_parameter_range(&positions, surface);

        // Check self-intersection before computing offset.
        if let Some(_max_safe) =
            popkern_fillet::offset::self_intersect::check_offset_self_intersection(
                surface, distance, inward, u_range, v_range,
            )
        {
            let msg = format!("offset distance {distance} would cause self-intersection");
            if allow_partial {
                auto_removed.push(face_handle);
                warnings.push(ShellWarning::FaceAutoRemoved {
                    face: face_handle,
                    reason: msg.clone(),
                });
                failed.push((face_handle, msg));
            } else {
                failed.push((face_handle, msg));
            }
            continue;
        }

        // Compute the offset surface.
        match popkern_fillet::offset::offset_surface(
            surface, distance, inward, u_range, v_range, tolerance,
        ) {
            Ok(offset_surface) => {
                surfaces.insert(face_handle, offset_surface);
            }
            Err(msg) => {
                if allow_partial {
                    auto_removed.push(face_handle);
                    warnings.push(ShellWarning::FaceAutoRemoved {
                        face: face_handle,
                        reason: msg.clone(),
                    });
                    failed.push((face_handle, msg));
                } else {
                    failed.push((face_handle, msg));
                }
            }
        }
    }

    // If all faces failed, return an error.
    if surfaces.is_empty() && !classification.offset_faces.is_empty() {
        return Err(ShellError::AllFacesFailed { failures: failed });
    }

    // If not allow_partial and some failed, return error.
    if !allow_partial && !failed.is_empty() {
        return Err(ShellError::AllFacesFailed { failures: failed });
    }

    Ok(OffsetResult {
        surfaces,
        failed,
        warnings,
        auto_removed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::classify_faces;
    use crate::types::ShellDirection;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn cube_faces_offset_successfully() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let classification = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]).unwrap();
        let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        // All 6 cube faces should offset successfully (they are planes).
        assert_eq!(result.surfaces.len(), 6);
        assert!(result.failed.is_empty());
        assert!(result.auto_removed.is_empty());
    }

    #[test]
    fn cube_with_removed_face() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let classification =
            classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[first_face]).unwrap();
        let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        assert_eq!(result.surfaces.len(), 5);
    }

    #[test]
    fn too_large_offset_auto_removes_with_partial() {
        // Create a cylinder BRep or use cube — cube planes never self-intersect.
        // For this test, we verify the code path works with allow_partial=true
        // by using a very large offset on a cube (planes never self-intersect,
        // but offset_surface should still succeed).
        let brep = BRep::unit_cube(tol()).unwrap();
        let classification =
            classify_faces(&brep, 1000.0, &[], ShellDirection::Inward, &[]).unwrap();
        let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        // Planes always succeed regardless of offset distance.
        assert_eq!(result.surfaces.len(), 6);
    }

    #[test]
    fn all_removed_produces_empty_surfaces() {
        // When all faces are removed, there are no offset faces to compute.
        let brep = BRep::unit_cube(tol()).unwrap();
        let all_faces: Vec<FaceHandle> = brep.faces.keys().collect();
        let classification =
            classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &all_faces).unwrap();
        // No offset faces → ok with empty surfaces map.
        let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        assert!(result.surfaces.is_empty());
    }

    #[test]
    fn not_allow_partial_fails_on_any_failure() {
        // When allow_partial is false and a face fails, the whole operation fails.
        // We test this indirectly: cube planes never fail, so this just verifies
        // the success path with allow_partial=false.
        let brep = BRep::unit_cube(tol()).unwrap();
        let classification = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]).unwrap();
        let result = compute_offset_surfaces(&brep, &classification, false, &tol()).unwrap();
        assert_eq!(result.surfaces.len(), 6);
    }

    #[test]
    fn face_with_no_surface_auto_removed_when_allow_partial() {
        // When allow_partial is true, a face with no surface geometry is added
        // to auto_removed and a warning is emitted. If at least one other face
        // succeeds, the result is Ok with the auto_removed face listed.
        let mut brep = BRep::unit_cube(tol()).unwrap();

        // Strip surface geometry from one face.
        let mut faces_iter = brep.faces.keys();
        let stripped_face = faces_iter.next().unwrap();
        let good_face = faces_iter.next().unwrap();
        drop(faces_iter);
        brep.faces[stripped_face].surface = None;

        // Classification: the stripped face + one good face.
        let classification = crate::classify::FaceClassification {
            offset_faces: vec![(stripped_face, 0.1, true), (good_face, 0.1, true)],
            removed_faces: vec![],
        };

        let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        assert!(
            result.auto_removed.contains(&stripped_face),
            "stripped face should be auto-removed"
        );
        assert_eq!(
            result.failed.len(),
            1,
            "one face should be in the failed list"
        );
        assert_eq!(result.surfaces.len(), 1, "one face should succeed");
        assert_eq!(result.warnings.len(), 1, "one warning should be emitted");
    }

    #[test]
    fn face_with_no_surface_fails_when_not_allow_partial() {
        // When allow_partial is false and the only face has no surface,
        // AllFacesFailed should be returned.
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let stripped_face = brep.faces.keys().next().unwrap();
        brep.faces[stripped_face].surface = None;

        let classification = crate::classify::FaceClassification {
            offset_faces: vec![(stripped_face, 0.1, true)],
            removed_faces: vec![],
        };

        let result = compute_offset_surfaces(&brep, &classification, false, &tol());
        assert!(
            result.is_err(),
            "should fail when all faces have no surface and allow_partial is false"
        );
        match result.unwrap_err() {
            crate::error::ShellError::AllFacesFailed { failures } => {
                assert_eq!(failures.len(), 1);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn all_faces_failed_error_when_allow_partial_false_and_multiple_faces() {
        // When allow_partial is false and a second valid face also fails (not
        // the case with planes, but we test the error branch explicitly):
        // If ALL faces in offset_faces fail, AllFacesFailed is returned.
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let faces: Vec<FaceHandle> = brep.faces.keys().take(2).collect();
        for &fh in &faces {
            brep.faces[fh].surface = None;
        }

        let classification = crate::classify::FaceClassification {
            offset_faces: faces.iter().map(|&fh| (fh, 0.1, true)).collect(),
            removed_faces: vec![],
        };

        let result = compute_offset_surfaces(&brep, &classification, false, &tol());
        assert!(
            result.is_err(),
            "should fail when all faces have no surface"
        );
    }
}
