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
    // Phase 1 (parallel, PK-329): each face's offset is a pure function
    // of its own surface and boundary — dispatch over the deterministic
    // `offset_faces` list with an ordered collect, then replay the
    // decisions sequentially so warning order, failure order, and
    // first-face-wins early declines are positionally identical to the
    // sequential loop at any thread count.
    let solved = crate::parallel_dispatch::dispatch_ordered(
        &classification.offset_faces,
        false,
        |&(face_handle, distance, inward)| {
            let face = &brep.faces[face_handle];
            let surface = match &face.surface {
                Some(sof) => &sof.geometry,
                None => return None,
            };
            // Get UV range from face boundary vertices.
            let positions = face_vertex_positions(brep, face_handle);
            let (u_range, v_range) = bounded_parameter_range(&positions, surface);
            // Compute the offset surface with an honest quality record.
            // Existence declines (radius through zero, spindle torus)
            // and fold declines on the approximation path come back
            // typed from the dispatch and are mapped below.
            Some(popkern_fillet::offset::offset_surface_measured(
                surface, distance, inward, u_range, v_range, tolerance,
            ))
        },
    );

    // Phase 2 (sequential): fold the per-face outcomes in list order.
    let mut surfaces = HashMap::new();
    let mut failed = Vec::new();
    let mut warnings = Vec::new();
    let mut auto_removed = Vec::new();

    for (&(face_handle, distance, _inward), outcome) in
        classification.offset_faces.iter().zip(solved)
    {
        let Some(outcome) = outcome else {
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
        };
        match outcome {
            Ok(outcome) => {
                if let Some(achieved) = outcome.approximation_deviation {
                    if achieved > tolerance.linear() {
                        // The finest sampling grid still missed the stage
                        // tolerance. Under allow_partial the face follows
                        // the existing failed-face contract (auto-removed
                        // with a warning); otherwise the operation
                        // declines with the measured gap.
                        if allow_partial {
                            let msg = format!(
                                "offset approximation deviation {achieved:.3e} exceeds \
                                 tolerance {:.3e}",
                                tolerance.linear()
                            );
                            auto_removed.push(face_handle);
                            warnings.push(ShellWarning::FaceAutoRemoved {
                                face: face_handle,
                                reason: msg.clone(),
                            });
                            failed.push((face_handle, msg));
                        } else {
                            return Err(ShellError::ApproximationOutOfTolerance {
                                face: face_handle,
                                achieved,
                                requested: tolerance.linear(),
                            });
                        }
                        continue;
                    }
                    warnings.push(ShellWarning::FaceOffsetApproximated {
                        face: face_handle,
                        achieved_deviation: achieved,
                    });
                }
                surfaces.insert(face_handle, outcome.surface);
            }
            Err(popkern_fillet::offset::OffsetError::SpindleTorus { new_minor, .. }) => {
                // Spindled torus (tube swallows the axis): the offset
                // SURFACE exists as an immersed spindle; construct it
                // and let the trim stage's pinch arm serve the apple
                // (PK-327 §6). The surface-level bound decline was the
                // pre-trim containment.
                let Some(SurfaceGeometry::Torus {
                    center,
                    axis,
                    major_radius,
                    ..
                }) = brep.faces[face_handle]
                    .surface
                    .as_ref()
                    .map(|s| &s.geometry)
                else {
                    return Err(ShellError::Topology(
                        "spindle decline on a non-torus surface".into(),
                    ));
                };
                surfaces.insert(
                    face_handle,
                    SurfaceGeometry::Torus {
                        center: *center,
                        axis: *axis,
                        major_radius: *major_radius,
                        minor_radius: new_minor,
                    },
                );
            }
            Err(e) => {
                // Existence/fold declines are properties of the REQUEST
                // (surface class × distance × direction), not of one
                // face's health — auto-removing the face would punch an
                // opening the user never asked for, so the operation
                // declines outright with the safe bound and the binding
                // face named.
                if let Some(max_safe) = e.max_safe() {
                    return Err(ShellError::OffsetTooLarge {
                        requested: distance,
                        max_safe,
                        binding: format!("face {face_handle:?}: {e}"),
                    });
                }
                let msg = e.to_string();
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
    fn all_removed_declines_at_classification() {
        // Removing EVERY face leaves nothing to offset: classification
        // declines up front (PK-347 — the empty result used to sail
        // through the pipeline and serve a 0-face body clean).
        let brep = BRep::unit_cube(tol()).unwrap();
        let all_faces: Vec<FaceHandle> = brep.faces.keys().collect();
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &all_faces);
        assert!(matches!(result, Err(crate::error::ShellError::NoFaces)));
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
        let (stripped_face, good_face) = {
            let mut face_keys = brep.faces.keys();
            (face_keys.next().unwrap(), face_keys.next().unwrap())
        };
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
