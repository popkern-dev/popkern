//! Offset solid operations.
//!
//! Offsets all faces of a solid by a uniform or per-face distance,
//! producing a new solid that is larger (positive offset) or smaller
//! (negative offset) than the original.

use crate::ShellWarning;
use crate::classify::classify_faces;
use crate::error::ShellError;
use crate::inner_skin::{apply_offset_geometry, build_inner_skin};
use crate::surface_offset::compute_offset_surfaces;
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Honest per-run stage count of [`offset_impl`] for progress totals
/// (PK-276).
pub(crate) const OFFSET_STAGE_COUNT: usize = 7;

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
    offset_solid_impl(
        brep,
        distance,
        &[],
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )
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
    offset_solid_impl(
        brep,
        default_distance,
        face_distances,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// Shared offset ladder for [`offset_solid`] and
/// [`offset_solid_variable`] (identical stage sequence; the variable
/// entry adds per-face distances).
///
/// `tracker` threads PK-276 interactive control: a cancellation poll
/// before every stage (`ShellError::Cancelled`, model untouched) and
/// one progress event per stage.
pub(crate) fn offset_solid_impl(
    brep: &BRep,
    default_distance: f64,
    face_distances: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<ShellResult, ShellError> {
    use std::time::Instant;
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(ShellError::Cancelled);
            }
        };
    }
    // Directions are MATERIAL directions (classify resolves the
    // component material sense): positive distance grows the solid —
    // every face moves out of the material — and negative shrinks it.
    // (The historical mapping compensated for an inverted flag chain
    // that moved "Inward" skins outward; both lies died together.)
    let direction = if default_distance >= 0.0 {
        ShellDirection::Outward
    } else {
        ShellDirection::Inward
    };
    let abs_default = default_distance.abs();
    let abs_faces: Vec<(FaceHandle, f64)> = face_distances
        .iter()
        .map(|&(fh, d)| (fh, d.abs()))
        .collect();

    stage_gate!();
    let t0 = Instant::now();
    let classification = classify_faces(brep, abs_default, &abs_faces, direction, &[])?;
    tracker.finish_stage("offset.classify_faces", t0, brep.faces.len());

    stage_gate!();
    let t0 = Instant::now();
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);
    let offset_result =
        compute_offset_surfaces(brep, &classification, options.allow_partial, &stage_tol)?;
    tracker.finish_stage("offset.offset_surfaces", t0, offset_result.surfaces.len());

    stage_gate!();
    let t0 = Instant::now();
    let (mut result, skin_map) = build_inner_skin(brep, &[])?;
    tracker.finish_stage("offset.skin", t0, result.faces.len());

    stage_gate!();
    let t0 = Instant::now();
    let join_plan = crate::join_classification::classify_joins(
        &result,
        &skin_map,
        &offset_result.surfaces,
        &stage_tol,
        options.join_style,
    )?;
    crate::vertex_recompute::recompute_vertices_with_gap_corners(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
        &join_plan.gap_corners,
    )?;
    tracker.finish_stage("offset.vertex_recompute", t0, result.vertices.len());

    stage_gate!();
    let t0 = Instant::now();
    let undone = crate::edge_recompute::recompute_edges_with_skips(
        &mut result,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
        &join_plan.arc_skip_edges(),
    )?;
    crate::join_classification::intersection_gap_postcheck(&join_plan, &undone)?;
    // Collapse gate (part of the recompute stage): a shrink past the
    // local feature size inverts the skin.
    crate::collapse_check::check_skin_collapse(
        brep,
        &result,
        &skin_map,
        default_distance,
        &tolerance,
    )?;
    tracker.finish_stage("offset.edge_recompute", t0, result.edges.len());

    stage_gate!();
    let t0 = Instant::now();
    apply_offset_geometry(&mut result, &offset_result.surfaces, &skin_map)?;
    // Join construction (PK-325): rolling-ball strips + sphere corner
    // windows, via Euler operators only.
    crate::arc_join::construct_arc_joins(
        &mut result,
        &join_plan,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    tracker.finish_stage("offset.apply_offset_geometry", t0, result.faces.len());

    // The clone carries the ORIGINAL pcurves. Most exact offsets are
    // parametrization-preserving (planes, cylinders, spheres, tori) so
    // those stay valid, but the exact cone arm SHIFTS the axial
    // parameter (apex translation) and the inherited pcurves/trim
    // windows go geometrically stale. The recompute's survival rule
    // keeps every still-accurate pcurve and re-authors only the stale
    // ones.
    stage_gate!();
    let t0 = Instant::now();
    result
        .recompute_all_pcurves_winding_repaired()
        .map_err(|e| ShellError::Topology(e.to_string()))?;
    tracker.finish_stage("offset.pcurve_recompute", t0, result.faces.len());

    // Orientation flags are DERIVED state (a function of finalized
    // loop winding + surface geometry): reconcile them at the serving
    // boundary. Join faces on charts whose handedness opposes their
    // construction-time material measurement (the cone-rim torus
    // strips) get their honest value here (PK-325).
    result.reconcile_orientation_flags();

    // Serve gate: never hand back a body a validator would reject.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    let mut warnings = offset_result.warnings;
    warnings.extend(join_plan.warnings.iter().cloned());
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

    // ---------------------------------------------------------------
    // offset_solid — existing tests
    // ---------------------------------------------------------------

    #[test]
    fn offset_cube_outward() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
        // Minkowski-faithful default (PK-325): the rounded box.
        assert_eq!(result.brep.faces.len(), 26);
        // Every corner-sphere vertex sits at an original corner
        // displaced by 0.1 along one face normal: coordinates come
        // from {-0.1, 0, 1, 1.1}.
        for (_, v) in &result.brep.vertices {
            let coords = [v.position.x, v.position.y, v.position.z];
            for &c in &coords {
                assert!(
                    [-0.1, 0.0, 1.0, 1.1]
                        .iter()
                        .any(|expect| (c - expect).abs() < 1e-9),
                    "vertex coord {c} unexpected after outward arc offset"
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

    // ---------------------------------------------------------------
    // offset_solid_variable — dedicated tests
    // ---------------------------------------------------------------

    #[test]
    fn variable_empty_brep_errors() {
        let brep = BRep::new(tol());
        let result = offset_solid_variable(&brep, 0.1, &[], &ShellOptions::default(), tol());
        assert!(
            matches!(result, Err(ShellError::NoFaces)),
            "variable offset on empty BRep should return NoFaces"
        );
    }

    #[test]
    fn variable_with_no_overrides_matches_uniform() {
        // offset_solid_variable with empty face_distances should behave
        // identically to offset_solid with the same default distance.
        let brep = BRep::unit_cube(tol()).unwrap();
        let uniform = offset_solid(&brep, 0.15, &ShellOptions::default(), tol()).unwrap();
        let variable =
            offset_solid_variable(&brep, 0.15, &[], &ShellOptions::default(), tol()).unwrap();
        assert_eq!(uniform.brep.faces.len(), variable.brep.faces.len());
        assert_eq!(uniform.brep.vertices.len(), variable.brep.vertices.len());
    }

    #[test]
    fn variable_negative_default_shrinks_inward() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result =
            offset_solid_variable(&brep, -0.1, &[], &ShellOptions::default(), tol()).unwrap();
        assert_eq!(result.brep.faces.len(), 6);
        // Negative default should shrink: corners at {0.1, 0.9}.
        for (_, v) in &result.brep.vertices {
            let coords = [v.position.x, v.position.y, v.position.z];
            for &c in &coords {
                assert!(
                    (c - 0.1).abs() < 0.05 || (c - 0.9).abs() < 0.05,
                    "vertex coord {c} unexpected after negative variable offset"
                );
            }
        }
    }

    #[test]
    fn variable_per_face_override_differs_from_default() {
        // When one face gets a larger offset, the result's vertex positions
        // should differ from a uniform offset.
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();

        let uniform = offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
        let variable = offset_solid_variable(
            &brep,
            0.1,
            &[(first_face, 0.3)],
            &ShellOptions::default(),
            tol(),
        )
        .unwrap();

        // The uniform run serves the rounded box (arc default); the
        // variable run has no uniform rolling ball and contains to the
        // sharp join (E4 downgrade) with its 6 faces preserved.
        assert_eq!(uniform.brep.faces.len(), 26);
        assert_eq!(variable.brep.faces.len(), 6);

        // At least some vertices should differ between uniform and variable
        // offset results. Collect all vertex positions as sorted tuples for
        // comparison.
        let positions = |brep: &BRep| -> Vec<[i64; 3]> {
            brep.vertices
                .iter()
                .map(|(_, v)| {
                    [
                        (v.position.x * 1000.0).round() as i64,
                        (v.position.y * 1000.0).round() as i64,
                        (v.position.z * 1000.0).round() as i64,
                    ]
                })
                .collect()
        };
        let uni_pos = positions(&uniform.brep);
        let var_pos = positions(&variable.brep);
        assert_ne!(
            uni_pos, var_pos,
            "variable offset with a per-face override should produce \
             different vertex positions than uniform offset"
        );
    }

    #[test]
    fn variable_multiple_overrides() {
        // Provide per-face overrides for two faces.
        let brep = BRep::unit_cube(tol()).unwrap();
        let faces: Vec<FaceHandle> = brep.faces.keys().collect();
        let overrides = vec![(faces[0], 0.2), (faces[1], 0.05)];
        let result =
            offset_solid_variable(&brep, 0.1, &overrides, &ShellOptions::default(), tol()).unwrap();
        assert_eq!(result.brep.faces.len(), 6);
    }

    #[test]
    fn variable_zero_default_with_overrides() {
        // Default distance is 0 — only faces with explicit overrides move.
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = offset_solid_variable(
            &brep,
            0.0,
            &[(first_face, 0.15)],
            &ShellOptions::default(),
            tol(),
        )
        .unwrap();
        assert_eq!(result.brep.faces.len(), 6);
    }

    #[test]
    fn variable_preserves_face_count() {
        // No matter what per-face distances we choose (within reason),
        // the face count should stay at 6 for a cube.
        let brep = BRep::unit_cube(tol()).unwrap();
        let faces: Vec<FaceHandle> = brep.faces.keys().collect();
        let overrides: Vec<(FaceHandle, f64)> = faces
            .iter()
            .enumerate()
            .map(|(i, &fh)| (fh, 0.05 + 0.02 * i as f64))
            .collect();
        let result =
            offset_solid_variable(&brep, 0.1, &overrides, &ShellOptions::default(), tol()).unwrap();
        assert_eq!(result.brep.faces.len(), 6);
    }

    #[test]
    fn variable_negative_override_with_positive_default() {
        // face_distances uses .abs(), so a negative override is treated
        // as its absolute value. Verify it doesn't error.
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = offset_solid_variable(
            &brep,
            0.1,
            &[(first_face, -0.2)],
            &ShellOptions::default(),
            tol(),
        )
        .unwrap();
        assert_eq!(result.brep.faces.len(), 6);
    }

    #[test]
    fn variable_allow_partial_false_still_succeeds_on_cube() {
        // A cube with small offsets should succeed even with
        // allow_partial = false.
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let opts = ShellOptions {
            allow_partial: false,
            ..ShellOptions::default()
        };
        let result = offset_solid_variable(&brep, 0.1, &[(first_face, 0.2)], &opts, tol()).unwrap();
        assert_eq!(result.brep.faces.len(), 6);
    }

    #[test]
    fn variable_small_offset_vertices_stay_close() {
        // A very small offset should produce vertices close to the original.
        let brep = BRep::unit_cube(tol()).unwrap();
        let result =
            offset_solid_variable(&brep, 0.01, &[], &ShellOptions::default(), tol()).unwrap();
        for (_, v) in &result.brep.vertices {
            let coords = [v.position.x, v.position.y, v.position.z];
            for &c in &coords {
                assert!(
                    (c - (-0.01)).abs() < 0.05
                        || (c - 0.0).abs() < 0.05
                        || (c - 1.0).abs() < 0.05
                        || (c - 1.01).abs() < 0.05,
                    "vertex coord {c} too far from expected range \
                     after small outward offset"
                );
            }
        }
    }

    #[test]
    fn variable_result_has_no_warnings_on_clean_cube() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result =
            offset_solid_variable(&brep, 0.1, &[], &ShellOptions::default(), tol()).unwrap();
        assert!(
            result.warnings.is_empty(),
            "a clean cube offset should produce no warnings, got {} warnings",
            result.warnings.len()
        );
    }
}
