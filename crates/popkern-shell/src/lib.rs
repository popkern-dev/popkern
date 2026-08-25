#![warn(missing_docs)]
//! popkern-shell: Shell and offset solid operations for the popkern BRep kernel.
//!
//! Provides the [`shell`] and [`shell_variable`] functions for creating
//! thin-walled solids, plus [`offset_solid`] and [`thicken`] operations.

pub mod assemble;
pub mod classify;
pub mod collapse_check;
mod curved_recompute;
mod curved_recompute_similarity;
pub mod edge_recompute;
pub mod error;
pub mod inner_skin;
pub mod interactive;
mod midplane;
pub mod offset;
pub mod side_walls;
mod skin_ladder;
pub mod skin_map;
pub mod surface_offset;
pub mod thicken;
pub mod types;
pub mod vertex_recompute;

pub use error::{ShellError, ShellWarning};
pub use interactive::{
    offset_solid_variable_with_control, offset_solid_with_control, shell_variable_with_control,
    shell_with_control, thicken_with_control,
};
pub use types::*;

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Collected face data for ancestry matching during shell naming.
pub(crate) struct InputFaceData {
    /// The entity name (or placeholder).
    name: EntityName,
    /// Boundary centroid.
    centroid: Point3,
    /// Surface normal at the centroid.
    normal: Vec3,
}

/// Create a thin-walled solid by offsetting faces of a BRep body.
///
/// Removes the specified faces and offsets the remaining faces inward (or
/// outward, depending on `direction`) by `thickness`, producing a hollow
/// shell solid.  When faces are removed, side wall faces are built to close
/// the gap between the outer and inner skins.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, all faces fail to
/// offset, or the offset distance is too large.
#[allow(clippy::too_many_arguments)]
pub fn shell(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<ShellResult, ShellError> {
    let input_faces = if tracking.is_some() {
        collect_shell_face_data(brep)
    } else {
        Vec::new()
    };
    let mut result = shell_impl(
        brep,
        thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, &input_faces, &brep.tolerance);
    }
    Ok(result)
}

/// Create a shell with per-face variable thickness.
///
/// Like [`shell`], but the `face_thicknesses` parameter provides
/// per-face thickness overrides. Faces not listed use `default_thickness`.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces, all faces fail to
/// offset, or the offset distance is too large.
#[allow(clippy::too_many_arguments)]
pub fn shell_variable(
    brep: &BRep,
    default_thickness: f64,
    face_thicknesses: &[(FaceHandle, f64)],
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    options: &ShellOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<ShellResult, ShellError> {
    let input_faces = if tracking.is_some() {
        collect_shell_face_data(brep)
    } else {
        Vec::new()
    };
    let mut result = shell_impl(
        brep,
        default_thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
        &mut popkern_core::progress::StageTracker::inert(),
    )?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, &input_faces, &brep.tolerance);
    }
    Ok(result)
}

/// Honest per-run stage count of [`shell_impl`] for progress totals
/// (PK-276).
pub(crate) const SHELL_STAGE_COUNT: usize = 8;

/// Shared implementation for [`shell`] and [`shell_variable`].
///
/// `tracker` threads PK-276 interactive control: a cancellation poll
/// before every stage (`ShellError::Cancelled`, model untouched — the
/// working skins are private) and one progress event per stage.
#[allow(clippy::too_many_arguments)]
pub(crate) fn shell_impl(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
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
    // Midplane moves BOTH skins; it runs its own two-ladder pipeline
    // (PK-330) instead of the single-skin ladder below.
    if matches!(direction, ShellDirection::Midplane) {
        return midplane::midplane_shell_impl(
            brep,
            thickness,
            removed_faces,
            face_thicknesses,
            options,
            tolerance,
            tracker,
        );
    }

    // Stage 1: classify faces.
    stage_gate!();
    let t0 = Instant::now();
    let classification =
        classify::classify_faces(brep, thickness, face_thicknesses, direction, removed_faces)?;
    tracker.finish_stage("shell.classify_faces", t0, brep.faces.len());

    // Stage 2: compute offset surfaces.
    stage_gate!();
    let t0 = Instant::now();
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);
    let offset_result = surface_offset::compute_offset_surfaces(
        brep,
        &classification,
        options.allow_partial,
        &stage_tol,
    )?;
    tracker.finish_stage("shell.offset_surfaces", t0, offset_result.surfaces.len());

    // Merge auto-removed faces into the removed list.
    let mut all_removed = classification.removed_faces.clone();
    all_removed.extend_from_slice(&offset_result.auto_removed);

    // Stage 3: build inner skin.
    stage_gate!();
    let t0 = Instant::now();
    let (mut inner, skin_map) = inner_skin::build_inner_skin(brep, &all_removed)?;
    tracker.finish_stage("shell.inner_skin", t0, inner.faces.len());

    // Stage 4: recompute inner skin geometry.
    stage_gate!();
    let t0 = Instant::now();
    vertex_recompute::recompute_vertices(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    tracker.finish_stage("shell.vertex_recompute", t0, inner.vertices.len());
    stage_gate!();
    let t0 = Instant::now();
    edge_recompute::recompute_edges(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    // Collapse gate (part of the recompute stage, O(E) and sub-stage
    // cheap): a thickness past the local feature size inverts skin
    // edges; the fold would serve structurally-clean wrong geometry.
    collapse_check::check_skin_collapse(brep, &inner, &skin_map, thickness, &tolerance)?;
    tracker.finish_stage("shell.edge_recompute", t0, inner.edges.len());
    stage_gate!();
    let t0 = Instant::now();
    inner_skin::apply_offset_geometry(&mut inner, &offset_result.surfaces, &skin_map)?;
    tracker.finish_stage("shell.apply_offset_geometry", t0, inner.faces.len());

    // Mirror the inner skin's material side honestly: winding and
    // orientation flags flip together, so the void skin keeps the
    // interior-left invariant (a flag-only flip ships the PK-244
    // winding-dishonest mirror class).
    inner.reverse_body_orientation();

    // Stage 5: assemble outer + inner + side walls.
    stage_gate!();
    let t0 = Instant::now();
    let mut result = assemble::assemble_shell(brep, &inner, &all_removed, &skin_map, tolerance)?;
    tracker.finish_stage("shell.assemble", t0, result.faces.len());

    // Stage 6: recompute pcurves for all faces.
    stage_gate!();
    let t0 = Instant::now();
    result
        .recompute_all_pcurves_winding_repaired()
        .map_err(|e| ShellError::Topology(e.to_string()))?;
    tracker.finish_stage("shell.pcurve_recompute", t0, result.faces.len());

    // Serve gate: never hand back a body a validator would reject. The
    // structural collection is the affordable in-pipeline form (PK-273);
    // tests and acceptance batteries additionally hold results to
    // validate_all.
    let validation_errors = result.validate_collect_structural();
    if !validation_errors.is_empty() {
        return Err(ShellError::ValidationFailed {
            errors: validation_errors,
        });
    }

    // Collect warnings.
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

/// Collect face names, centroids, and normals from a BRep for ancestry matching.
pub(crate) fn collect_shell_face_data(brep: &BRep) -> Vec<InputFaceData> {
    use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
    use popkern_core::naming_criteria::capture_face;

    let mut result = Vec::new();
    for (i, (fh, _face)) in brep.faces.iter().enumerate() {
        let Some(fc) = capture_face(brep, fh) else {
            continue;
        };
        let name = if let Some(name_id) = brep.faces[fh].name {
            if let Some(n) = brep.name_pool.get(name_id) {
                n.clone()
            } else {
                EntityName {
                    operation: OperationId(0),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: i as u32,
                    },
                    ancestors: Vec::new(),
                }
            }
        } else {
            EntityName {
                operation: OperationId(0),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            }
        };
        result.push(InputFaceData {
            name,
            centroid: fc.centroid,
            normal: fc.normal,
        });
    }
    result
}

/// Assign stable names to entities produced by a shell operation.
///
/// Matches output faces to input faces by centroid proximity and normal
/// alignment. Outer skin faces that match an input face (same-direction
/// normal) get `Modified` with the input face as ancestor. Inner skin
/// faces that match an input face (opposite normal) also get `Modified`
/// with the input face as ancestor. Side wall faces (no match) get
/// `Generated`.
pub(crate) fn assign_shell_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    input_faces: &[InputFaceData],
    tol: &Tolerance,
) {
    use popkern_core::naming::{EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let post_faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;
    let generous_dist = tol.linear() * 100.0;

    for fh in &post_faces {
        let Some(fc) = capture_face(brep, *fh) else {
            continue;
        };
        let centroid = fc.centroid;
        let normal = fc.normal;

        // Match against input faces: allow same-direction (outer) or
        // opposite-direction (inner/offset) normals.
        let best_match = input_faces
            .iter()
            .filter(|inp| {
                centroid.distance_to(&inp.centroid) < generous_dist
                    && normal.dot(&inp.normal).abs() > 0.9
            })
            .min_by(|a, b| {
                centroid
                    .distance_to(&a.centroid)
                    .partial_cmp(&centroid.distance_to(&b.centroid))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some(input) = best_match {
            let ancestor_id = brep.name_pool.intern(input.name.clone(), None);
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors: vec![ancestor_id],
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_modified(input.name.clone(), new_name);
            mod_idx += 1;
        } else {
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
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
            None,
        );
        assert!(matches!(result, Err(ShellError::NoFaces)));
    }

    #[test]
    fn shell_cube_no_removal_produces_result() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let result = shell(
            &brep,
            0.1,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            Tolerance::default(),
            None,
        )
        .unwrap();
        // Outer (6) + inner (6) = 12 faces.
        assert_eq!(result.brep.faces.len(), 12);
    }

    #[test]
    fn shell_cube_one_face_removed() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = shell(
            &brep,
            0.1,
            ShellDirection::Inward,
            &[first_face],
            &[],
            &ShellOptions::default(),
            Tolerance::default(),
            None,
        )
        .unwrap();
        // Outer (5) + inner (5) + side walls (4) = 14 faces.
        assert_eq!(result.brep.faces.len(), 14);
    }

    #[test]
    fn shell_variable_produces_result() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = shell_variable(
            &brep,
            0.1,
            &[(first_face, 0.2)],
            ShellDirection::Inward,
            &[],
            &ShellOptions::default(),
            Tolerance::default(),
            None,
        )
        .unwrap();
        assert_eq!(result.brep.faces.len(), 12);
    }
}
