#![warn(missing_docs)]
//! popkern-shell: Shell and offset solid operations for the popkern BRep kernel.
//!
//! Provides the [`shell`] and [`shell_variable`] functions for creating
//! thin-walled solids, plus [`offset_solid`] and [`thicken`] operations.

pub mod assemble;
pub mod classify;
pub mod collision;
pub mod edge_recompute;
pub mod error;
pub mod inner_skin;
pub mod offset;
pub mod side_walls;
pub mod skin_map;
pub mod surface_offset;
pub mod thicken;
pub mod types;
pub mod vertex_recompute;

pub use error::{ShellError, ShellWarning};
pub use types::*;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

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
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();
    let mut result = shell_impl(
        brep,
        thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
    )?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, &pre_faces);
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
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();
    let mut result = shell_impl(
        brep,
        default_thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
    )?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, &pre_faces);
    }
    Ok(result)
}

/// Shared implementation for [`shell`] and [`shell_variable`].
fn shell_impl(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    // Stage 1: classify faces.
    let classification =
        classify::classify_faces(brep, thickness, face_thicknesses, direction, removed_faces)?;

    // Stage 2: compute offset surfaces.
    let stage_tol = Tolerance::default().with_linear(tolerance.linear / 4.0);
    let offset_result = surface_offset::compute_offset_surfaces(
        brep,
        &classification,
        options.allow_partial,
        &stage_tol,
    )?;

    // Merge auto-removed faces into the removed list.
    let mut all_removed = classification.removed_faces.clone();
    all_removed.extend_from_slice(&offset_result.auto_removed);

    // Stage 3: build inner skin.
    let (mut inner, skin_map) = inner_skin::build_inner_skin(brep, &all_removed)?;

    // Stage 4: recompute inner skin geometry.
    vertex_recompute::recompute_vertices(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    edge_recompute::recompute_edges(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        brep,
        &stage_tol,
    )?;
    inner_skin::apply_offset_geometry(&mut inner, &offset_result.surfaces, &skin_map)?;

    // Stage 5: assemble outer + inner + side walls.
    let mut result = assemble::assemble_shell(brep, &inner, &all_removed, &skin_map, tolerance)?;

    // Stage 6: recompute pcurves for all faces.
    result
        .recompute_all_pcurves()
        .map_err(|e| ShellError::Topology(e.to_string()))?;

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

/// Assign stable names to entities produced by a shell operation.
///
/// Faces from the pre-operation set are `Modified` (outer skin faces that
/// were offset). New faces are `Generated` (inner skin and side walls).
fn assign_shell_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    pre_faces: &std::collections::HashSet<FaceHandle>,
) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let post_faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;

    for fh in &post_faces {
        if pre_faces.contains(fh) {
            let ancestors: Vec<_> = brep.faces[*fh].name.into_iter().collect();
            let name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors,
            };
            let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
            let name_id = brep.name_pool.intern(name.clone(), criteria);
            brep.faces[*fh].name = Some(name_id);
            ctx.record_generated(name);
            mod_idx += 1;
        } else {
            let name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
            let name_id = brep.name_pool.intern(name.clone(), criteria);
            brep.faces[*fh].name = Some(name_id);
            ctx.record_generated(name);
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
