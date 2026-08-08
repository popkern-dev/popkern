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

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Collected face data for ancestry matching during shell naming.
struct InputFaceData {
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
    )?;
    if let Some(ctx) = tracking {
        assign_shell_names(&mut result.brep, ctx, &input_faces, &brep.tolerance);
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
    let stage_tol = Tolerance::default().with_linear(tolerance.linear() / 4.0);
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

    // Mirror the inner skin's material side honestly: winding and
    // orientation flags flip together, so the void skin keeps the
    // interior-left invariant (a flag-only flip ships the PK-244
    // winding-dishonest mirror class).
    inner.reverse_body_orientation();

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

/// Collect face names, centroids, and normals from a BRep for ancestry matching.
fn collect_shell_face_data(brep: &BRep) -> Vec<InputFaceData> {
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
fn assign_shell_names(
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
