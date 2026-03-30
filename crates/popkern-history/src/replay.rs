//! Feature execution engine — dispatches each [`Feature`] variant to the
//! appropriate popkern operation and records naming history.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::feature::{Feature, ShellDirectionKind};
use crate::resolve::{
    resolve_direction, resolve_edges, resolve_face, resolve_faces, resolve_sketch_plane,
};

/// Execute a single feature against the current BRep.
///
/// Dispatches on the [`Feature`] variant, calling the appropriate popkern
/// operation. For operations that return a new BRep (fillet, chamfer, shell,
/// draft), the current `brep` is replaced. For operations that mutate in
/// place (extrude), the BRep is updated directly.
///
/// # Errors
///
/// Returns [`ReplayError`] if entity resolution fails or the underlying
/// operation fails.
pub(crate) fn execute_feature(
    brep: &mut BRep,
    feature: &Feature,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    match feature {
        Feature::Sketch { plane, sketch_data } => {
            execute_sketch(brep, plane, sketch_data, tracking, tolerance)
        }
        Feature::Extrude {
            face,
            direction,
            distance,
        } => {
            let face_handle = resolve_face(face, brep)?;
            let dir = resolve_direction(direction, brep)?;
            brep.extrude(face_handle, dir, *distance, Some(tracking))
                .map_err(|e| ReplayError::FeatureFailed {
                    label: "extrude".into(),
                    error: e.to_string(),
                })?;
            Ok(())
        }
        Feature::Fillet { edges, radius } => {
            let edge_handles = resolve_edges(edges, brep)?;
            let specs: Vec<popkern_fillet::FilletSpec> = edge_handles
                .into_iter()
                .map(|edge| popkern_fillet::FilletSpec {
                    edge,
                    radius: popkern_fillet::FilletRadius::Constant(*radius),
                })
                .collect();
            let result = popkern_fillet::fillet_edges(
                brep,
                &specs,
                &popkern_fillet::FilletOptions::default(),
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "fillet".into(),
                error: e.to_string(),
            })?;
            *brep = result.brep;
            Ok(())
        }
        Feature::Chamfer { edges, distance } => {
            let edge_handles = resolve_edges(edges, brep)?;
            let specs: Vec<popkern_fillet::ChamferSpec> = edge_handles
                .into_iter()
                .map(|edge| popkern_fillet::ChamferSpec {
                    edge,
                    kind: popkern_fillet::ChamferKind::Equal(*distance),
                })
                .collect();
            let result = popkern_fillet::chamfer_edges(
                brep,
                &specs,
                &popkern_fillet::ChamferOptions::default(),
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "chamfer".into(),
                error: e.to_string(),
            })?;
            *brep = result.brep;
            Ok(())
        }
        Feature::Boolean { .. } => Err(ReplayError::FeatureFailed {
            label: "boolean".into(),
            error: "multi-body not yet wired".into(),
        }),
        Feature::Shell {
            faces_to_remove,
            thickness,
            direction,
        } => {
            let face_handles = resolve_faces(faces_to_remove, brep)?;
            let shell_dir = to_shell_direction(direction);
            let result = popkern_shell::shell(
                brep,
                *thickness,
                shell_dir,
                &face_handles,
                &[],
                &popkern_shell::ShellOptions::default(),
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "shell".into(),
                error: e.to_string(),
            })?;
            *brep = result.brep;
            Ok(())
        }
        Feature::Draft {
            faces,
            angle,
            neutral_plane,
        } => {
            let face_handles = resolve_faces(faces, brep)?;
            let neutral_face = resolve_face(neutral_plane, brep)?;
            let (neutral_origin, neutral_normal) = extract_plane_from_face(brep, neutral_face)?;
            // Pull direction is the neutral plane normal.
            let pull_direction = neutral_normal;
            let specs: Vec<popkern_draft::DraftSpec> = face_handles
                .into_iter()
                .map(|face| popkern_draft::DraftSpec {
                    face,
                    angle: *angle,
                    pull_direction,
                    neutral: popkern_draft::NeutralReference::Plane {
                        origin: neutral_origin,
                        normal: neutral_normal,
                    },
                    angle_variation: popkern_draft::AngleVariation::Constant,
                })
                .collect();
            let result = popkern_draft::draft_faces(
                brep,
                &specs,
                &popkern_draft::DraftOptions::default(),
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "draft".into(),
                error: e.to_string(),
            })?;
            *brep = result.brep;
            Ok(())
        }
    }
}

/// Execute a sketch feature: deserialize the sketch, place it on a plane,
/// assign stable names, solve, and convert to BRep.
fn execute_sketch(
    brep: &mut BRep,
    plane_desc: &crate::feature::SketchPlaneDesc,
    sketch_data: &str,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    let plane = resolve_sketch_plane(plane_desc, brep)?;

    let mut sketch: popkern_constraint::Sketch =
        serde_json::from_str(sketch_data).map_err(|e| ReplayError::FeatureFailed {
            label: "sketch".into(),
            error: format!("failed to deserialize sketch: {e}"),
        })?;

    // Re-solve the sketch so parametric changes take effect on replay.
    sketch
        .solve()
        .map_err(|e| ReplayError::FeatureFailed {
            label: "sketch".into(),
            error: format!("constraint solve failed: {e}"),
        })?;

    let mut sketch_on_plane = popkern_sketch::SketchOnPlane {
        sketch,
        plane,
        element_names: std::collections::HashMap::new(),
    };
    sketch_on_plane.assign_element_names(tracking.operation);

    let result = sketch_on_plane
        .to_brep(*tolerance, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "sketch".into(),
            error: e.to_string(),
        })?;

    *brep = result.brep;
    Ok(())
}

/// Convert a serialisable [`ShellDirectionKind`] to the runtime
/// [`popkern_shell::ShellDirection`].
fn to_shell_direction(kind: &ShellDirectionKind) -> popkern_shell::ShellDirection {
    match kind {
        ShellDirectionKind::Inward => popkern_shell::ShellDirection::Inward,
        ShellDirectionKind::Outward => popkern_shell::ShellDirection::Outward,
        ShellDirectionKind::Midplane => popkern_shell::ShellDirection::Midplane,
    }
}

/// Extract the plane origin and normal from a face's planar surface.
///
/// Returns an error if the face has no surface or its surface is not a plane.
fn extract_plane_from_face(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Result<(Point3, Vec3), ReplayError> {
    let sof = brep.faces[face]
        .surface
        .as_ref()
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "draft".into(),
            error: "neutral plane face has no surface geometry".into(),
        })?;
    match &sof.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let n = if sof.orientation { *normal } else { -*normal };
            Ok((*origin, n))
        }
        _ => Err(ReplayError::FeatureFailed {
            label: "draft".into(),
            error: "neutral plane face is not planar".into(),
        }),
    }
}

#[cfg(test)]
#[path = "replay_tests.rs"]
mod tests;
