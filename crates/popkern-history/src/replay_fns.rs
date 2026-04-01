//! Helper functions for the feature execution engine.
//!
//! Extracted from `replay.rs` to keep each source file under 500 lines.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use std::collections::HashSet;

use crate::error::ReplayError;
use crate::feature::{
    BoolOpKind, ExtrudeEnd, ExtrudeSecondDir, ShellDirectionKind, SketchPlaneDesc, Vec3Desc,
};
use crate::parameter::{ParameterId, ParameterTable};
use crate::resolve::{
    resolve_direction, resolve_extrude_distance, resolve_face, resolve_sketch_plane,
};

/// Resolve a [`ParameterId`] to its `f64` value via the parameter table.
pub(super) fn resolve_param(table: &ParameterTable, id: ParameterId) -> Result<f64, ReplayError> {
    table
        .get(id)
        .map(|p| p.value)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "parameter".into(),
            error: format!("parameter {id:?} not found"),
        })
}

/// Execute a sketch feature: deserialize the sketch, place it on a plane,
/// assign stable names, solve, and convert to BRep.
pub(super) fn execute_sketch(
    brep: &mut BRep,
    plane_desc: &SketchPlaneDesc,
    sketch_data: &str,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<(), ReplayError> {
    let plane = resolve_sketch_plane(plane_desc, brep, param_table)?;

    let mut sketch: popkern_constraint::Sketch =
        serde_json::from_str(sketch_data).map_err(|e| ReplayError::FeatureFailed {
            label: "sketch".into(),
            error: format!("failed to deserialize sketch: {e}"),
        })?;

    // Re-solve the sketch so parametric changes take effect on replay.
    sketch.solve().map_err(|e| ReplayError::FeatureFailed {
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
pub(super) fn to_shell_direction(kind: &ShellDirectionKind) -> popkern_shell::ShellDirection {
    match kind {
        ShellDirectionKind::Inward => popkern_shell::ShellDirection::Inward,
        ShellDirectionKind::Outward => popkern_shell::ShellDirection::Outward,
        ShellDirectionKind::Midplane => popkern_shell::ShellDirection::Midplane,
    }
}

/// Convert a serialisable [`BoolOpKind`] to the runtime
/// [`popkern_bool::BoolOperation`].
pub(super) fn to_bool_operation(kind: BoolOpKind) -> popkern_bool::BoolOperation {
    match kind {
        BoolOpKind::Union => popkern_bool::BoolOperation::Union,
        BoolOpKind::Intersection => popkern_bool::BoolOperation::Intersection,
        BoolOpKind::Subtraction => popkern_bool::BoolOperation::Subtraction,
    }
}

/// Extract a [`NurbsCurve`] from an edge's curve geometry.
pub(super) fn extract_edge_curve(brep: &BRep, edge: EdgeHandle) -> Result<NurbsCurve, ReplayError> {
    let curve_on_edge =
        brep.edges[edge]
            .curve
            .as_ref()
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: "sweep".into(),
                error: "path edge has no curve geometry".into(),
            })?;
    Ok(curve_on_edge
        .geometry
        .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end)))
}

/// Extract the outer loop of a face as a single [`NurbsCurve`].
///
/// Uses the first edge's curve from the face's outer loop. For simple
/// cross-section faces (single-curve boundary) this is sufficient; complex
/// multi-edge boundaries would need curve concatenation (future enhancement).
pub(super) fn extract_face_outer_curve(
    brep: &BRep,
    face: FaceHandle,
) -> Result<NurbsCurve, ReplayError> {
    let outer_loop = brep.faces[face].outer_loop;
    let first_he = brep.loops[outer_loop].half_edge;
    let first_edge = brep.half_edges[first_he].edge;
    extract_edge_curve(brep, first_edge)
}

/// Extract the plane origin and normal from a face's planar surface.
///
/// Returns an error if the face has no surface or its surface is not a plane.
pub(super) fn extract_plane_from_face(
    brep: &BRep,
    face: FaceHandle,
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

/// Compute the centroid of a face's outer loop vertices.
fn face_centroid(brep: &BRep, face_handle: FaceHandle) -> Point3 {
    let positions = brep.face_outer_loop_positions(face_handle);
    if positions.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let sum_x: f64 = positions.iter().map(|p| p.x).sum();
    let sum_y: f64 = positions.iter().map(|p| p.y).sum();
    let sum_z: f64 = positions.iter().map(|p| p.z).sum();
    let n = positions.len() as f64;
    Point3::new(sum_x / n, sum_y / n, sum_z / n)
}

/// Apply a draft angle to the side faces produced by an extrude.
fn apply_draft_to_side_faces(
    brep: &mut BRep,
    side_faces: &[FaceHandle],
    angle: f64,
    pull_direction: Vec3,
    neutral_origin: Point3,
    tolerance: &Tolerance,
    tracking: &mut TrackingContext,
) -> Result<(), ReplayError> {
    let specs: Vec<popkern_draft::DraftSpec> = side_faces
        .iter()
        .map(|&face| popkern_draft::DraftSpec {
            face,
            angle,
            pull_direction,
            neutral: popkern_draft::NeutralReference::Plane {
                origin: neutral_origin,
                normal: pull_direction,
            },
            angle_variation: popkern_draft::AngleVariation::Constant,
        })
        .collect();

    if !specs.is_empty() {
        let draft_result = popkern_draft::draft_faces(
            brep,
            &specs,
            &popkern_draft::DraftOptions::default(),
            *tolerance,
            Some(tracking),
        )
        .map_err(|e| ReplayError::FeatureFailed {
            label: String::new(),
            error: format!("draft angle failed: {e}"),
        })?;
        *brep = draft_result.brep;
    }
    Ok(())
}

/// Execute an extrude operation with optional draft angle and bidirectional
/// support.
///
/// Steps:
/// 1. Resolve face, direction, distance.
/// 2. Compute centroid (for draft neutral plane and distance resolution).
/// 3. Call `brep.extrude()` and capture `ExtrudeResult`.
/// 4. If `draft_angle` is `Some`, apply draft to side faces.
/// 5. If `second_direction` is `Some`, handle bidirectional extrude.
#[allow(clippy::too_many_arguments)]
pub(super) fn execute_extrude(
    brep: &mut BRep,
    face: &EntityName,
    direction: Option<&Vec3Desc>,
    end_condition: &ExtrudeEnd,
    draft_angle: Option<&ParameterId>,
    second_direction: Option<&ExtrudeSecondDir>,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<(), ReplayError> {
    let face_handle = resolve_face(face, brep)?;
    let dir = match direction {
        Some(d) => resolve_direction(d, brep)?,
        None => {
            // Default to the face normal when no direction is specified.
            let sof = brep.faces[face_handle].surface.as_ref().ok_or_else(|| {
                ReplayError::FeatureFailed {
                    label: "extrude".into(),
                    error: "face has no surface geometry".into(),
                }
            })?;
            match &sof.geometry {
                SurfaceGeometry::Plane { normal, .. } => {
                    if sof.orientation {
                        *normal
                    } else {
                        -*normal
                    }
                }
                _ => {
                    return Err(ReplayError::FeatureFailed {
                        label: "extrude".into(),
                        error: "direction=None only supported for planar faces".into(),
                    });
                }
            }
        }
    };

    // Compute centroid for distance resolution and draft neutral plane.
    let centroid = face_centroid(brep, face_handle);

    let dist = resolve_extrude_distance(
        end_condition,
        centroid,
        dir,
        None, // host_brep not available here; auto-merge supplies it
        brep,
        param_table,
        tolerance,
    )?;

    // Primary extrude.
    let extrude_result = brep
        .extrude(face_handle, dir, dist, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "extrude".into(),
            error: e.to_string(),
        })?;

    // Apply draft angle to primary direction side faces.
    if let Some(draft_param) = draft_angle {
        let angle = resolve_param(param_table, *draft_param)?;
        if angle.abs() > 1e-10 {
            apply_draft_to_side_faces(
                brep,
                &extrude_result.side_faces,
                angle,
                dir,
                centroid,
                tolerance,
                tracking,
            )?;
        }
    }

    // Bidirectional extrude.
    if let Some(second) = second_direction {
        let dir2 = resolve_direction(&second.direction, brep)?;
        let dist2 = resolve_extrude_distance(
            &second.end_condition,
            centroid,
            dir2,
            None,
            brep,
            param_table,
            tolerance,
        )?;

        // Find the base face (the original face, not among cap or side faces).
        let mut accounted: HashSet<FaceHandle> =
            extrude_result.side_faces.iter().copied().collect();
        accounted.insert(extrude_result.cap_face);
        let base_face = brep
            .faces
            .keys()
            .find(|fh| !accounted.contains(fh))
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: String::new(),
                error: "could not find base face for bidirectional extrude".into(),
            })?;

        let result2 = brep
            .extrude(base_face, dir2, dist2, Some(tracking))
            .map_err(|e| ReplayError::FeatureFailed {
                label: "extrude".into(),
                error: format!("second direction extrude failed: {e}"),
            })?;

        // Apply draft to second direction's side faces if requested.
        if let Some(draft2_param) = second.draft_angle {
            let angle2 = resolve_param(param_table, draft2_param)?;
            if angle2.abs() > 1e-10 {
                apply_draft_to_side_faces(
                    brep,
                    &result2.side_faces,
                    angle2,
                    dir2,
                    centroid,
                    tolerance,
                    tracking,
                )?;
            }
        }
    }

    Ok(())
}
