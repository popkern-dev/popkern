//! Feature execution engine — dispatches each [`Feature`] variant to the
//! appropriate popkern operation and records naming history.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use std::collections::HashMap;

use crate::error::ReplayError;
use crate::feature::{BoolOpKind, Feature, ShellDirectionKind};
use crate::node::BodyId;
use crate::parameter::{ParameterId, ParameterTable};

/// Side effects produced by a successful feature execution.
///
/// Returned from [`execute_feature`] so that the caller can apply deferred
/// mutations (e.g. removing consumed tool bodies) only on success.
pub(crate) struct ExecuteSideEffects {
    /// Body IDs that were consumed by this feature (e.g. Boolean tool body).
    pub consumed_bodies: Vec<BodyId>,
}

impl ExecuteSideEffects {
    /// No side effects.
    pub fn empty() -> Self {
        Self {
            consumed_bodies: Vec::new(),
        }
    }
}
use crate::resolve::{
    resolve_direction, resolve_edge, resolve_edges, resolve_face, resolve_faces,
    resolve_sketch_plane,
};

/// Execute a single feature against the current BRep.
///
/// Dispatches on the [`Feature`] variant, calling the appropriate popkern
/// operation. For operations that return a new BRep (fillet, chamfer, shell,
/// draft), the current `brep` is replaced. For operations that mutate in
/// place (extrude), the BRep is updated directly.
///
/// The `body_breps` map holds the BReps for all other bodies. Most features
/// ignore it. Boolean features use it to look up (and consume) the tool
/// body's BRep.
///
/// # Errors
///
/// Returns [`ReplayError`] if entity resolution fails or the underlying
/// operation fails.
/// Resolve a [`ParameterId`] to its `f64` value via the parameter table.
fn resolve_param(table: &ParameterTable, id: ParameterId) -> Result<f64, ReplayError> {
    table
        .get(id)
        .map(|p| p.value)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "parameter".into(),
            error: format!("parameter {id:?} not found"),
        })
}

pub(crate) fn execute_feature(
    brep: &mut BRep,
    feature: &Feature,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    body_breps: &mut HashMap<BodyId, BRep>,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    match feature {
        Feature::Sketch { plane, sketch_data } => {
            execute_sketch(brep, plane, sketch_data, tracking, tolerance, param_table)?;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Extrude {
            face,
            direction,
            distance,
        } => {
            let face_handle = resolve_face(face, brep)?;
            let dir = resolve_direction(direction, brep)?;
            let dist = resolve_param(param_table, *distance)?;
            brep.extrude(face_handle, dir, dist, Some(tracking))
                .map_err(|e| ReplayError::FeatureFailed {
                    label: "extrude".into(),
                    error: e.to_string(),
                })?;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Fillet { edges, radius } => {
            let r = resolve_param(param_table, *radius)?;
            let edge_handles = resolve_edges(edges, brep)?;
            let specs: Vec<popkern_fillet::FilletSpec> = edge_handles
                .into_iter()
                .map(|edge| popkern_fillet::FilletSpec {
                    edge,
                    radius: popkern_fillet::FilletRadius::Constant(r),
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
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Chamfer { edges, distance } => {
            let dist = resolve_param(param_table, *distance)?;
            let edge_handles = resolve_edges(edges, brep)?;
            let specs: Vec<popkern_fillet::ChamferSpec> = edge_handles
                .into_iter()
                .map(|edge| popkern_fillet::ChamferSpec {
                    edge,
                    kind: popkern_fillet::ChamferKind::Equal(dist),
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
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Boolean {
            operation,
            tool_body,
        } => {
            let tool_brep =
                body_breps
                    .get(tool_body)
                    .cloned()
                    .ok_or_else(|| ReplayError::FeatureFailed {
                        label: "boolean".into(),
                        error: format!("tool body {tool_body:?} not found"),
                    })?;

            let bool_op = to_bool_operation(*operation);
            let result =
                popkern_bool::boolean(bool_op, brep.clone(), tool_brep, *tolerance, Some(tracking))
                    .map_err(|e| ReplayError::FeatureFailed {
                        label: "boolean".into(),
                        error: e.to_string(),
                    })?;
            *brep = result;
            Ok(ExecuteSideEffects {
                consumed_bodies: vec![*tool_body],
            })
        }
        Feature::Shell {
            faces_to_remove,
            thickness,
            direction,
        } => {
            let thick = resolve_param(param_table, *thickness)?;
            let face_handles = resolve_faces(faces_to_remove, brep)?;
            let shell_dir = to_shell_direction(direction);
            let result = popkern_shell::shell(
                brep,
                thick,
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
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Draft {
            faces,
            angle,
            neutral_plane,
        } => {
            let ang = resolve_param(param_table, *angle)?;
            let face_handles = resolve_faces(faces, brep)?;
            let neutral_face = resolve_face(neutral_plane, brep)?;
            let (neutral_origin, neutral_normal) = extract_plane_from_face(brep, neutral_face)?;
            // Pull direction is the neutral plane normal.
            let pull_direction = neutral_normal;
            let specs: Vec<popkern_draft::DraftSpec> = face_handles
                .into_iter()
                .map(|face| popkern_draft::DraftSpec {
                    face,
                    angle: ang,
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
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Sweep {
            face,
            path_edge,
            twist,
        } => {
            let tw = resolve_param(param_table, *twist)?;
            let face_handle = resolve_face(face, brep)?;
            let edge_handle = resolve_edge(path_edge, brep)?;
            let path_curve = extract_edge_curve(brep, edge_handle)?;
            let profile = popkern_sweep::SweepProfile::Face {
                brep: Box::new(brep.clone()),
                face: face_handle,
            };
            let options = popkern_sweep::SweepOptions {
                twist: tw,
                ..popkern_sweep::SweepOptions::default()
            };
            let result = popkern_sweep::sweep(&profile, &path_curve, &options, *tolerance)
                .map_err(|e| ReplayError::FeatureFailed {
                    label: "sweep".into(),
                    error: e.to_string(),
                })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Loft { sections, closed } => {
            let mut loft_sections = Vec::new();
            for section_name in sections {
                let face_handle = resolve_face(section_name, brep)?;
                let curve = extract_face_outer_curve(brep, face_handle)?;
                loft_sections.push(popkern_sweep::LoftSection {
                    curve,
                    position: None,
                });
            }
            let options = popkern_sweep::LoftOptions {
                closed: *closed,
                ..popkern_sweep::LoftOptions::default()
            };
            let result =
                popkern_sweep::loft(&loft_sections, &options, *tolerance).map_err(|e| {
                    ReplayError::FeatureFailed {
                        label: "loft".into(),
                        error: e.to_string(),
                    }
                })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Pipe { path_edge, radius } => {
            let r = resolve_param(param_table, *radius)?;
            let edge_handle = resolve_edge(path_edge, brep)?;
            let path_curve = extract_edge_curve(brep, edge_handle)?;
            let result = popkern_sweep::pipe(
                &path_curve,
                &popkern_sweep::PipeRadius::Constant(r),
                &popkern_sweep::SweepOptions::default(),
                *tolerance,
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "pipe".into(),
                error: e.to_string(),
            })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Revolve {
            face,
            axis_origin,
            axis_direction,
            angle,
        } => {
            let face_handle = resolve_face(face, brep)?;
            let angle_val = resolve_param(param_table, *angle)?;

            let origin = Point3::new(axis_origin[0], axis_origin[1], axis_origin[2]);
            let direction = Vec3::new(axis_direction[0], axis_direction[1], axis_direction[2]);

            let profile = popkern_sweep::SweepProfile::Face {
                brep: Box::new(brep.clone()),
                face: face_handle,
            };

            let result =
                popkern_sweep::revolve(&profile, &origin, &direction, angle_val, *tolerance)
                    .map_err(|e| ReplayError::FeatureFailed {
                        label: "revolve".into(),
                        error: e.to_string(),
                    })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
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
fn to_shell_direction(kind: &ShellDirectionKind) -> popkern_shell::ShellDirection {
    match kind {
        ShellDirectionKind::Inward => popkern_shell::ShellDirection::Inward,
        ShellDirectionKind::Outward => popkern_shell::ShellDirection::Outward,
        ShellDirectionKind::Midplane => popkern_shell::ShellDirection::Midplane,
    }
}

/// Convert a serialisable [`BoolOpKind`] to the runtime
/// [`popkern_bool::BoolOperation`].
fn to_bool_operation(kind: BoolOpKind) -> popkern_bool::BoolOperation {
    match kind {
        BoolOpKind::Union => popkern_bool::BoolOperation::Union,
        BoolOpKind::Intersection => popkern_bool::BoolOperation::Intersection,
        BoolOpKind::Subtraction => popkern_bool::BoolOperation::Subtraction,
    }
}

/// Extract a [`NurbsCurve`] from an edge's curve geometry.
fn extract_edge_curve(brep: &BRep, edge: EdgeHandle) -> Result<NurbsCurve, ReplayError> {
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
fn extract_face_outer_curve(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Result<NurbsCurve, ReplayError> {
    let outer_loop = brep.faces[face].outer_loop;
    let first_he = brep.loops[outer_loop].half_edge;
    let first_edge = brep.half_edges[first_he].edge;
    extract_edge_curve(brep, first_edge)
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
