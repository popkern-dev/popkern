//! Feature execution engine — dispatches each [`Feature`] variant to the
//! appropriate popkern operation and records naming history.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use std::collections::HashMap;

use crate::error::ReplayError;
use crate::feature::{AxisDesc, ExtrudeEnd, Feature, RevolveEnd};
use crate::node::BodyId;
use crate::parameter::ParameterTable;

#[path = "replay_fns.rs"]
mod replay_fns;
use replay_fns::{
    execute_sketch, extract_edge_curve, extract_face_outer_curve, extract_plane_from_face,
    resolve_param, to_bool_operation, to_shell_direction,
};

use crate::resolve::{resolve_direction, resolve_edge, resolve_edges, resolve_face, resolve_faces};

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
            end_condition,
            ..
        } => {
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
                        popkern_core::geometry::SurfaceGeometry::Plane { normal, .. } => {
                            if sof.orientation { *normal } else { -*normal }
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
            let dist = match end_condition {
                ExtrudeEnd::Blind { distance } => resolve_param(param_table, *distance)?,
                ExtrudeEnd::MidPlane { half_distance } => {
                    resolve_param(param_table, *half_distance)? * 2.0
                }
                ExtrudeEnd::ThroughAll
                | ExtrudeEnd::UpToSurface { .. }
                | ExtrudeEnd::OffsetFromSurface { .. }
                | ExtrudeEnd::UpToVertex { .. } => {
                    todo!("end condition {:?} not yet implemented", end_condition)
                }
            };
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
            ..
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
            let result =
                popkern_sweep::sweep(&profile, &path_curve, &options, *tolerance, Some(tracking))
                    .map_err(|e| ReplayError::FeatureFailed {
                    label: "sweep".into(),
                    error: e.to_string(),
                })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Loft {
            sections, closed, ..
        } => {
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
            let result = popkern_sweep::loft(&loft_sections, &options, *tolerance, Some(tracking))
                .map_err(|e| ReplayError::FeatureFailed {
                    label: "loft".into(),
                    error: e.to_string(),
                })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Pipe {
            path_edge, radius, ..
        } => {
            let r = resolve_param(param_table, *radius)?;
            let edge_handle = resolve_edge(path_edge, brep)?;
            let path_curve = extract_edge_curve(brep, edge_handle)?;
            let result = popkern_sweep::pipe(
                &path_curve,
                &popkern_sweep::PipeRadius::Constant(r),
                &popkern_sweep::SweepOptions::default(),
                *tolerance,
                Some(tracking),
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
            axis,
            end_condition,
            ..
        } => {
            let face_handle = resolve_face(face, brep)?;
            let angle_val = match end_condition {
                RevolveEnd::Blind { angle } => resolve_param(param_table, *angle)?,
                RevolveEnd::FullCircle => std::f64::consts::TAU,
                RevolveEnd::MidPlane { half_angle } => {
                    resolve_param(param_table, *half_angle)? * 2.0
                }
                RevolveEnd::UpToSurface { .. } => {
                    todo!("RevolveEnd::UpToSurface not yet implemented")
                }
            };

            let (origin, direction) = match axis {
                AxisDesc::Literal {
                    origin: o,
                    direction: d,
                } => (Point3::new(o[0], o[1], o[2]), Vec3::new(d[0], d[1], d[2])),
                AxisDesc::Edge(name) => {
                    let edge_handle = resolve_edge(name, brep)?;
                    let curve_on_edge =
                        brep.edges[edge_handle].curve.as_ref().ok_or_else(|| {
                            ReplayError::FeatureFailed {
                                label: "revolve".into(),
                                error: "axis edge has no curve geometry".into(),
                            }
                        })?;
                    let t = curve_on_edge.t_start;
                    let tangent = curve_on_edge.geometry.tangent_at(t);
                    let pt = curve_on_edge.geometry.point_at(t);
                    (pt, tangent)
                }
            };

            let profile = popkern_sweep::SweepProfile::Face {
                brep: Box::new(brep.clone()),
                face: face_handle,
            };

            let result = popkern_sweep::revolve(
                &profile,
                &origin,
                &direction,
                angle_val,
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "revolve".into(),
                error: e.to_string(),
            })?;
            *brep = result.brep;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Reference { .. } => {
            // Reference features don't modify the BRep. Measurement evaluation
            // and parameter update happen at the tree-replay level.
            Ok(ExecuteSideEffects::empty())
        }
        Feature::LinearPattern { .. }
        | Feature::CircularPattern { .. }
        | Feature::CurvePattern { .. }
        | Feature::Mirror { .. } => {
            todo!("pattern execution not yet implemented")
        }
    }
}

#[cfg(test)]
#[path = "replay_tests.rs"]
mod tests;
