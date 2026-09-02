//! Feature execution engine — dispatches each [`Feature`] variant to the
//! appropriate popkern operation and records naming history.

use popkern_core::geometry::Point3;
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use std::collections::HashMap;

use crate::error::ReplayError;
use crate::feature::Feature;
use crate::node::BodyId;
use crate::parameter::ParameterTable;
use crate::replay_report::{FeatureEvidence, face_census};

#[path = "replay_fns.rs"]
pub(crate) mod replay_fns;
use replay_fns::{
    chamfer_specs, execute_extrude, execute_sketch, extract_edge_curve, extract_face_outer_curve,
    fillet_specs, resolve_param, to_bool_operation, to_shell_direction, to_shell_options,
};

use crate::resolve::{
    resolve_axis, resolve_edge, resolve_edges, resolve_face, resolve_faces, resolve_revolve_angle,
};

/// Side effects produced by a successful feature execution.
///
pub(crate) use crate::replay_side_effects::ExecuteSideEffects;

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
#[allow(clippy::too_many_arguments)]
pub(crate) fn execute_feature(
    brep: &mut BRep,
    feature: &Feature,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    body_breps: &mut HashMap<BodyId, BRep>,
    param_table: &ParameterTable,
    exact_profile: Option<&crate::replay_exact_extrude::ExactProfileWire>,
    merge_host: Option<BodyId>,
) -> Result<ExecuteSideEffects, ReplayError> {
    match feature {
        Feature::Sketch { plane, sketch_data } => {
            let host_brep = merge_host.and_then(|h| body_breps.get(&h));
            execute_sketch(
                brep,
                plane,
                sketch_data,
                host_brep,
                tracking,
                tolerance,
                param_table,
            )?;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::RawBody { snapshot, .. } => {
            crate::feature_snapshot::materialize_snapshot(brep, snapshot, tracking, tolerance)?;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Extrude {
            face,
            direction,
            end_condition,
            draft_angle,
            second_direction,
            ..
        } => {
            // Analytic profile + blind end: the exact wire sweep
            // `popkern-part` uses (one code path, one topology).
            let host_brep = merge_host.and_then(|h| body_breps.get(&h));
            if let Some(exact) = exact_profile
                && crate::replay_exact_extrude::try_execute_exact(
                    brep,
                    exact,
                    face,
                    direction.as_ref(),
                    end_condition,
                    draft_angle.as_ref(),
                    second_direction.as_ref(),
                    host_brep,
                    tracking,
                    tolerance,
                    param_table,
                )?
            {
                return Ok(ExecuteSideEffects::empty());
            }
            let mut evidence = FeatureEvidence::empty();
            if exact_profile.is_some_and(|e| e.circular) {
                // A full-circle profile on the generic face extrude
                // splits its lateral into co-cylindrical (or with
                // draft, co-conical) halves — recorded, never silent
                // (PK-387/PK-395; reached only by combinations the
                // exact sweep declines, e.g. draft + second direction).
                evidence
                    .fallbacks
                    .push(crate::replay_report::FallbackRecord {
                        source: crate::replay_report::RecordSource::Sweep,
                        mechanism: "exact periodic profile sweep -> generic face extrude"
                            .to_string(),
                        reason: "the exact sweep declined this combination; the full-circle \
                             profile extrudes as split laterals"
                            .to_string(),
                    });
            }
            execute_extrude(
                brep,
                face,
                direction.as_ref(),
                end_condition,
                draft_angle.as_ref(),
                second_direction.as_ref(),
                host_brep,
                tracking,
                tolerance,
                param_table,
            )?;
            Ok(ExecuteSideEffects::reported(evidence))
        }
        Feature::Fillet {
            edges,
            radius,
            recipe,
        } => {
            let edge_handles = resolve_edges(edges, brep)?;
            let edges_requested = edge_handles.len();
            let specs = fillet_specs(edge_handles, *radius, recipe, param_table)?;
            let options = recipe.to_kernel_options();
            let pre_faces = face_census(brep);
            let result =
                popkern_blend::fillet_edges(brep, &specs, &options, *tolerance, Some(tracking))
                    .map_err(|e| ReplayError::FeatureFailed {
                        label: "fillet".into(),
                        error: crate::failure::FeatureFailure::blend(&e, Some(brep)).into(),
                    })?;
            let evidence = FeatureEvidence::blend(
                &result,
                &pre_faces,
                edges_requested,
                Some(options.continuity),
            );
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
        }
        Feature::Chamfer {
            edges,
            distance,
            kind,
        } => {
            let edge_handles = resolve_edges(edges, brep)?;
            let edges_requested = edge_handles.len();
            let specs = chamfer_specs(brep, edge_handles, *distance, kind, param_table)?;
            let pre_faces = face_census(brep);
            let result = popkern_blend::chamfer_edges(
                brep,
                &specs,
                &popkern_blend::ChamferOptions::default(),
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "chamfer".into(),
                error: crate::failure::FeatureFailure::blend(&e, Some(brep)).into(),
            })?;
            let evidence = FeatureEvidence::blend(&result, &pre_faces, edges_requested, None);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
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
                        error: format!("tool body {tool_body:?} not found").into(),
                    })?;

            let bool_op = to_bool_operation(*operation);
            // Drain stale registry entries so the firings taken after
            // the call belong to this boolean.
            let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
            let result =
                popkern_bool::boolean(bool_op, brep.clone(), tool_brep, *tolerance, Some(tracking))
                    .map_err(|e| ReplayError::FeatureFailed {
                        label: "boolean".into(),
                        error: crate::failure::FeatureFailure::boolean(
                            &format!("{operation:?}").to_lowercase(),
                            &e,
                        )
                        .into(),
                    })?;
            let mut evidence = FeatureEvidence::empty();
            evidence.record_boolean_fallbacks(
                popkern_bool::merge_fallback::take_merge_fallback_events(),
            );
            *brep = result;
            Ok(ExecuteSideEffects {
                consumed_bodies: vec![*tool_body],
                evidence,
            })
        }
        Feature::Shell {
            faces_to_remove,
            thickness,
            direction,
            join_style,
            trim_policy,
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
                &to_shell_options(join_style, trim_policy),
                *tolerance,
                Some(tracking),
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "shell".into(),
                error: crate::failure::FeatureFailure::shell("shell", &e).into(),
            })?;
            let evidence = FeatureEvidence::shell(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
        }
        Feature::Transform { motion } => {
            let t = motion.resolve(brep, param_table)?;
            *brep = brep.transform(&t).map_err(|e| ReplayError::FeatureFailed {
                label: "transform".into(),
                error: crate::failure::FeatureFailure::other(e.to_string()).into(),
            })?;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::Thicken {
            thickness,
            join_style,
            trim_policy,
        } => {
            let thick = resolve_param(param_table, *thickness)?;
            let result = popkern_shell::thicken::thicken(
                brep,
                thick,
                &to_shell_options(join_style, trim_policy),
                *tolerance,
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "thicken".into(),
                error: crate::failure::FeatureFailure::shell("thicken", &e).into(),
            })?;
            let evidence = FeatureEvidence::shell(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
        }
        Feature::OffsetSolid {
            distance,
            join_style,
            trim_policy,
        } => {
            let dist = resolve_param(param_table, *distance)?;
            let result = popkern_shell::offset::offset_solid(
                brep,
                dist,
                &to_shell_options(join_style, trim_policy),
                *tolerance,
            )
            .map_err(|e| ReplayError::FeatureFailed {
                label: "offset_solid".into(),
                error: crate::failure::FeatureFailure::shell("offset_solid", &e).into(),
            })?;
            let evidence = FeatureEvidence::shell(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
        }
        Feature::Draft {
            faces,
            angle,
            neutral_plane,
        } => {
            let ang = resolve_param(param_table, *angle)?;
            let face_handles = resolve_faces(faces, brep)?;
            let (neutral_origin, neutral_normal) = neutral_plane.resolve(brep)?;
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
                error: crate::failure::FeatureFailure::draft(&e).into(),
            })?;
            let evidence = FeatureEvidence::draft(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
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
                    error: crate::failure::FeatureFailure::sweep("sweep", &e).into(),
                })?;
            let evidence = FeatureEvidence::sweep(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
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
                    error: crate::failure::FeatureFailure::sweep("loft", &e).into(),
                })?;
            let evidence = FeatureEvidence::sweep(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
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
                error: crate::failure::FeatureFailure::sweep("pipe", &e).into(),
            })?;
            let evidence = FeatureEvidence::sweep(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
        }
        Feature::Revolve {
            face,
            axis,
            end_condition,
            ..
        } => {
            let face_handle = resolve_face(face, brep)?;
            let (origin, direction) = resolve_axis(axis, brep)?;

            // Build context for UpToSurface resolution (centroid of profile face).
            let positions = brep.face_outer_loop_positions(face_handle);
            let profile_centroid = if positions.is_empty() {
                Point3::new(0.0, 0.0, 0.0)
            } else {
                let n = positions.len() as f64;
                Point3::new(
                    positions.iter().map(|p| p.x).sum::<f64>() / n,
                    positions.iter().map(|p| p.y).sum::<f64>() / n,
                    positions.iter().map(|p| p.z).sum::<f64>() / n,
                )
            };
            let revolve_ctx = crate::resolve::RevolveContext {
                axis_origin: origin,
                axis_dir: direction,
                profile_centroid,
                brep,
                tolerance,
            };
            let angle_val = resolve_revolve_angle(end_condition, param_table, Some(&revolve_ctx))?;

            // An analytic sketch revolves its EXACT wire (arcs stay
            // arcs); the compiled face's loop is only a vertex polyline
            // (PK-408: a circle profile lost its curvature entirely).
            let profile = match exact_profile {
                Some(exact) => {
                    let placement =
                        exact.plane.u_axis * exact.center[0] + exact.plane.v_axis * exact.center[1];
                    popkern_sweep::SweepProfile::Wire(popkern_core::geometry::NurbsCurve {
                        control_points: exact
                            .wire
                            .control_points
                            .iter()
                            .map(|p| *p + placement)
                            .collect(),
                        ..exact.wire.clone()
                    })
                }
                None => popkern_sweep::SweepProfile::Face {
                    brep: Box::new(brep.clone()),
                    face: face_handle,
                },
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
                error: crate::failure::FeatureFailure::sweep("revolve", &e).into(),
            })?;
            let evidence = FeatureEvidence::sweep(&result.warnings);
            *brep = result.brep;
            Ok(ExecuteSideEffects::reported(evidence))
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
            // Pattern features are executed by the pattern_execute module, not
            // by execute_feature.  Return empty side effects so the replay loop
            // can record a success result.
            Ok(ExecuteSideEffects::empty())
        }
    }
}

#[cfg(test)]
#[path = "replay_tests.rs"]
mod tests;
