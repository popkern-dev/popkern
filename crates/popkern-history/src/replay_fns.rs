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
use crate::feature_sketch::SketchPayload;
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
            error: format!("parameter {id:?} not found").into(),
        })
}

/// Execute a sketch feature in either dialect: place it on its plane,
/// assign stable names, and convert to a BRep lamina. The constraint
/// arm is deserialised and re-solved so parametric changes take
/// effect; the analytic arm has its bindings resolved against the
/// parameter table and is compiled exactly (never solved).
pub(super) fn execute_sketch(
    brep: &mut BRep,
    plane_desc: &SketchPlaneDesc,
    payload: &SketchPayload,
    host_brep: Option<&BRep>,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<(), ReplayError> {
    // An on-face plane lives on the host body (this sketch is its own
    // tool body); base planes resolve anywhere.
    let plane = resolve_sketch_plane(plane_desc, host_brep.unwrap_or(brep), param_table)?;
    let failed = |error: String| ReplayError::FeatureFailed {
        label: "sketch".into(),
        error: error.into(),
    };

    let result = match payload {
        crate::SketchPayload::Constraint(sketch_data) => {
            let mut sketch: popkern_constraint::Sketch = serde_json::from_str(sketch_data)
                .map_err(|e| failed(format!("failed to deserialize sketch: {e}")))?;

            // Re-solve the sketch so parametric changes take effect on replay.
            sketch
                .solve()
                .map_err(|e| failed(format!("constraint solve failed: {e}")))?;

            let mut sketch_on_plane = popkern_sketch::SketchOnPlane {
                sketch,
                plane,
                element_names: std::collections::HashMap::new(),
            };
            sketch_on_plane.assign_element_names(tracking.operation);
            sketch_on_plane
                .to_brep(*tolerance, Some(tracking))
                .map_err(|e| failed(e.to_string()))?
        }
        SketchPayload::Analytic(profile) => {
            let profile = crate::feature_sketch::resolve_profile(profile, param_table)?;
            profile
                .to_brep(plane, *tolerance, Some(tracking))
                .map_err(|e| failed(format!("analytic profile {}: {e}", profile.describe())))?
        }
    };

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

/// Convert the serialisable shell option kinds to runtime
/// [`popkern_shell::ShellOptions`] (PK-328: history features carry
/// every knob so replay is parameter-complete).
pub(super) fn to_shell_options(
    join_style: &crate::feature::JoinStyleKind,
    trim_policy: &crate::feature::TrimPolicyKind,
) -> popkern_shell::ShellOptions {
    let join = match join_style {
        crate::feature::JoinStyleKind::Arc => popkern_shell::JoinStyle::Arc,
        crate::feature::JoinStyleKind::Intersection => popkern_shell::JoinStyle::Intersection,
    };
    let trim = match trim_policy {
        crate::feature::TrimPolicyKind::Trim => popkern_shell::TrimPolicy::Trim,
        crate::feature::TrimPolicyKind::Decline => popkern_shell::TrimPolicy::Decline,
    };
    popkern_shell::ShellOptions::default()
        .with_join_style(join)
        .with_trim_policy(trim)
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
pub(crate) fn extract_plane_from_face(
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
pub(crate) fn face_centroid(brep: &BRep, face_handle: FaceHandle) -> Point3 {
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
            error: format!("draft angle failed: {e}").into(),
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
    host_brep: Option<&BRep>,
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
        host_brep,
        brep,
        param_table,
        tolerance,
    )?;

    // Primary extrude.
    let extrude_result = brep
        .extrude(face_handle, dir, dist, Some(tracking))
        .map_err(|e| ReplayError::FeatureFailed {
            label: "extrude".into(),
            error: crate::failure::FeatureFailure::sweep("extrude", &e).into(),
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
    // Mid-plane on the generic path: the resolved distance is the half
    // extent; grow the other half from the base face (was a one-way
    // extrude of half the height before PK-387).
    if second_direction.is_none() && matches!(end_condition, ExtrudeEnd::MidPlane { .. }) {
        let mut accounted: HashSet<FaceHandle> =
            extrude_result.side_faces.iter().copied().collect();
        accounted.insert(extrude_result.cap_face);
        let base_face = brep
            .faces
            .keys()
            .find(|fh| !accounted.contains(fh))
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: String::new(),
                error: "could not find base face for mid-plane extrude".into(),
            })?;
        brep.extrude(base_face, -dir, dist, Some(tracking))
            .map_err(|e| ReplayError::FeatureFailed {
                label: "extrude".into(),
                error: crate::failure::FeatureFailure::sweep(
                    "extrude",
                    &format!("mid-plane second half extrude failed: {e}"),
                )
                .into(),
            })?;
    }

    if let Some(second) = second_direction {
        let dir2 = resolve_direction(&second.direction, brep)?;
        let dist2 = resolve_extrude_distance(
            &second.end_condition,
            centroid,
            dir2,
            host_brep,
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
                error: crate::failure::FeatureFailure::sweep(
                    "extrude",
                    &format!("second direction extrude failed: {e}"),
                )
                .into(),
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

/// Kernel fillet specs for `edges`: the feature's `radius` parameter
/// combined with the recipe's radius law, every value resolved through
/// the parameter table.
pub(super) fn fillet_specs(
    edges: Vec<popkern_core::topology::EdgeHandle>,
    radius: crate::parameter::ParameterId,
    recipe: &crate::feature_blend::FilletRecipeDesc,
    param_table: &ParameterTable,
) -> Result<Vec<popkern_blend::FilletSpec>, ReplayError> {
    use crate::feature_blend::FilletRadiusLaw;
    let start = resolve_param(param_table, radius)?;
    let law = match &recipe.radius_law {
        FilletRadiusLaw::Constant => popkern_blend::FilletRadius::Constant(start),
        FilletRadiusLaw::Linear { end } => popkern_blend::FilletRadius::Variable {
            start,
            end: resolve_param(param_table, *end)?,
        },
        FilletRadiusLaw::ByStations(stations) => {
            let points = stations
                .iter()
                .map(|s| Ok((s.t, resolve_param(param_table, s.radius)?)))
                .collect::<Result<Vec<(f64, f64)>, ReplayError>>()?;
            popkern_blend::FilletRadius::ByPoints(points)
        }
    };
    Ok(edges
        .into_iter()
        .map(|edge| popkern_blend::FilletSpec {
            edge,
            radius: law.clone(),
        })
        .collect())
}

/// Kernel chamfer specs for `edges`: the feature's `distance`
/// parameter combined with the chamfer kind, reference faces resolved
/// by stable name and every value through the parameter table.
pub(super) fn chamfer_specs(
    brep: &BRep,
    edges: Vec<popkern_core::topology::EdgeHandle>,
    distance: crate::parameter::ParameterId,
    kind: &crate::feature_blend::ChamferKindDesc,
    param_table: &ParameterTable,
) -> Result<Vec<popkern_blend::ChamferSpec>, ReplayError> {
    use crate::feature_blend::ChamferKindDesc;
    let dist = resolve_param(param_table, distance)?;
    let kernel_kind = match kind {
        ChamferKindDesc::Equal => popkern_blend::ChamferKind::Equal(dist),
        ChamferKindDesc::TwoDistance { face, distance_2 } => {
            popkern_blend::ChamferKind::TwoDistance {
                face_1: crate::resolve::resolve_face(face, brep)?,
                dist_1: dist,
                dist_2: resolve_param(param_table, *distance_2)?,
            }
        }
        ChamferKindDesc::DistanceAngle { face, angle } => {
            popkern_blend::ChamferKind::DistanceAngle {
                face: crate::resolve::resolve_face(face, brep)?,
                distance: dist,
                angle: resolve_param(param_table, *angle)?,
            }
        }
    };
    Ok(edges
        .into_iter()
        .map(|edge| popkern_blend::ChamferSpec {
            edge,
            kind: kernel_kind.clone(),
        })
        .collect())
}
