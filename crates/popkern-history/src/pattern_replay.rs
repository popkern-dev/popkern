//! Feature-replay engine for pattern instances.
//!
//! Provides helpers to clone and remap seed features for a specific pattern
//! instance, and to replay a chain of seed features against a fresh BRep to
//! produce a tool body.

use std::collections::HashMap;

use popkern_core::history::TrackingContext;
use popkern_core::naming::OperationId;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::feature::Feature;
use crate::node::BodyId;
use crate::parameter::ParameterTable;
use crate::replay;

/// Clone a seed feature for a specific pattern instance.
///
/// Walks every [`popkern_core::naming::EntityName`] reference inside the
/// feature and replaces `name.operation` values that appear in `op_remap`
/// with the mapped instance-specific [`OperationId`].
///
/// For `Sketch` features with `OnFace` planes, the face reference is kept
/// as-is (the host face is the same for all instances; the different position
/// comes from the auto-merge target).  This is a simplification that works
/// for the initial implementation.
pub(crate) fn clone_seed_feature(
    feature: &Feature,
    op_remap: &HashMap<OperationId, OperationId>,
) -> Feature {
    let mut cloned = feature.clone();
    remap_feature_entity_names(&mut cloned, op_remap);
    cloned
}

/// Replay a chain of seed features to produce a tool body BRep.
///
/// Creates a fresh [`BRep`] and executes each (feature, operation_id) pair in
/// order, returning the resulting BRep.
///
/// # Errors
///
/// Returns [`ReplayError`] if any feature in the chain fails.
pub(crate) fn replay_seed_chain(
    seed_features: &[(Feature, OperationId)],
    tolerance: &Tolerance,
    param_table: &ParameterTable,
    body_breps: &mut HashMap<BodyId, BRep>,
) -> Result<BRep, ReplayError> {
    let mut brep = BRep::new(*tolerance);
    // The chain's on-face sketch names its host; the chain replays into
    // a fresh tool body, so the host is where planes, exact profiles and
    // up-to targets resolve.
    let host = seed_features.iter().find_map(|(f, _)| match f {
        crate::feature::Feature::Sketch { plane, .. } => {
            crate::auto_merge::host_for_sketch_plane(plane, body_breps)
        }
        _ => None,
    });
    for (feature, op_id) in seed_features {
        let mut ctx = TrackingContext::new(*op_id);
        let exact_profile = match feature {
            crate::feature::Feature::Extrude { face, .. } => {
                crate::replay_exact_extrude::exact_profile_for(
                    seed_features.iter().map(|(f, op)| (f, *op)),
                    face,
                    host.and_then(|h| body_breps.get(&h)).unwrap_or(&brep),
                    param_table,
                )
            }
            _ => None,
        };
        replay::execute_feature(
            &mut brep,
            feature,
            &mut ctx,
            tolerance,
            body_breps,
            param_table,
            exact_profile.as_ref(),
            host,
        )?;
    }
    Ok(brep)
}

// ---------------------------------------------------------------------------
// EntityName remapping
// ---------------------------------------------------------------------------

/// Remap all [`popkern_core::naming::EntityName`] references inside a feature.
///
/// For each `EntityName` whose `operation` field appears in `op_remap`, the
/// operation is replaced with the mapped value.
fn remap_feature_entity_names(feature: &mut Feature, op_remap: &HashMap<OperationId, OperationId>) {
    match feature {
        Feature::Sketch { plane, .. } => {
            remap_sketch_plane(plane, op_remap);
        }
        Feature::Extrude {
            face,
            direction,
            end_condition,
            second_direction,
            ..
        } => {
            remap_name(face, op_remap);
            if let Some(d) = direction {
                remap_vec3_desc(d, op_remap);
            }
            remap_extrude_end(end_condition, op_remap);
            if let Some(sd) = second_direction {
                remap_vec3_desc(&mut sd.direction, op_remap);
                remap_extrude_end(&mut sd.end_condition, op_remap);
            }
        }
        Feature::Fillet { edges, .. } => {
            for name in edges.iter_mut() {
                remap_name(name, op_remap);
            }
        }
        Feature::Chamfer { edges, kind, .. } => {
            for name in edges.iter_mut() {
                remap_name(name, op_remap);
            }
            match kind {
                crate::feature_blend::ChamferKindDesc::Equal => {}
                crate::feature_blend::ChamferKindDesc::TwoDistance { face, .. }
                | crate::feature_blend::ChamferKindDesc::DistanceAngle { face, .. } => {
                    remap_name(face, op_remap);
                }
            }
        }
        Feature::Boolean { .. }
        | Feature::Thicken { .. }
        | Feature::OffsetSolid { .. }
        | Feature::RawBody { .. } => {}
        Feature::Transform { motion } => {
            if let Some(name) = motion.entity_name_mut() {
                remap_name(name, op_remap);
            }
        }
        Feature::Shell {
            faces_to_remove, ..
        } => {
            for name in faces_to_remove.iter_mut() {
                remap_name(name, op_remap);
            }
        }
        Feature::Draft {
            faces,
            neutral_plane,
            ..
        } => {
            for name in faces.iter_mut() {
                remap_name(name, op_remap);
            }
            if let Some(name) = neutral_plane.entity_name_mut() {
                remap_name(name, op_remap);
            }
        }
        Feature::Sweep {
            face, path_edge, ..
        } => {
            remap_name(face, op_remap);
            remap_name(path_edge, op_remap);
        }
        Feature::Loft { sections, .. } => {
            for name in sections.iter_mut() {
                remap_name(name, op_remap);
            }
        }
        Feature::Pipe { path_edge, .. } => {
            remap_name(path_edge, op_remap);
        }
        Feature::Revolve {
            face,
            axis,
            end_condition,
            ..
        } => {
            remap_name(face, op_remap);
            remap_axis_desc(axis, op_remap);
            remap_revolve_end(end_condition, op_remap);
        }
        Feature::Reference { query, .. } => {
            remap_measure_query(query, op_remap);
        }
        Feature::LinearPattern { .. }
        | Feature::CircularPattern { .. }
        | Feature::CurvePattern { .. }
        | Feature::Mirror { .. } => {
            // Pattern features reference seed FeatureIds, not entity names
            // that need remapping.
        }
    }
}

/// Remap a single [`popkern_core::naming::EntityName`].
fn remap_name(
    name: &mut popkern_core::naming::EntityName,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    if let Some(new_op) = op_remap.get(&name.operation) {
        name.operation = *new_op;
    }
}

/// Remap entity names inside a [`crate::feature::SketchPlaneDesc`].
fn remap_sketch_plane(
    plane: &mut crate::feature::SketchPlaneDesc,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    match plane {
        crate::feature::SketchPlaneDesc::OnFace(name) => remap_name(name, op_remap),
        crate::feature::SketchPlaneDesc::OnFaceOffset(name, _) => remap_name(name, op_remap),
        crate::feature::SketchPlaneDesc::XY
        | crate::feature::SketchPlaneDesc::XZ
        | crate::feature::SketchPlaneDesc::YZ
        | crate::feature::SketchPlaneDesc::Custom { .. }
        | crate::feature::SketchPlaneDesc::BaseOffset { .. } => {}
    }
}

/// Remap entity names inside a [`crate::feature::Vec3Desc`].
fn remap_vec3_desc(
    desc: &mut crate::feature::Vec3Desc,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    match desc {
        crate::feature::Vec3Desc::Literal(_) => {}
        crate::feature::Vec3Desc::FaceNormal(name)
        | crate::feature::Vec3Desc::EdgeDirection(name) => remap_name(name, op_remap),
        crate::feature::Vec3Desc::Negated(inner) => remap_vec3_desc(inner, op_remap),
    }
}

/// Remap entity names inside an [`crate::feature::ExtrudeEnd`].
fn remap_extrude_end(
    end: &mut crate::feature::ExtrudeEnd,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    match end {
        crate::feature::ExtrudeEnd::UpToSurface { face } => remap_name(face, op_remap),
        crate::feature::ExtrudeEnd::OffsetFromSurface { face, .. } => {
            remap_name(face, op_remap);
        }
        crate::feature::ExtrudeEnd::UpToVertex { vertex } => remap_name(vertex, op_remap),
        crate::feature::ExtrudeEnd::Blind { .. }
        | crate::feature::ExtrudeEnd::ThroughAll
        | crate::feature::ExtrudeEnd::MidPlane { .. } => {}
    }
}

/// Remap entity names inside an [`crate::feature::AxisDesc`].
fn remap_axis_desc(
    axis: &mut crate::feature::AxisDesc,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    match axis {
        crate::feature::AxisDesc::Edge(name) => remap_name(name, op_remap),
        crate::feature::AxisDesc::Literal { .. } => {}
    }
}

/// Remap entity names inside a [`crate::feature::RevolveEnd`].
fn remap_revolve_end(
    end: &mut crate::feature::RevolveEnd,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    match end {
        crate::feature::RevolveEnd::UpToSurface { face } => remap_name(face, op_remap),
        crate::feature::RevolveEnd::Blind { .. }
        | crate::feature::RevolveEnd::FullCircle
        | crate::feature::RevolveEnd::MidPlane { .. } => {}
    }
}

/// Remap entity names inside a [`crate::measure_query::MeasureQuery`].
fn remap_measure_query(
    query: &mut crate::measure_query::MeasureQuery,
    op_remap: &HashMap<OperationId, OperationId>,
) {
    use crate::measure_query::MeasureQuery;
    match query {
        MeasureQuery::EdgeLength(n)
        | MeasureQuery::EdgeTangent(n, _)
        | MeasureQuery::FaceArea(n)
        | MeasureQuery::FacePerimeter(n)
        | MeasureQuery::FaceNormal(n, _)
        | MeasureQuery::FaceCentroid(n, _)
        | MeasureQuery::FaceRadius(n)
        | MeasureQuery::FaceMinCurvature(n)
        | MeasureQuery::FaceMaxCurvature(n)
        | MeasureQuery::VertexPosition(n, _) => {
            remap_name(n, op_remap);
        }
        MeasureQuery::MinDistance(a, b)
        | MeasureQuery::DihedralAngle(a, b)
        | MeasureQuery::EdgeAngle(a, b) => {
            remap_name(a, op_remap);
            remap_name(b, op_remap);
        }
        MeasureQuery::BodyVolume
        | MeasureQuery::BodySurfaceArea
        | MeasureQuery::BodyCentroid(_)
        | MeasureQuery::BodyBound(_)
        | MeasureQuery::BodyMomentOfInertia(_) => {}
    }
}
