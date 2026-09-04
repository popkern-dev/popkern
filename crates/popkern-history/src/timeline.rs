//! Timeline management and dependency inference for the feature tree.
//!
//! This module provides:
//! - Automatic dependency inference by walking [`Feature`] variants for
//!   [`EntityName`] references and mapping them back to the features that
//!   produced them (via [`OperationId`]).
//! - Timeline insertion and reordering with dependency validation.

use std::collections::HashMap;

use popkern_core::naming::{EntityName, OperationId};
use slotmap::SlotMap;

use crate::body::BodyManager;
use crate::feature::{AxisDesc, ExtrudeEnd, Feature, RevolveEnd, SketchPlaneDesc, Vec3Desc};
use crate::feature_params::feature_param_ids;
use crate::node::{BodyId, FeatureId, FeatureNode};
use crate::parameter::{ParamScope, ParameterTable};

// ---------------------------------------------------------------------------
// EntityName extraction
// ---------------------------------------------------------------------------

/// Collect [`EntityName`] references from a [`Vec3Desc`].
fn vec3_desc_entity_names(desc: &Vec3Desc) -> Vec<&EntityName> {
    match desc {
        Vec3Desc::Literal(_) => Vec::new(),
        Vec3Desc::FaceNormal(name) | Vec3Desc::EdgeDirection(name) => vec![name],
        Vec3Desc::Negated(inner) => vec3_desc_entity_names(inner),
    }
}

/// Collect [`EntityName`] references from an [`AxisDesc`].
fn axis_desc_entity_names(desc: &AxisDesc) -> Vec<&EntityName> {
    match desc {
        AxisDesc::Edge(name) => vec![name],
        AxisDesc::Literal { .. } => Vec::new(),
    }
}

/// Collect [`EntityName`] references from an [`ExtrudeEnd`].
fn extrude_end_entity_names(end: &ExtrudeEnd) -> Vec<&EntityName> {
    match end {
        ExtrudeEnd::UpToSurface { face } | ExtrudeEnd::OffsetFromSurface { face, .. } => {
            vec![face]
        }
        ExtrudeEnd::UpToVertex { vertex } => vec![vertex],
        ExtrudeEnd::Blind { .. } | ExtrudeEnd::ThroughAll | ExtrudeEnd::MidPlane { .. } => {
            Vec::new()
        }
    }
}

/// Collect [`EntityName`] references from a [`RevolveEnd`].
fn revolve_end_entity_names(end: &RevolveEnd) -> Vec<&EntityName> {
    match end {
        RevolveEnd::UpToSurface { face } => vec![face],
        RevolveEnd::Blind { .. } | RevolveEnd::FullCircle | RevolveEnd::MidPlane { .. } => {
            Vec::new()
        }
    }
}

/// Collect [`EntityName`] references from a [`SketchPlaneDesc`].
fn sketch_plane_entity_names(desc: &SketchPlaneDesc) -> Vec<&EntityName> {
    match desc {
        SketchPlaneDesc::XY
        | SketchPlaneDesc::XZ
        | SketchPlaneDesc::YZ
        | SketchPlaneDesc::BaseOffset { .. } => Vec::new(),
        SketchPlaneDesc::OnFace(name) => vec![name],
        SketchPlaneDesc::OnFaceOffset(name, _) => vec![name],
        SketchPlaneDesc::Custom { .. } => Vec::new(),
    }
}

/// Extract all [`EntityName`] references from a [`Feature`] variant.
///
/// This walks the Feature enum recursively, collecting every stable name that
/// the feature references.  Boolean features reference a [`BodyId`] rather than
/// entity names, so they return an empty list.
pub(crate) fn feature_entity_names(feature: &Feature) -> Vec<&EntityName> {
    match feature {
        Feature::Sketch { plane, .. } => sketch_plane_entity_names(plane),
        Feature::Extrude {
            face,
            direction,
            end_condition,
            second_direction,
            ..
        } => {
            let mut names = vec![face];
            if let Some(d) = direction {
                names.extend(vec3_desc_entity_names(d));
            }
            names.extend(extrude_end_entity_names(end_condition));
            if let Some(sd) = second_direction {
                names.extend(vec3_desc_entity_names(&sd.direction));
                names.extend(extrude_end_entity_names(&sd.end_condition));
            }
            names
        }
        Feature::Fillet { edges, .. } => edges.iter().collect(),
        Feature::Transform { motion } => motion.entity_names(),
        Feature::Chamfer { edges, kind, .. } => {
            let mut names: Vec<&EntityName> = edges.iter().collect();
            names.extend(kind.reference_face());
            names
        }
        Feature::Boolean { .. }
        | Feature::Thicken { .. }
        | Feature::OffsetSolid { .. }
        | Feature::RawBody { .. } => Vec::new(),
        Feature::Shell {
            faces_to_remove, ..
        } => faces_to_remove.iter().collect(),
        Feature::Draft {
            faces,
            neutral_plane,
            ..
        } => {
            let mut names: Vec<&EntityName> = faces.iter().collect();
            names.extend(neutral_plane.entity_name());
            names
        }
        Feature::Sweep {
            face, path_edge, ..
        } => vec![face, path_edge],
        Feature::Loft { sections, .. } => sections.iter().collect(),
        Feature::Pipe { path_edge, .. } => vec![path_edge],
        Feature::HelicalSweep { face, axis, .. } => {
            let mut names = vec![face];
            names.extend(axis_desc_entity_names(axis));
            names
        }
        Feature::Revolve {
            face,
            axis,
            end_condition,
            ..
        } => {
            let mut names = vec![face];
            names.extend(axis_desc_entity_names(axis));
            names.extend(revolve_end_entity_names(end_condition));
            names
        }
        Feature::Reference { query, .. } => query.entity_names(),
        Feature::LinearPattern {
            primary, secondary, ..
        } => {
            let mut names = vec3_desc_entity_names(&primary.direction);
            if let Some(sec) = secondary {
                names.extend(vec3_desc_entity_names(&sec.direction));
            }
            names
        }
        Feature::CircularPattern { axis, .. } => axis_desc_entity_names(axis),
        Feature::CurvePattern { curve, .. } => vec![curve],
        Feature::Mirror { plane, .. } => match plane {
            crate::pattern::MirrorPlaneDesc::OnFace(name) => vec![name],
            crate::pattern::MirrorPlaneDesc::XY
            | crate::pattern::MirrorPlaneDesc::XZ
            | crate::pattern::MirrorPlaneDesc::YZ
            | crate::pattern::MirrorPlaneDesc::Custom { .. } => Vec::new(),
        },
    }
}

// ---------------------------------------------------------------------------
// Dependency inference
// ---------------------------------------------------------------------------

/// Infer dependencies by mapping [`EntityName`] references to the
/// [`FeatureId`]s that produced them.
///
/// For each [`EntityName`] referenced by `feature`, the function extracts the
/// [`OperationId`] and looks it up in `op_index` for an O(1) mapping to the
/// originating [`FeatureId`].  Pattern features additionally record direct
/// `FeatureId` dependencies on their seed features.  The resulting list is
/// deduplicated.
pub(crate) fn infer_dependencies(
    feature: &Feature,
    op_index: &HashMap<OperationId, FeatureId>,
) -> Vec<FeatureId> {
    let entity_names = feature_entity_names(feature);
    let mut deps: Vec<FeatureId> = Vec::new();

    for name in entity_names {
        if let Some(&fid) = op_index.get(&name.operation)
            && !deps.contains(&fid)
        {
            deps.push(fid);
        }
    }

    // Add direct seed feature dependencies for pattern features.
    match feature {
        Feature::LinearPattern { seed_features, .. }
        | Feature::CircularPattern { seed_features, .. }
        | Feature::CurvePattern { seed_features, .. }
        | Feature::Mirror { seed_features, .. } => {
            for &seed_id in seed_features {
                if !deps.contains(&seed_id) {
                    deps.push(seed_id);
                }
            }
        }
        _ => {}
    }

    deps
}

// ---------------------------------------------------------------------------
// Body resolution
// ---------------------------------------------------------------------------

/// Whether a feature starts a new body when it has no dependencies
/// (a sketch on a base or custom plane, a raw body).
pub(crate) fn starts_new_body(feature: &Feature) -> bool {
    // Every sketch is its own body: a sketch executes by REPLACING its
    // body's BRep with the profile face, so a sketch placed in the body
    // it references (an on-face sketch inheriting the host through its
    // dependency) erased the host — the flagship plate vanished under a
    // same-body boss (PK-387). The host relation is recovered by the
    // Boss/Cut auto-merge through the on-face reference.
    //
    // A pipe's only dependency is its path edge; inheriting that body
    // would REPLACE the path solid with the tube (PK-440). The tube is
    // a new body; the path body is referenced, not consumed.
    matches!(
        feature,
        Feature::Sketch { .. } | Feature::RawBody { .. } | Feature::Pipe { .. }
    )
}

/// Determine the [`BodyId`] for a newly created feature.
///
/// Rules:
/// - A [`Feature::Sketch`] with no dependencies creates a new body.
/// - Pattern features ([`Feature::LinearPattern`], [`Feature::CircularPattern`],
///   [`Feature::CurvePattern`], [`Feature::Mirror`]) are always assigned
///   `BodyId(0)` — they operate on the host body and the concrete body is
///   resolved at replay time.
/// - All other features inherit the body of their first dependency.
/// - If there are no dependencies, fall back to `BodyId(0)`.
pub(super) fn resolve_body(
    feature: &Feature,
    depends_on: &[FeatureId],
    body_manager: &mut BodyManager,
    features: &SlotMap<FeatureId, FeatureNode>,
) -> BodyId {
    // A base sketch (no planar face dependency) and a body snapshot
    // always start a new body.
    if starts_new_body(feature) {
        return body_manager.create_body();
    }

    // Pattern features always land on BodyId(0); the real body is resolved at
    // replay time when the seed features are evaluated.
    if matches!(
        feature,
        Feature::LinearPattern { .. }
            | Feature::CircularPattern { .. }
            | Feature::CurvePattern { .. }
            | Feature::Mirror { .. }
    ) {
        if !body_manager.bodies().contains(&BodyId(0)) {
            body_manager.create_body();
        }
        return BodyId(0);
    }

    // Inherit from the first dependency's body.
    if let Some(&first_dep) = depends_on.first()
        && let Some(node) = features.get(first_dep)
    {
        return node.body;
    }

    // No dependencies: fall back to BodyId(0), creating it if needed.
    if !body_manager.bodies().contains(&BodyId(0)) {
        // Insert BodyId(0) explicitly without incrementing next_id past 0
        // (create_body would assign 0 if next_id == 0).
        body_manager.create_body();
    }
    BodyId(0)
}

// ---------------------------------------------------------------------------
// Parameter registration
// ---------------------------------------------------------------------------

/// Update parameter scopes and register qualified names after feature insertion.
pub(super) fn register_feature_params(
    table: &mut ParameterTable,
    feature: &Feature,
    feature_id: FeatureId,
    label: &str,
) {
    let param_fields = feature_param_ids(feature);
    let prefix = crate::labels::qualified_prefix(label);
    for (field_name, param_id) in param_fields {
        if let Some(param) = table.get_mut(param_id) {
            param.scope = ParamScope::Feature(feature_id);
            param.name = field_name.to_string();
        }
        // Update name_index with qualified name.
        let qualified = format!("{prefix}.{field_name}");
        // Remove old entry if one exists for this param_id.
        table.update_name_index(param_id, &qualified);
    }
}

#[path = "timeline_ops.rs"]
mod timeline_ops;

#[cfg(test)]
#[path = "timeline_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "timeline_pattern_tests.rs"]
mod pattern_tests;
