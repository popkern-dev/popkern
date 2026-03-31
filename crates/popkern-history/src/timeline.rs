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

use crate::FeatureTree;
use crate::body::BodyManager;
use crate::error::HistoryError;
use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc, feature_param_ids};
use crate::node::{BodyId, FeatureId, FeatureNode};
use crate::parameter::{ParamScope, ParameterTable};
use crate::undo::TreeCommand;

// ---------------------------------------------------------------------------
// EntityName extraction
// ---------------------------------------------------------------------------

/// Collect [`EntityName`] references from a [`Vec3Desc`].
fn vec3_desc_entity_names(desc: &Vec3Desc) -> Vec<&EntityName> {
    match desc {
        Vec3Desc::Literal(_) => Vec::new(),
        Vec3Desc::FaceNormal(name) => vec![name],
        Vec3Desc::Negated(inner) => vec3_desc_entity_names(inner),
    }
}

/// Collect [`EntityName`] references from a [`SketchPlaneDesc`].
fn sketch_plane_entity_names(desc: &SketchPlaneDesc) -> Vec<&EntityName> {
    match desc {
        SketchPlaneDesc::XY | SketchPlaneDesc::XZ | SketchPlaneDesc::YZ => Vec::new(),
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
            face, direction, ..
        } => {
            let mut names = vec![face];
            names.extend(vec3_desc_entity_names(direction));
            names
        }
        Feature::Fillet { edges, .. } | Feature::Chamfer { edges, .. } => edges.iter().collect(),
        Feature::Boolean { .. } => Vec::new(),
        Feature::Shell {
            faces_to_remove, ..
        } => faces_to_remove.iter().collect(),
        Feature::Draft {
            faces,
            neutral_plane,
            ..
        } => {
            let mut names: Vec<&EntityName> = faces.iter().collect();
            names.push(neutral_plane);
            names
        }
        Feature::Sweep {
            face, path_edge, ..
        } => vec![face, path_edge],
        Feature::Loft { sections, .. } => sections.iter().collect(),
        Feature::Pipe { path_edge, .. } => vec![path_edge],
        Feature::Revolve { face, .. } => vec![face],
        Feature::Reference { query, .. } => query.entity_names(),
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
/// originating [`FeatureId`].  The resulting list is deduplicated.
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

    deps
}

// ---------------------------------------------------------------------------
// Body resolution
// ---------------------------------------------------------------------------

/// Determine the [`BodyId`] for a newly created feature.
///
/// Rules:
/// - A [`Feature::Sketch`] with no dependencies creates a new body.
/// - All other features inherit the body of their first dependency.
/// - If there are no dependencies, fall back to `BodyId(0)`.
fn resolve_body(
    feature: &Feature,
    depends_on: &[FeatureId],
    body_manager: &mut BodyManager,
    features: &SlotMap<FeatureId, FeatureNode>,
) -> BodyId {
    // A base sketch (no planar face dependency) always starts a new body.
    if matches!(
        feature,
        Feature::Sketch {
            plane: SketchPlaneDesc::XY
                | SketchPlaneDesc::XZ
                | SketchPlaneDesc::YZ
                | SketchPlaneDesc::Custom { .. },
            ..
        }
    ) && depends_on.is_empty()
    {
        return body_manager.create_body();
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
fn register_feature_params(
    table: &mut ParameterTable,
    feature: &Feature,
    feature_id: FeatureId,
    label: &str,
) {
    let param_fields = feature_param_ids(feature);
    let prefix = label.to_lowercase().replace(' ', "_");
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

// ---------------------------------------------------------------------------
// FeatureTree methods
// ---------------------------------------------------------------------------

impl FeatureTree {
    /// Append a new feature at the end of the timeline.
    ///
    /// The feature's dependencies are automatically inferred from its
    /// [`EntityName`] references.  A new [`OperationId`] is assigned from the
    /// tree's internal counter.
    ///
    /// Body assignment follows these rules:
    /// - A [`Feature::Sketch`] with no dependencies creates a new body.
    /// - All other features (and sketches placed on a face) inherit the body
    ///   of their first dependency.  If there are no dependencies, they fall
    ///   back to `BodyId(0)`.
    pub fn add_feature(&mut self, feature: Feature, label: &str) -> FeatureId {
        let op_id = OperationId(self.next_op_id);
        self.next_op_id += 1;

        let depends_on = infer_dependencies(&feature, &self.op_index);

        let body = resolve_body(
            &feature,
            &depends_on,
            &mut self.body_manager,
            &self.features,
        );

        let node = FeatureNode {
            feature,
            operation_id: op_id,
            label: label.to_string(),
            depends_on,
            suppressed: false,
            body,
            suppress_if: None,
            last_result: None,
        };

        let feature_clone = node.feature.clone();
        let id = self.features.insert(node);
        self.timeline.push(id);
        self.body_manager.assign_feature(body, id);
        self.op_index.insert(op_id, id);

        // Register parameter scopes and qualified names.
        register_feature_params(&mut self.parameter_table, &feature_clone, id, label);

        let index = self.timeline.len() - 1;
        self.record_undo(TreeCommand::AddFeature {
            id,
            snapshot: self.features[id].clone(),
            index,
        });

        id
    }

    /// Return the timeline as a slice of [`FeatureId`]s in evaluation order.
    pub fn timeline(&self) -> &[FeatureId] {
        &self.timeline
    }

    /// Look up a feature node by its [`FeatureId`].
    pub fn feature(&self, id: FeatureId) -> Option<&FeatureNode> {
        self.features.get(id)
    }

    /// Insert a new feature at a specific position in the timeline.
    ///
    /// All auto-inferred dependencies must appear before `index` in the
    /// timeline.  Returns [`HistoryError::DependencyAfterInsert`] if any
    /// dependency would be at or after the insertion point, and
    /// [`HistoryError::InvalidPosition`] if `index` is out of range.
    pub fn insert_feature_at(
        &mut self,
        index: usize,
        feature: Feature,
        label: &str,
    ) -> Result<FeatureId, HistoryError> {
        if index > self.timeline.len() {
            return Err(HistoryError::InvalidPosition(index));
        }

        let op_id = OperationId(self.next_op_id);
        let depends_on = infer_dependencies(&feature, &self.op_index);

        // Validate: every dependency must appear before `index` in the timeline.
        let timeline_before = &self.timeline[..index];
        for &dep_id in &depends_on {
            if !timeline_before.contains(&dep_id) {
                let dep_label = self
                    .features
                    .get(dep_id)
                    .map(|n| n.label.clone())
                    .unwrap_or_else(|| String::from("<unknown>"));
                return Err(HistoryError::DependencyAfterInsert(dep_label));
            }
        }

        self.next_op_id += 1;

        let body = resolve_body(
            &feature,
            &depends_on,
            &mut self.body_manager,
            &self.features,
        );

        let node = FeatureNode {
            feature,
            operation_id: op_id,
            label: label.to_string(),
            depends_on,
            suppressed: false,
            body,
            suppress_if: None,
            last_result: None,
        };

        let feature_clone = node.feature.clone();
        let id = self.features.insert(node);
        self.timeline.insert(index, id);
        self.body_manager.assign_feature(body, id);
        self.op_index.insert(op_id, id);

        // Register parameter scopes and qualified names.
        register_feature_params(&mut self.parameter_table, &feature_clone, id, label);

        self.record_undo(TreeCommand::AddFeature {
            id,
            snapshot: self.features[id].clone(),
            index,
        });

        Ok(id)
    }

    /// Insert a new feature at the current rollback position.
    ///
    /// Equivalent to [`insert_feature_at`](Self::insert_feature_at) with
    /// the current [`rollback_position`](Self::rollback_position) (or the
    /// end of the timeline if no rollback bar is set).
    ///
    /// # Errors
    ///
    /// Same as [`insert_feature_at`](Self::insert_feature_at).
    pub fn insert_feature_at_rollback(
        &mut self,
        feature: Feature,
        label: &str,
    ) -> Result<FeatureId, HistoryError> {
        let index = self.rollback_index.unwrap_or(self.timeline.len());
        self.insert_feature_at(index, feature, label)
    }

    /// Move an existing feature to a new position in the timeline.
    ///
    /// Validates that:
    /// 1. All dependencies of the moved feature appear before `new_index`.
    /// 2. No feature that depends on the moved feature would end up after it.
    ///
    /// Returns [`HistoryError::FeatureNotFound`] if `id` is not in the tree,
    /// [`HistoryError::InvalidPosition`] if `new_index` is out of range, or
    /// [`HistoryError::DependencyAfterInsert`] if the move would violate
    /// dependency ordering.
    pub fn reorder(&mut self, id: FeatureId, new_index: usize) -> Result<(), HistoryError> {
        // Find and remove from current position.
        let current_pos = self
            .timeline
            .iter()
            .position(|&fid| fid == id)
            .ok_or(HistoryError::FeatureNotFound)?;

        if new_index >= self.timeline.len() {
            return Err(HistoryError::InvalidPosition(new_index));
        }

        // Build the candidate timeline with the feature removed.
        let mut candidate = self.timeline.clone();
        candidate.remove(current_pos);

        // Clamp new_index to the valid range after removal.
        let insert_at = new_index.min(candidate.len());

        // Check 1: all dependencies of `id` must appear before `insert_at`.
        if let Some(node) = self.features.get(id) {
            let before = &candidate[..insert_at];
            for dep_id in &node.depends_on {
                if !before.contains(dep_id) {
                    let dep_label = self
                        .features
                        .get(*dep_id)
                        .map(|n| n.label.clone())
                        .unwrap_or_else(|| String::from("<unknown>"));
                    return Err(HistoryError::DependencyAfterInsert(dep_label));
                }
            }
        }

        // Check 2: any feature that depends on `id` must not appear before `insert_at`.
        for &fid in &candidate[..insert_at] {
            if let Some(node) = self.features.get(fid)
                && node.depends_on.contains(&id)
            {
                let dep_label = node.label.clone();
                return Err(HistoryError::DependencyAfterInsert(dep_label));
            }
        }

        // Record the undo command before applying.
        // `insert_at` after removal mirrors the clamped logic below.
        let mut candidate_len = self.timeline.len();
        candidate_len -= 1; // account for the remove
        let final_index = new_index.min(candidate_len);

        self.record_undo(TreeCommand::ReorderFeature {
            id,
            old_index: current_pos,
            new_index: final_index,
        });

        // Apply the reorder.
        self.timeline.remove(current_pos);
        let insert_at = new_index.min(self.timeline.len());
        self.timeline.insert(insert_at, id);

        Ok(())
    }

    /// Suppress a feature, skipping it during replay.
    ///
    /// If the feature is already suppressed this is a no-op.
    pub fn suppress(&mut self, id: FeatureId) {
        if let Some(node) = self.features.get_mut(id)
            && !node.suppressed
        {
            node.suppressed = true;
            self.record_undo(TreeCommand::Suppress { id });
        }
    }

    /// Unsuppress a previously suppressed feature.
    ///
    /// If the feature is not suppressed this is a no-op.
    pub fn unsuppress(&mut self, id: FeatureId) {
        if let Some(node) = self.features.get_mut(id)
            && node.suppressed
        {
            node.suppressed = false;
            self.record_undo(TreeCommand::Unsuppress { id });
        }
    }

    /// Modify a feature's parameters.
    ///
    /// Records an undo command, re-infers dependencies, and invalidates the
    /// cache from this feature onward.
    pub fn modify_feature(&mut self, id: FeatureId, feature: Feature) {
        // Check existence first; bail early if not found.
        if !self.features.contains_key(id) {
            return;
        }

        let old = self.features[id].feature.clone();
        self.features[id].feature = feature.clone();

        // Re-infer dependencies using the op_index for O(1) lookup.
        let new_deps = infer_dependencies(&self.features[id].feature, &self.op_index);
        self.features[id].depends_on = new_deps;

        self.record_undo(TreeCommand::ModifyFeature {
            id,
            old,
            new: feature,
        });

        // Invalidate cache from this feature onward.
        if let Some(pos) = self.timeline.iter().position(|&fid| fid == id) {
            let timeline_snap = self.timeline.clone();
            self.cache.invalidate_from(pos, &timeline_snap);
        }
    }

    /// Rename a feature.
    ///
    /// Records an undo command.  Also updates the qualified names of all
    /// parameters scoped to this feature so that `"old_label.field"` becomes
    /// `"new_label.field"`.
    pub fn rename_feature(&mut self, id: FeatureId, label: &str) {
        if let Some(node) = self.features.get_mut(id) {
            let old_label = node.label.clone();
            node.label = label.to_string();
            // Update parameter qualified names for this feature.
            self.parameter_table
                .rename_feature_params(id, &old_label, label);
            self.record_undo(TreeCommand::RenameFeature {
                id,
                old_label,
                new_label: label.to_string(),
            });
        }
    }
}

#[cfg(test)]
#[path = "timeline_tests.rs"]
mod tests;
