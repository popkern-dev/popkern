//! [`FeatureTree`] mutation methods for the parametric feature timeline.
//!
//! Extracted from `timeline.rs` to keep each source file under 500 lines.
//! All items here are part of the `crate::timeline` module.

use popkern_core::naming::OperationId;

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::feature::Feature;
use crate::node::{FeatureId, FeatureNode};
use crate::undo::TreeCommand;

use super::{infer_dependencies, register_feature_params, resolve_body};

impl FeatureTree {
    /// Append a new feature at the end of the timeline.
    ///
    /// The feature's dependencies are automatically inferred from its
    /// [`crate::naming::EntityName`] references.  A new [`OperationId`] is
    /// assigned from the tree's internal counter.
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
