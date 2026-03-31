//! Delete operations for the parametric feature tree.
//!
//! Provides [`DeleteBehavior`] to control what happens to features that depend
//! on the deleted feature, and [`DeleteResult`] reporting which features were
//! affected.

use serde::{Deserialize, Serialize};

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::node::FeatureId;
use crate::undo::TreeCommand;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Behavior when deleting a feature that has dependents.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum DeleteBehavior {
    /// Delete this feature and all transitively dependent features.
    Cascade,
    /// Reparent direct dependents to this feature's own dependencies.
    Reparent,
}

/// Result of a delete operation.
#[derive(Debug)]
pub struct DeleteResult {
    /// Features that were deleted (removed from the timeline).
    pub deleted: Vec<FeatureId>,
    /// Features that were reparented (their dependency list was rewritten).
    pub reparented: Vec<FeatureId>,
    /// Warnings about potential issues after the operation.
    pub warnings: Vec<String>,
}

// ---------------------------------------------------------------------------
// FeatureTree impl
// ---------------------------------------------------------------------------

impl FeatureTree {
    /// Delete a feature from the timeline, handling dependents according to
    /// `behavior`.
    ///
    /// Features are **not** removed from the internal [`SlotMap`] so that undo
    /// can restore them by re-inserting their [`FeatureId`]s into the timeline.
    /// Only the `timeline` vector (and, for `Cascade`, child entries) is
    /// modified.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::FeatureNotFound`] if `id` is not present in the
    /// tree.
    ///
    /// [`SlotMap`]: slotmap::SlotMap
    pub fn delete_feature(
        &mut self,
        id: FeatureId,
        behavior: DeleteBehavior,
    ) -> Result<DeleteResult, HistoryError> {
        if !self.features.contains_key(id) {
            return Err(HistoryError::FeatureNotFound);
        }

        match behavior {
            DeleteBehavior::Cascade => self.delete_cascade(id),
            DeleteBehavior::Reparent => self.delete_reparent(id),
        }
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Delete `id` and every feature that transitively depends on it.
    fn delete_cascade(&mut self, id: FeatureId) -> Result<DeleteResult, HistoryError> {
        // Collect all transitive dependents (including `id`), in timeline order.
        let to_delete = self.transitive_dependents(id);

        // Snapshot each deleted node together with its FeatureId and current
        // timeline position, so the undo handler can restore them at exactly
        // the right index.
        let mut snapshots: Vec<(crate::node::FeatureId, crate::node::FeatureNode, usize)> =
            Vec::new();
        for &did in &to_delete {
            if let Some(pos) = self.timeline.iter().position(|&fid| fid == did) {
                // clone here before we mutate the timeline below
                if let Some(node) = self.features.get(did) {
                    snapshots.push((did, node.clone(), pos));
                }
            }
        }

        // Collect param IDs BEFORE removing them (needed for config cleanup).
        let mut all_param_ids: Vec<crate::parameter::ParameterId> = Vec::new();
        for &did in &to_delete {
            all_param_ids.extend(self.parameter_table.feature_params(did));
        }

        // Remove from timeline (remove in reverse-index order to keep indices stable).
        for &did in to_delete.iter().rev() {
            self.timeline.retain(|&fid| fid != did);
        }

        // Remove parameters owned by all deleted features.
        for &did in &to_delete {
            self.parameter_table.remove_feature_params(did);
        }

        // Clean up config overrides that reference deleted features or params.
        for (_, config) in &mut self.configs {
            for &did in &to_delete {
                config.suppression.remove(&did);
            }
            for &pid in &all_param_ids {
                config.overrides.remove(&pid);
            }
        }

        self.record_undo(TreeCommand::DeleteFeature { snapshots });
        self.cache.clear();

        Ok(DeleteResult {
            deleted: to_delete,
            reparented: Vec::new(),
            warnings: Vec::new(),
        })
    }

    /// Delete `id` and reparent its direct dependents to `id`'s own parents.
    fn delete_reparent(&mut self, id: FeatureId) -> Result<DeleteResult, HistoryError> {
        // Direct dependents: features whose `depends_on` contains `id` and
        // that are present in the timeline.
        let direct_deps: Vec<FeatureId> = self
            .features
            .iter()
            .filter(|&(fid, node)| node.depends_on.contains(&id) && self.timeline.contains(&fid))
            .map(|(fid, _)| fid)
            .collect();

        // The deleted feature's own parents — dependents will be reparented here.
        let parent_deps = self
            .features
            .get(id)
            .map(|n| n.depends_on.clone())
            .unwrap_or_default();

        // Rewrite each direct dependent's dependency list.
        let mut warnings = Vec::new();
        for &dep_id in &direct_deps {
            if let Some(node) = self.features.get_mut(dep_id) {
                // Remove `id` from the list, then add the grandparent deps.
                node.depends_on.retain(|&did| did != id);
                for &pd in &parent_deps {
                    if !node.depends_on.contains(&pd) {
                        node.depends_on.push(pd);
                    }
                }
                // Preserve a stable order by deduplicating (order not
                // semantically significant, just needs to be unique).
                node.depends_on.dedup();
                warnings.push(format!(
                    "Feature '{}' reparented — references may not resolve",
                    node.label
                ));
            }
        }

        // Collect param IDs BEFORE removing them (needed for config cleanup).
        let param_ids: Vec<crate::parameter::ParameterId> = self.parameter_table.feature_params(id);

        // Snapshot the single deleted feature with its FeatureId.
        let pos = self.timeline.iter().position(|&fid| fid == id).unwrap_or(0);
        let snapshot = self
            .features
            .get(id)
            .cloned()
            .ok_or(HistoryError::FeatureNotFound)?;

        self.timeline.retain(|&fid| fid != id);

        // Remove parameters owned by the deleted feature.
        self.parameter_table.remove_feature_params(id);

        // Clean up config overrides that reference the deleted feature or params.
        for (_, config) in &mut self.configs {
            config.suppression.remove(&id);
            for &pid in &param_ids {
                config.overrides.remove(&pid);
            }
        }

        self.record_undo(TreeCommand::DeleteFeature {
            snapshots: vec![(id, snapshot, pos)],
        });
        self.cache.clear();

        Ok(DeleteResult {
            deleted: vec![id],
            reparented: direct_deps,
            warnings,
        })
    }

    /// Return `id` followed by all features that transitively depend on `id`,
    /// in timeline order.
    ///
    /// The returned [`Vec`] always begins with `id` itself (if it is in the
    /// timeline), followed by dependents in the order they appear in the
    /// timeline.
    pub(crate) fn transitive_dependents(&self, id: FeatureId) -> Vec<FeatureId> {
        use std::collections::HashSet;

        // Seed with the root.
        let mut to_delete: HashSet<FeatureId> = HashSet::new();
        to_delete.insert(id);

        // Walk the timeline in order, collecting anything that depends on an
        // already-marked feature.  A single pass is sufficient because the
        // timeline is in dependency order (dependencies always appear before
        // their dependents).
        for &fid in &self.timeline {
            if let Some(node) = self.features.get(fid)
                && node.depends_on.iter().any(|dep| to_delete.contains(dep))
            {
                to_delete.insert(fid);
            }
        }

        // Return in timeline order.
        self.timeline
            .iter()
            .filter(|&&fid| to_delete.contains(&fid))
            .copied()
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "delete_tests.rs"]
mod tests;
