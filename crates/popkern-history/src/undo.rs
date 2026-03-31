//! Undo/redo command log for the parametric feature tree.
//!
//! [`TreeCommand`] represents a single reversible tree mutation.
//! [`reverse_command`] and [`apply_command`] directly mutate the
//! [`crate::FeatureTree`]'s internal state **without** recording new undo
//! entries — avoiding infinite recursion.

use crate::FeatureTree;
use crate::expression::Expression;
use crate::feature::Feature;
use crate::node::{FeatureId, FeatureNode};
use crate::parameter::{Parameter, ParameterId};
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// TreeCommand
// ---------------------------------------------------------------------------

/// A single reversible tree mutation recorded in the undo log.
#[derive(Debug, Clone)]
#[allow(dead_code)] // DeleteFeature wired in Task 7; snapshot kept for spec completeness
pub(crate) enum TreeCommand {
    /// A feature was added at a specific timeline index.
    AddFeature {
        /// The ID that was assigned to the new feature.
        id: FeatureId,
        /// Full snapshot of the node at the time of insertion.
        snapshot: FeatureNode,
        /// The timeline index at which the feature was inserted.
        index: usize,
    },
    /// One or more features were deleted.
    DeleteFeature {
        /// Each deleted feature's ID, node snapshot, and former timeline index.
        snapshots: Vec<(FeatureId, FeatureNode, usize)>,
    },
    /// A feature's parameters were modified.
    ModifyFeature {
        /// The ID of the modified feature.
        id: FeatureId,
        /// The previous feature parameters.
        old: Feature,
        /// The new feature parameters.
        new: Feature,
    },
    /// A feature was moved in the timeline.
    ReorderFeature {
        /// The ID of the moved feature.
        id: FeatureId,
        /// The timeline index before the move.
        old_index: usize,
        /// The timeline index after the move.
        new_index: usize,
    },
    /// A feature was suppressed.
    Suppress {
        /// The ID of the suppressed feature.
        id: FeatureId,
    },
    /// A feature was unsuppressed.
    Unsuppress {
        /// The ID of the unsuppressed feature.
        id: FeatureId,
    },
    /// The rollback bar was moved.
    MoveRollback {
        /// The previous rollback position.
        old: Option<usize>,
        /// The new rollback position.
        new: Option<usize>,
    },
    /// A feature was renamed.
    RenameFeature {
        /// The ID of the renamed feature.
        id: FeatureId,
        /// The label before the rename.
        old_label: String,
        /// The label after the rename.
        new_label: String,
    },
    /// An expression was set or cleared on a parameter.
    SetExpression {
        /// The parameter whose expression changed.
        param: ParameterId,
        /// The expression before the change (`None` = was literal).
        old: Option<Expression>,
        /// The expression after the change (`None` = cleared).
        new: Option<Expression>,
    },
    /// A parameter's value and/or unit was changed.
    SetParamValue {
        /// The parameter whose value changed.
        param: ParameterId,
        /// The value before the change (SI units).
        old_value: f64,
        /// The unit before the change.
        old_unit: Unit,
        /// The value after the change (SI units).
        new_value: f64,
        /// The unit after the change.
        new_unit: Unit,
    },
    /// Parameter bounds were changed or cleared.
    SetBounds {
        /// The parameter whose bounds changed.
        param: ParameterId,
        /// The bounds before the change.
        old: Option<(f64, f64)>,
        /// The bounds after the change.
        new: Option<(f64, f64)>,
    },
    /// A global parameter was added.
    AddGlobalParam {
        /// The stable ID assigned to the new parameter.
        id: ParameterId,
        /// Full snapshot of the parameter at insertion time.
        snapshot: Parameter,
    },
    /// A global parameter was removed.
    RemoveGlobalParam {
        /// The ID of the removed parameter.
        id: ParameterId,
        /// Full snapshot of the parameter before removal.
        snapshot: Parameter,
    },
    /// A parameter was renamed.
    RenameParam {
        /// The parameter that was renamed.
        id: ParameterId,
        /// The name before the rename.
        old_name: String,
        /// The name after the rename.
        new_name: String,
    },
    /// A suppress_if expression was set or cleared on a feature.
    SetSuppressIf {
        /// The feature whose suppress_if changed.
        feature: FeatureId,
        /// The expression before the change.
        old: Option<Expression>,
        /// The expression after the change.
        new: Option<Expression>,
    },
}

// ---------------------------------------------------------------------------
// reverse_command (undo)
// ---------------------------------------------------------------------------

/// Apply the reverse of a command (undo).
///
/// Directly mutates `tree` without recording a new undo entry.
pub(crate) fn reverse_command(tree: &mut FeatureTree, cmd: &TreeCommand) {
    match cmd {
        TreeCommand::AddFeature { id, index, .. } => {
            // Remove from timeline only; leave the slotmap entry intact so that
            // redo can re-insert it without needing to generate a new key.
            if let Some(pos) = tree.timeline.iter().position(|&fid| fid == *id) {
                tree.timeline.remove(pos);
            }
            // Invalidate any cached state that depended on this feature.
            tree.cache.invalidate_from(*index, &tree.timeline.clone());
        }

        TreeCommand::DeleteFeature { snapshots } => {
            // Re-insert each feature ID at its former timeline position.
            // The SlotMap entries are still intact (delete only removes from
            // the timeline), so we just restore the timeline ordering.
            let mut sorted = snapshots.clone();
            sorted.sort_by_key(|(_, _, idx)| *idx);
            for (id, node, idx) in sorted {
                let insert_at = idx.min(tree.timeline.len());
                tree.timeline.insert(insert_at, id);
                tree.body_manager.assign_feature(node.body, id);
            }
        }

        TreeCommand::ModifyFeature { id, old, .. } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.feature = old.clone();
            }
            // Invalidate cache from this feature onward.
            if let Some(pos) = tree.timeline.iter().position(|&fid| fid == *id) {
                let timeline_snap = tree.timeline.clone();
                tree.cache.invalidate_from(pos, &timeline_snap);
            }
        }

        TreeCommand::ReorderFeature {
            id,
            old_index,
            new_index,
        } => {
            // Move from new_index back to old_index.
            move_in_timeline(tree, *id, *new_index, *old_index);
        }

        TreeCommand::Suppress { id } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.suppressed = false;
            }
        }

        TreeCommand::Unsuppress { id } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.suppressed = true;
            }
        }

        TreeCommand::MoveRollback { old, .. } => {
            tree.rollback_index = *old;
        }

        TreeCommand::RenameFeature { id, old_label, .. } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.label = old_label.clone();
            }
        }

        TreeCommand::SetExpression { param, old, .. } => match old {
            Some(expr) => {
                let _ = tree.parameter_table.set_expression(*param, expr.clone());
            }
            None => {
                tree.parameter_table.clear_expression(*param);
            }
        },

        TreeCommand::SetParamValue {
            param,
            old_value,
            old_unit,
            ..
        } => {
            let _ = tree
                .parameter_table
                .set_value(*param, *old_value, *old_unit);
        }

        TreeCommand::SetBounds { param, old, .. } => match old {
            Some((min, max)) => tree.parameter_table.set_bounds(*param, *min, *max),
            None => tree.parameter_table.clear_bounds(*param),
        },

        TreeCommand::AddGlobalParam { id, .. } => {
            // On undo, remove from name_index only — leave the SlotMap entry
            // intact so that redo can re-register the name without a key
            // conflict.
            tree.parameter_table.deregister_name(*id);
        }

        TreeCommand::RemoveGlobalParam { id, snapshot } => {
            // On undo of a remove, re-register the name in the index.
            // The SlotMap entry was never actually removed (same pattern as
            // AddGlobalParam undo above).
            tree.parameter_table.reregister_name(*id, &snapshot.name);
        }

        TreeCommand::RenameParam { id, old_name, .. } => {
            let _ = tree.parameter_table.rename(*id, old_name);
        }

        TreeCommand::SetSuppressIf { feature, old, .. } => {
            if let Some(node) = tree.features.get_mut(*feature) {
                node.suppress_if = old.clone();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// apply_command (redo)
// ---------------------------------------------------------------------------

/// Re-apply a command (redo).
///
/// Directly mutates `tree` without recording a new undo entry.
pub(crate) fn apply_command(tree: &mut FeatureTree, cmd: &TreeCommand) {
    match cmd {
        TreeCommand::AddFeature { id, index, .. } => {
            // The slotmap entry still exists from the original insertion.
            // Simply re-insert the ID into the timeline at the recorded index.
            let insert_at = (*index).min(tree.timeline.len());
            tree.timeline.insert(insert_at, *id);
        }

        TreeCommand::DeleteFeature { snapshots } => {
            // Remove each feature from the timeline only — the SlotMap
            // entries are kept so undo can restore them.
            for (id, _, _) in snapshots {
                tree.timeline.retain(|&tid| tid != *id);
            }
        }

        TreeCommand::ModifyFeature { id, new, .. } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.feature = new.clone();
            }
            if let Some(pos) = tree.timeline.iter().position(|&fid| fid == *id) {
                let timeline_snap = tree.timeline.clone();
                tree.cache.invalidate_from(pos, &timeline_snap);
            }
        }

        TreeCommand::ReorderFeature {
            id,
            old_index,
            new_index,
        } => {
            // Move from old_index forward to new_index.
            move_in_timeline(tree, *id, *old_index, *new_index);
        }

        TreeCommand::Suppress { id } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.suppressed = true;
            }
        }

        TreeCommand::Unsuppress { id } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.suppressed = false;
            }
        }

        TreeCommand::MoveRollback { new, .. } => {
            tree.rollback_index = *new;
        }

        TreeCommand::RenameFeature { id, new_label, .. } => {
            if let Some(node) = tree.features.get_mut(*id) {
                node.label = new_label.clone();
            }
        }

        TreeCommand::SetExpression { param, new, .. } => match new {
            Some(expr) => {
                let _ = tree.parameter_table.set_expression(*param, expr.clone());
            }
            None => {
                tree.parameter_table.clear_expression(*param);
            }
        },

        TreeCommand::SetParamValue {
            param,
            new_value,
            new_unit,
            ..
        } => {
            let _ = tree
                .parameter_table
                .set_value(*param, *new_value, *new_unit);
        }

        TreeCommand::SetBounds { param, new, .. } => match new {
            Some((min, max)) => tree.parameter_table.set_bounds(*param, *min, *max),
            None => tree.parameter_table.clear_bounds(*param),
        },

        TreeCommand::AddGlobalParam { id, snapshot } => {
            // Redo: re-register the name in the index.
            tree.parameter_table.reregister_name(*id, &snapshot.name);
        }

        TreeCommand::RemoveGlobalParam { id, .. } => {
            // Redo: remove from name_index only, leave SlotMap intact.
            tree.parameter_table.deregister_name(*id);
        }

        TreeCommand::RenameParam { id, new_name, .. } => {
            let _ = tree.parameter_table.rename(*id, new_name);
        }

        TreeCommand::SetSuppressIf { feature, new, .. } => {
            if let Some(node) = tree.features.get_mut(*feature) {
                node.suppress_if = new.clone();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// FeatureTree public undo/redo API
// ---------------------------------------------------------------------------

impl FeatureTree {
    /// Undo the last tree mutation. Returns `true` if something was undone.
    pub fn undo(&mut self) -> bool {
        if self.undo_cursor == 0 {
            return false;
        }
        self.undo_cursor -= 1;
        let cmd = self.undo_log[self.undo_cursor].clone();
        reverse_command(self, &cmd);
        true
    }

    /// Redo the last undone mutation. Returns `true` if something was redone.
    pub fn redo(&mut self) -> bool {
        if self.undo_cursor >= self.undo_log.len() {
            return false;
        }
        let cmd = self.undo_log[self.undo_cursor].clone();
        apply_command(self, &cmd);
        self.undo_cursor += 1;
        true
    }

    /// Whether undo is available.
    pub fn can_undo(&self) -> bool {
        self.undo_cursor > 0
    }

    /// Whether redo is available.
    pub fn can_redo(&self) -> bool {
        self.undo_cursor < self.undo_log.len()
    }

    /// Record a command in the undo log. Truncates any redo history.
    pub(crate) fn record_undo(&mut self, cmd: TreeCommand) {
        self.undo_log.truncate(self.undo_cursor);
        self.undo_log.push(cmd);
        self.undo_cursor = self.undo_log.len();
    }
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Move feature `id` from `from_index` to `to_index` in the timeline, without
/// any dependency validation and without recording undo commands.
///
/// Uses a remove-then-insert approach, clamping `to_index` to the valid range
/// after removal.
fn move_in_timeline(tree: &mut FeatureTree, id: FeatureId, _from_index: usize, to_index: usize) {
    // Find the feature by ID rather than trusting the caller's index
    // directly — positions may have shifted.
    let actual_pos = match tree.timeline.iter().position(|&fid| fid == id) {
        Some(p) => p,
        None => return,
    };
    tree.timeline.remove(actual_pos);
    let insert_at = to_index.min(tree.timeline.len());
    tree.timeline.insert(insert_at, id);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "undo_tests.rs"]
mod tests;
