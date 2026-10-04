//! Undo/redo command log for the parametric feature tree.
//!
//! [`TreeCommand`] represents a single reversible tree mutation.
//! [`reverse_command`] and [`apply_command`] directly mutate the
//! [`crate::FeatureTree`]'s internal state **without** recording new undo
//! entries — avoiding infinite recursion.
//!
//! The [`TreeCommand`] enum is defined in [`crate::undo_command`]; config
//! command handlers live in [`crate::undo_config`].

use crate::FeatureTree;
use crate::undo_config;

pub(crate) use crate::undo_command::TreeCommand;

// ---------------------------------------------------------------------------
// reverse_command (undo)
// ---------------------------------------------------------------------------

/// Apply the reverse of a command (undo).
///
/// Directly mutates `tree` without recording a new undo entry.
pub(crate) fn reverse_command(tree: &mut FeatureTree, cmd: &TreeCommand) {
    match cmd {
        TreeCommand::Compound(cmds) => {
            // Undo in reverse order.
            for c in cmds.iter().rev() {
                reverse_command(tree, c);
            }
        }
        TreeCommand::AddFeature { id, .. } => {
            // Remove from timeline only; leave the slotmap entry intact so that
            // redo can re-insert it without needing to generate a new key.
            if let Some(pos) = tree.timeline.iter().position(|&fid| fid == *id) {
                tree.timeline.remove(pos);
            }
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

        TreeCommand::RenameFeature {
            id,
            old_label,
            new_label,
        } => {
            tree.relabel_feature(*id, new_label, old_label);
        }

        TreeCommand::SetEntityLabel { label, old, .. } => {
            tree.labels.bind_entity(label, old.clone());
        }

        TreeCommand::SetBodyLabel { label, old, .. } => match old {
            Some(body) => {
                tree.labels.bodies.insert(label.clone(), *body);
            }
            None => {
                tree.labels.bodies.remove(label);
            }
        },

        TreeCommand::SetBuildTolerance { old, .. } => {
            tree.build_tolerance = *old;
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

        // Config commands — delegated to undo_config helpers.
        TreeCommand::AddConfig { id, .. } | TreeCommand::DuplicateConfig { id, .. } => {
            undo_config::reverse_add_config(tree, *id);
        }
        TreeCommand::DeleteConfig {
            id,
            snapshot,
            reparented,
        } => {
            undo_config::reverse_delete_config(tree, *id, snapshot, reparented);
        }
        TreeCommand::RenameConfig { id, old_name, .. } => {
            undo_config::reverse_rename_config(tree, *id, old_name);
        }
        TreeCommand::SetConfigDescription { id, old, .. } => {
            undo_config::reverse_set_config_description(tree, *id, old);
        }
        TreeCommand::SetConfigOverride {
            config, param, old, ..
        } => {
            undo_config::reverse_set_config_override(tree, *config, *param, old.as_ref());
        }
        TreeCommand::SetConfigSuppression {
            config,
            feature,
            old,
            ..
        } => {
            undo_config::reverse_set_config_suppression(tree, *config, *feature, *old);
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
        TreeCommand::Compound(cmds) => {
            // Redo in forward order.
            for c in cmds {
                apply_command(tree, c);
            }
        }
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

        TreeCommand::RenameFeature {
            id,
            old_label,
            new_label,
        } => {
            tree.relabel_feature(*id, old_label, new_label);
        }

        TreeCommand::SetEntityLabel { label, new, .. } => {
            tree.labels.bind_entity(label, new.clone());
        }

        TreeCommand::SetBodyLabel { label, new, .. } => match new {
            Some(body) => {
                tree.labels.bodies.insert(label.clone(), *body);
            }
            None => {
                tree.labels.bodies.remove(label);
            }
        },

        TreeCommand::SetBuildTolerance { new, .. } => {
            tree.build_tolerance = *new;
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

        // Config commands — delegated to undo_config helpers.
        TreeCommand::AddConfig { id, snapshot } | TreeCommand::DuplicateConfig { id, snapshot } => {
            undo_config::apply_add_config(tree, *id, snapshot);
        }
        TreeCommand::DeleteConfig { id, reparented, .. } => {
            undo_config::apply_delete_config(tree, *id, reparented);
        }
        TreeCommand::RenameConfig { id, new_name, .. } => {
            undo_config::apply_rename_config(tree, *id, new_name);
        }
        TreeCommand::SetConfigDescription { id, new, .. } => {
            undo_config::apply_set_config_description(tree, *id, new);
        }
        TreeCommand::SetConfigOverride {
            config, param, new, ..
        } => {
            undo_config::apply_set_config_override(tree, *config, *param, new.as_ref());
        }
        TreeCommand::SetConfigSuppression {
            config,
            feature,
            new,
            ..
        } => {
            undo_config::apply_set_config_suppression(tree, *config, *feature, *new);
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
fn move_in_timeline(
    tree: &mut FeatureTree,
    id: crate::node::FeatureId,
    _from_index: usize,
    to_index: usize,
) {
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
