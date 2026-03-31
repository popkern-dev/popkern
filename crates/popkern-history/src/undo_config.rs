//! Config-related undo/redo helpers for [`crate::undo`].
//!
//! These functions handle the reverse and apply logic for all
//! [`crate::undo::TreeCommand`] variants that operate on configurations.
//! They are called from [`crate::undo::reverse_command`] and
//! [`crate::undo::apply_command`] to keep that file within the 500-line limit.

use crate::FeatureTree;
use crate::config::{ConfigId, Configuration, ParamOverride};
use crate::node::FeatureId;
use crate::parameter::ParameterId;

/// Undo handler for `AddConfig` / `DuplicateConfig`.
///
/// Soft-deletes the config by clearing its name so [`FeatureTree::configs`]
/// filters it out.
pub(crate) fn reverse_add_config(tree: &mut FeatureTree, id: ConfigId) {
    if let Some(c) = tree.configs.get_mut(id) {
        c.name.clear();
    }
}

/// Undo handler for `DeleteConfig`.
///
/// Restores the config from its snapshot and re-applies the original parent
/// references for any children that were reparented during deletion.
pub(crate) fn reverse_delete_config(
    tree: &mut FeatureTree,
    id: ConfigId,
    snapshot: &Configuration,
    reparented: &[(ConfigId, Option<ConfigId>)],
) {
    if let Some(c) = tree.configs.get_mut(id) {
        *c = snapshot.clone();
    }
    for (child_id, old_parent) in reparented {
        if let Some(child) = tree.configs.get_mut(*child_id) {
            child.parent = *old_parent;
        }
    }
}

/// Undo handler for `RenameConfig`.
pub(crate) fn reverse_rename_config(tree: &mut FeatureTree, id: ConfigId, old_name: &str) {
    if let Some(c) = tree.configs.get_mut(id) {
        c.name = old_name.to_string();
    }
}

/// Undo handler for `SetConfigDescription`.
pub(crate) fn reverse_set_config_description(tree: &mut FeatureTree, id: ConfigId, old: &str) {
    if let Some(c) = tree.configs.get_mut(id) {
        c.description = old.to_string();
    }
}

/// Undo handler for `SetConfigOverride`.
pub(crate) fn reverse_set_config_override(
    tree: &mut FeatureTree,
    config: ConfigId,
    param: ParameterId,
    old: Option<&ParamOverride>,
) {
    if let Some(c) = tree.configs.get_mut(config) {
        match old {
            Some(ov) => {
                c.overrides.insert(param, ov.clone());
            }
            None => {
                c.overrides.remove(&param);
            }
        }
    }
}

/// Undo handler for `SetConfigSuppression`.
pub(crate) fn reverse_set_config_suppression(
    tree: &mut FeatureTree,
    config: ConfigId,
    feature: FeatureId,
    old: Option<bool>,
) {
    if let Some(c) = tree.configs.get_mut(config) {
        match old {
            Some(v) => {
                c.suppression.insert(feature, v);
            }
            None => {
                c.suppression.remove(&feature);
            }
        }
    }
}

/// Redo handler for `AddConfig` / `DuplicateConfig`.
///
/// Restores the config snapshot (re-activates the soft-deleted entry).
pub(crate) fn apply_add_config(tree: &mut FeatureTree, id: ConfigId, snapshot: &Configuration) {
    if let Some(c) = tree.configs.get_mut(id) {
        *c = snapshot.clone();
    }
}

/// Redo handler for `DeleteConfig`.
///
/// Re-applies the delete by reparenting children and soft-deleting the entry.
pub(crate) fn apply_delete_config(
    tree: &mut FeatureTree,
    id: ConfigId,
    reparented: &[(ConfigId, Option<ConfigId>)],
) {
    let new_parent = tree.configs.get(id).and_then(|c| c.parent);
    for (child_id, _) in reparented {
        if let Some(child) = tree.configs.get_mut(*child_id) {
            child.parent = new_parent;
        }
    }
    if let Some(c) = tree.configs.get_mut(id) {
        c.name.clear();
        c.overrides.clear();
        c.suppression.clear();
    }
}

/// Redo handler for `RenameConfig`.
pub(crate) fn apply_rename_config(tree: &mut FeatureTree, id: ConfigId, new_name: &str) {
    if let Some(c) = tree.configs.get_mut(id) {
        c.name = new_name.to_string();
    }
}

/// Redo handler for `SetConfigDescription`.
pub(crate) fn apply_set_config_description(tree: &mut FeatureTree, id: ConfigId, new: &str) {
    if let Some(c) = tree.configs.get_mut(id) {
        c.description = new.to_string();
    }
}

/// Redo handler for `SetConfigOverride`.
pub(crate) fn apply_set_config_override(
    tree: &mut FeatureTree,
    config: ConfigId,
    param: ParameterId,
    new: Option<&ParamOverride>,
) {
    if let Some(c) = tree.configs.get_mut(config) {
        match new {
            Some(ov) => {
                c.overrides.insert(param, ov.clone());
            }
            None => {
                c.overrides.remove(&param);
            }
        }
    }
}

/// Redo handler for `SetConfigSuppression`.
pub(crate) fn apply_set_config_suppression(
    tree: &mut FeatureTree,
    config: ConfigId,
    feature: FeatureId,
    new: Option<bool>,
) {
    if let Some(c) = tree.configs.get_mut(config) {
        match new {
            Some(v) => {
                c.suppression.insert(feature, v);
            }
            None => {
                c.suppression.remove(&feature);
            }
        }
    }
}
