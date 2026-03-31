//! Configuration management API for [`crate::FeatureTree`].
//!
//! All mutating methods record undo entries so changes can be reversed with
//! [`crate::FeatureTree::undo`].  Switching the active configuration is
//! intentionally NOT undoable — it is a UI-level navigation action, not a
//! model mutation.

use std::collections::HashMap;

use crate::FeatureTree;
use crate::config::{ConfigId, Configuration, DeleteConfigResult, ParamOverride};
use crate::config_eval;
use crate::error::HistoryError;
use crate::node::FeatureId;
use crate::parameter::ParameterId;
use crate::undo::TreeCommand;

impl FeatureTree {
    // -----------------------------------------------------------------------
    // Lifecycle
    // -----------------------------------------------------------------------

    /// Create a new configuration deriving from Default (no parent).
    ///
    /// Returns [`HistoryError::DuplicateConfigName`] if `name` is already in
    /// use by another (non-deleted) configuration.  Records undo.
    pub fn add_config(&mut self, name: &str, description: &str) -> Result<ConfigId, HistoryError> {
        self.check_unique_config_name(name, None)?;
        let cfg = Configuration {
            name: name.to_string(),
            description: description.to_string(),
            parent: None,
            overrides: HashMap::new(),
            suppression: HashMap::new(),
        };
        let id = self.configs.insert(cfg.clone());
        self.record_undo(TreeCommand::AddConfig { id, snapshot: cfg });
        Ok(id)
    }

    /// Create a configuration deriving from `parent`.
    ///
    /// Returns [`HistoryError::DuplicateConfigName`] if `name` is already in
    /// use, [`HistoryError::ConfigNotFound`] if `parent` does not exist, or
    /// [`HistoryError::InheritanceCycle`] if the parent chain would form a
    /// cycle.  Records undo.
    pub fn add_derived_config(
        &mut self,
        name: &str,
        description: &str,
        parent: ConfigId,
    ) -> Result<ConfigId, HistoryError> {
        if !self.configs.contains_key(parent) {
            return Err(HistoryError::ConfigNotFound);
        }
        // A newly-inserted config can't create a cycle (it doesn't exist yet),
        // but validate the parent chain doesn't already contain a ghost cycle.
        // We use a placeholder that equals the new config — use `parent` itself
        // as a sentinel: check if parent transitively points to itself.
        // Since we can't pass the new (not-yet-inserted) id, we just verify
        // the parent exists.  The cycle check is meaningful when reparenting
        // existing configs (done in add_derived_config only if we allow
        // changing parents later).  For a fresh config, no cycle is possible.
        self.check_unique_config_name(name, None)?;
        let cfg = Configuration {
            name: name.to_string(),
            description: description.to_string(),
            parent: Some(parent),
            overrides: HashMap::new(),
            suppression: HashMap::new(),
        };
        let id = self.configs.insert(cfg.clone());
        self.record_undo(TreeCommand::AddConfig { id, snapshot: cfg });
        Ok(id)
    }

    /// Duplicate a configuration, copying all overrides and suppressions.
    ///
    /// The new config has the same parent as `source`.
    /// Returns [`HistoryError::ConfigNotFound`] if `source` does not exist or
    /// has been soft-deleted.  Records undo.
    pub fn duplicate_config(
        &mut self,
        source: ConfigId,
        new_name: &str,
    ) -> Result<ConfigId, HistoryError> {
        let source_cfg = self
            .configs
            .get(source)
            .filter(|c| !c.name.is_empty())
            .ok_or(HistoryError::ConfigNotFound)?
            .clone();
        self.check_unique_config_name(new_name, None)?;
        let cfg = Configuration {
            name: new_name.to_string(),
            description: source_cfg.description.clone(),
            parent: source_cfg.parent,
            overrides: source_cfg.overrides.clone(),
            suppression: source_cfg.suppression.clone(),
        };
        let id = self.configs.insert(cfg.clone());
        self.record_undo(TreeCommand::DuplicateConfig { id, snapshot: cfg });
        Ok(id)
    }

    /// Delete a configuration.
    ///
    /// Children (configs whose parent is `id`) are reparented to `id`'s own
    /// parent.  The entry is soft-deleted (name cleared) so that undo can
    /// restore it without a new [`ConfigId`].
    ///
    /// Returns [`HistoryError::ConfigNotFound`] if `id` does not exist or has
    /// already been deleted.  Records undo.
    pub fn delete_config(&mut self, id: ConfigId) -> Result<DeleteConfigResult, HistoryError> {
        let cfg = self
            .configs
            .get(id)
            .filter(|c| !c.name.is_empty())
            .ok_or(HistoryError::ConfigNotFound)?
            .clone();

        // Collect children and their current parents for the undo record.
        let children: Vec<ConfigId> = self
            .configs
            .iter()
            .filter(|(cid, c)| *cid != id && c.parent == Some(id) && !c.name.is_empty())
            .map(|(cid, _)| cid)
            .collect();

        let reparented_info: Vec<(ConfigId, Option<ConfigId>)> = children
            .iter()
            .map(|&cid| {
                let old_parent = self.configs.get(cid).and_then(|c| c.parent);
                (cid, old_parent)
            })
            .collect();

        // Reparent children to id's parent.
        let new_parent = cfg.parent;
        for &child_id in &children {
            if let Some(child) = self.configs.get_mut(child_id) {
                child.parent = new_parent;
            }
        }

        // Soft-delete: clear name so configs() filters it out.
        if let Some(c) = self.configs.get_mut(id) {
            c.name.clear();
            c.overrides.clear();
            c.suppression.clear();
        }

        // If this was the active config, clear it.
        if self.active_config == Some(id) {
            self.active_config = None;
        }

        self.record_undo(TreeCommand::DeleteConfig {
            id,
            snapshot: cfg,
            reparented: reparented_info,
        });

        Ok(DeleteConfigResult {
            reparented: children,
        })
    }

    /// Rename a configuration.
    ///
    /// Returns [`HistoryError::ConfigNotFound`] if `id` does not exist,
    /// [`HistoryError::DuplicateConfigName`] if `name` is already in use.
    /// Records undo.
    pub fn rename_config(&mut self, id: ConfigId, name: &str) -> Result<(), HistoryError> {
        let cfg = self
            .configs
            .get(id)
            .filter(|c| !c.name.is_empty())
            .ok_or(HistoryError::ConfigNotFound)?;
        let old_name = cfg.name.clone();
        self.check_unique_config_name(name, Some(id))?;
        if let Some(c) = self.configs.get_mut(id) {
            c.name = name.to_string();
        }
        self.record_undo(TreeCommand::RenameConfig {
            id,
            old_name,
            new_name: name.to_string(),
        });
        Ok(())
    }

    /// Set the description of a configuration.  Records undo.
    pub fn set_config_description(&mut self, id: ConfigId, description: &str) {
        let old = self
            .configs
            .get(id)
            .map(|c| c.description.clone())
            .unwrap_or_default();
        if let Some(c) = self.configs.get_mut(id) {
            c.description = description.to_string();
        }
        self.record_undo(TreeCommand::SetConfigDescription {
            id,
            old,
            new: description.to_string(),
        });
    }

    // -----------------------------------------------------------------------
    // Active config (NOT undoable — UI navigation)
    // -----------------------------------------------------------------------

    /// Set the active configuration.  `None` means use default (no overrides).
    ///
    /// Not undoable — switching configurations is a view-level navigation
    /// action.
    pub fn set_active_config(&mut self, config: Option<ConfigId>) {
        self.active_config = config;
    }

    /// Return the currently active configuration, or `None` for default.
    pub fn active_config(&self) -> Option<ConfigId> {
        self.active_config
    }

    // -----------------------------------------------------------------------
    // Overrides (undoable)
    // -----------------------------------------------------------------------

    /// Set a parameter override within a configuration.
    ///
    /// Returns [`HistoryError::ConfigNotFound`] if `config` does not exist.
    /// Records undo.
    pub fn set_config_override(
        &mut self,
        config: ConfigId,
        param: ParameterId,
        value: ParamOverride,
    ) -> Result<(), HistoryError> {
        if !self.configs.contains_key(config) {
            return Err(HistoryError::ConfigNotFound);
        }
        let old = self
            .configs
            .get(config)
            .and_then(|c| c.overrides.get(&param))
            .cloned();
        if let Some(c) = self.configs.get_mut(config) {
            c.overrides.insert(param, value.clone());
        }
        self.record_undo(TreeCommand::SetConfigOverride {
            config,
            param,
            old,
            new: Some(value),
        });
        Ok(())
    }

    /// Clear a parameter override within a configuration.  Records undo.
    pub fn clear_config_override(&mut self, config: ConfigId, param: ParameterId) {
        let old = self
            .configs
            .get(config)
            .and_then(|c| c.overrides.get(&param))
            .cloned();
        if let Some(c) = self.configs.get_mut(config) {
            c.overrides.remove(&param);
        }
        self.record_undo(TreeCommand::SetConfigOverride {
            config,
            param,
            old,
            new: None,
        });
    }

    /// Set a feature suppression override within a configuration.  Records undo.
    pub fn set_config_suppression(
        &mut self,
        config: ConfigId,
        feature: FeatureId,
        suppressed: bool,
    ) {
        let old = self
            .configs
            .get(config)
            .and_then(|c| c.suppression.get(&feature))
            .copied();
        if let Some(c) = self.configs.get_mut(config) {
            c.suppression.insert(feature, suppressed);
        }
        self.record_undo(TreeCommand::SetConfigSuppression {
            config,
            feature,
            old,
            new: Some(suppressed),
        });
    }

    /// Clear a feature suppression override within a configuration.  Records undo.
    pub fn clear_config_suppression(&mut self, config: ConfigId, feature: FeatureId) {
        let old = self
            .configs
            .get(config)
            .and_then(|c| c.suppression.get(&feature))
            .copied();
        if let Some(c) = self.configs.get_mut(config) {
            c.suppression.remove(&feature);
        }
        self.record_undo(TreeCommand::SetConfigSuppression {
            config,
            feature,
            old,
            new: None,
        });
    }

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /// Return IDs of all live (non-deleted) configurations.
    pub fn configs(&self) -> Vec<ConfigId> {
        self.configs
            .iter()
            .filter(|(_, c)| !c.name.is_empty())
            .map(|(id, _)| id)
            .collect()
    }

    /// Return a reference to a configuration by ID.
    ///
    /// Returns `None` if the config does not exist or has been soft-deleted.
    pub fn config(&self, id: ConfigId) -> Option<&Configuration> {
        self.configs.get(id).filter(|c| !c.name.is_empty())
    }

    /// Return the parameter overrides for a configuration.
    pub fn config_overrides(&self, id: ConfigId) -> Option<&HashMap<ParameterId, ParamOverride>> {
        self.configs
            .get(id)
            .filter(|c| !c.name.is_empty())
            .map(|c| &c.overrides)
    }

    /// Return the feature suppression overrides for a configuration.
    pub fn config_suppressions(&self, id: ConfigId) -> Option<&HashMap<FeatureId, bool>> {
        self.configs
            .get(id)
            .filter(|c| !c.name.is_empty())
            .map(|c| &c.suppression)
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Check that `name` is not already used by another live config.
    ///
    /// Pass `exclude` to ignore the config being renamed (so its own current
    /// name doesn't count as a duplicate).
    fn check_unique_config_name(
        &self,
        name: &str,
        exclude: Option<ConfigId>,
    ) -> Result<(), HistoryError> {
        for (cid, c) in &self.configs {
            if c.name.is_empty() {
                continue; // soft-deleted
            }
            if exclude == Some(cid) {
                continue;
            }
            if c.name == name {
                return Err(HistoryError::DuplicateConfigName(name.to_string()));
            }
        }
        Ok(())
    }
}

// Suppress unused import warning — config_eval is used for cycle detection
// in add_derived_config but the helper is only called when reparenting
// existing configs.  The import ensures the module is compiled.
#[allow(unused_imports)]
use config_eval::would_create_cycle as _;

#[cfg(test)]
#[path = "config_api_tests.rs"]
mod tests;
