//! [`Sketch`] methods for named configuration management.
//!
//! Configurations allow the sketch to be toggled between distinct design
//! variants (e.g. "open" vs "closed") without duplicating geometry.

use crate::config::{
    get_constraint_value, set_constraint_value_on, ConfigHandle, Configuration, ConstraintSnapshot,
};
use crate::constraint::ConstraintKind;
use crate::error::{ConstraintError, SolveResult};
use crate::Sketch;

impl Sketch {
    /// Create a new, empty configuration with the given name.
    ///
    /// Call [`Sketch::save_config`] to populate it with the current sketch
    /// state.
    pub fn create_config(&mut self, name: &str) -> ConfigHandle {
        self.configs.insert(Configuration {
            name: name.to_owned(),
            constraint_states: std::collections::HashMap::new(),
            group_states: std::collections::HashMap::new(),
            solved_params: None,
        })
    }

    /// Overwrite `config` with a snapshot of the current sketch state.
    ///
    /// Captures the suppressed/driven state and dimensional value of every
    /// constraint, the enabled flag of every group, and the current parameter
    /// vector as a warm-start seed. Returns
    /// [`ConstraintError::InvalidHandle`] if `config` is not present in the
    /// sketch.
    pub fn save_config(&mut self, config: ConfigHandle) -> Result<(), ConstraintError> {
        let cfg = self
            .configs
            .get_mut(config)
            .ok_or(ConstraintError::InvalidHandle)?;

        cfg.constraint_states.clear();
        for (handle, entry) in self.constraints.iter() {
            let suppressed = entry.kind == ConstraintKind::Driven;
            let value = get_constraint_value(&entry.constraint);
            cfg.constraint_states
                .insert(handle, ConstraintSnapshot { suppressed, value });
        }

        cfg.group_states.clear();
        for (handle, group) in self.groups.iter() {
            cfg.group_states.insert(handle, group.enabled);
        }

        cfg.solved_params = Some(self.params.clone());

        Ok(())
    }

    /// Restore the sketch state from `config` and re-solve.
    ///
    /// For each constraint handle stored in the configuration:
    /// - If the handle still exists in the sketch, restores the suppressed
    ///   state and dimensional value.
    /// - Handles that no longer exist are silently skipped.
    ///
    /// Group enabled flags are restored similarly. The parameter vector is
    /// replaced by the saved warm-start before solving. Returns
    /// [`ConstraintError::InvalidHandle`] if `config` is not present.
    pub fn activate_config(
        &mut self,
        config: ConfigHandle,
    ) -> Result<SolveResult, ConstraintError> {
        // Clone out what we need to avoid borrow-checker issues while mutating
        // constraints/groups/params.
        let cfg = self
            .configs
            .get(config)
            .ok_or(ConstraintError::InvalidHandle)?
            .clone();

        // Restore constraint states.
        for (handle, snapshot) in &cfg.constraint_states {
            if let Some(entry) = self.constraints.get_mut(*handle) {
                entry.kind = if snapshot.suppressed {
                    ConstraintKind::Driven
                } else {
                    ConstraintKind::Driving
                };
                if let Some(v) = snapshot.value {
                    set_constraint_value_on(&mut entry.constraint, v);
                }
            }
        }

        // Restore group states.
        for (handle, enabled) in &cfg.group_states {
            if let Some(g) = self.groups.get_mut(*handle) {
                g.enabled = *enabled;
            }
        }

        // Warm-start: replace params with saved vector if shapes still match.
        if let Some(saved) = &cfg.solved_params {
            if saved.len() == self.params.len() {
                self.params.clone_from(saved);
            }
        }

        self.solve()
    }

    /// Remove `config` from the sketch.
    ///
    /// Does nothing if the handle is not present.
    pub fn delete_config(&mut self, config: ConfigHandle) {
        self.configs.remove(config);
    }

    /// Return the list of `(handle, name)` pairs for all configurations.
    pub fn list_configs(&self) -> Vec<(ConfigHandle, &str)> {
        self.configs
            .iter()
            .map(|(h, c)| (h, c.name.as_str()))
            .collect()
    }
}

#[cfg(test)]
#[path = "sketch_config_tests.rs"]
mod tests;
