//! [`Sketch`] methods for constraint group management.
//!
//! Groups allow a named set of constraints to be enabled or disabled together
//! without removing them from the sketch.

use crate::Sketch;
use crate::constraint::ConstraintHandle;
use crate::error::ConstraintError;
use crate::groups::{ConstraintGroup, GroupHandle};

#[cfg(test)]
#[path = "sketch_groups_tests.rs"]
mod tests;

impl Sketch {
    /// Create a new constraint group with the given name.
    ///
    /// The group starts enabled. Returns a handle that can be used to add
    /// constraints to the group.
    pub fn create_group(&mut self, name: &str) -> GroupHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        self.groups.insert(ConstraintGroup {
            name: name.to_owned(),
            enabled: true,
        })
    }

    /// Assign `constraint` to `group`.
    ///
    /// A constraint belongs to at most one group; any previous group
    /// membership is silently replaced. Returns [`ConstraintError::InvalidHandle`]
    /// if either handle is not present in the sketch.
    pub fn add_to_group(
        &mut self,
        group: GroupHandle,
        constraint: ConstraintHandle,
    ) -> Result<(), ConstraintError> {
        if !self.groups.contains_key(group) {
            return Err(ConstraintError::InvalidHandle);
        }
        let entry = self
            .constraints
            .get_mut(constraint)
            .ok_or(ConstraintError::InvalidHandle)?;
        entry.group = Some(group);
        Ok(())
    }

    /// Remove `constraint` from its current group, if any.
    ///
    /// Does nothing if the handle is invalid or the constraint has no group.
    pub fn remove_from_group(&mut self, constraint: ConstraintHandle) {
        if let Some(entry) = self.constraints.get_mut(constraint) {
            entry.group = None;
        }
    }

    /// Enable all constraints in `group`.
    ///
    /// Does nothing if the handle is invalid.
    pub fn enable_group(&mut self, group: GroupHandle) {
        if let Some(g) = self.groups.get_mut(group) {
            g.enabled = true;
        }
    }

    /// Disable all constraints in `group`.
    ///
    /// Disabled constraints are excluded from the solver until the group is
    /// re-enabled. Does nothing if the handle is invalid.
    pub fn disable_group(&mut self, group: GroupHandle) {
        if let Some(g) = self.groups.get_mut(group) {
            g.enabled = false;
        }
    }

    /// Delete `group` and clear its membership from all member constraints.
    ///
    /// After this call the constraints that were in the group still exist in
    /// the sketch and are active (they no longer have a group). Does nothing
    /// if the handle is invalid.
    pub fn delete_group(&mut self, group: GroupHandle) {
        if self.groups.remove(group).is_none() {
            return;
        }
        for entry in self.constraints.values_mut() {
            if entry.group == Some(group) {
                entry.group = None;
            }
        }
    }

    /// Return the enabled state of `group`, or `None` if the handle is invalid.
    pub fn is_group_enabled(&self, group: GroupHandle) -> Option<bool> {
        Some(self.groups.get(group)?.enabled)
    }
}
