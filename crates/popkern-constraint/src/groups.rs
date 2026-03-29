//! Constraint group data model.
//!
//! A [`ConstraintGroup`] is a named, enable/disable toggle for a set of
//! constraints. Disabling a group removes those constraints from the solver
//! without deleting them, so they can be re-enabled later.

use serde::{Deserialize, Serialize};

slotmap::new_key_type! {
    /// Opaque handle to a [`ConstraintGroup`] inside a [`Sketch`](crate::Sketch).
    pub struct GroupHandle;
}

/// A named group of constraints that can be enabled or disabled as a unit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConstraintGroup {
    /// Human-readable name for this group.
    pub name: String,
    /// When `false`, all constraints in this group are excluded from the
    /// solver.
    pub enabled: bool,
}
