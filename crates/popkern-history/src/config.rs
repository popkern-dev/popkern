//! Configuration types for named design-table variants.
//!
//! A [`Configuration`] is a named variant of the model that overrides
//! specific parameter values and feature suppression states.  Configurations
//! form an inheritance tree: each configuration inherits from a parent (or
//! from the default state if `parent` is `None`), and only the *delta* from
//! the parent is stored.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use slotmap::new_key_type;

use crate::expression::Expression;
use crate::node::FeatureId;
use crate::parameter::ParameterId;

new_key_type! {
    /// Stable handle for a named configuration in the [`crate::FeatureTree`].
    pub struct ConfigId;
}

/// A parameter value override within a configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ParamOverride {
    /// Pin to a literal value (SI units). The default expression is ignored
    /// in this configuration.
    Literal(f64),
    /// Use a different expression in this configuration, replacing the
    /// default expression.
    Expression(Expression),
}

/// A named configuration — a variant of the model with specific
/// parameter and suppression overrides.
///
/// Configurations inherit from a parent (or from Default if `parent` is None).
/// Only the overrides stored here differ from the parent; all other values
/// are inherited.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Configuration {
    /// Display name (e.g., "Large", "Small-Steel").
    pub name: String,
    /// Human-readable description.
    pub description: String,
    /// Parent configuration for inheritance. `None` = derives from Default.
    pub parent: Option<ConfigId>,
    /// Parameter value/expression overrides. Only parameters that differ
    /// from the parent are stored here.
    pub overrides: HashMap<ParameterId, ParamOverride>,
    /// Feature suppression overrides. Only features whose suppression
    /// state differs from the parent are stored here.
    pub suppression: HashMap<FeatureId, bool>,
}

/// Result of deleting a configuration.
#[derive(Debug)]
pub struct DeleteConfigResult {
    /// Derived configurations that were reparented to the deleted
    /// config's parent.
    pub reparented: Vec<ConfigId>,
}

#[cfg(test)]
#[path = "config_tests.rs"]
mod tests;
