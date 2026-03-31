//! Parameter table for the parametric feature tree.
//!
//! A [`ParameterTable`] stores all named, typed parameters used by features.
//! Each parameter has a [`Unit`] for display conversion and a `value` in SI
//! base units.  Built-in read-only constants (`pi`, `e`, `tau`) are inserted
//! at construction time.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use slotmap::{SlotMap, new_key_type};

use crate::error::HistoryError;
use crate::expression::{Expression, expression_param_refs};
use crate::node::FeatureId;
use crate::unit::Unit;

new_key_type! {
    /// Stable handle for a parameter in the [`ParameterTable`].
    pub struct ParameterId;
}

/// Scope of a parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParamScope {
    /// User-defined global parameter.
    Global,
    /// Auto-generated parameter belonging to a feature.
    Feature(FeatureId),
}

/// Warning from parameter evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EvaluationWarning {
    /// A math error occurred (division by zero, sqrt of negative, etc.).
    MathError {
        /// Description of the math error.
        message: String,
    },
    /// The value was clamped to bounds.
    BoundsClamped {
        /// The original (unclamped) value.
        original: f64,
        /// The clamped value that was actually stored.
        clamped: f64,
    },
}

/// A single parameter in the parameter table.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Parameter {
    /// Short name (e.g., `"distance"`, `"width"`).
    pub name: String,
    /// Scope: global or owned by a feature.
    pub scope: ParamScope,
    /// Current resolved value in SI base units.
    pub value: f64,
    /// Display/input unit.
    pub unit: Unit,
    /// Expression driving this parameter. `None` = user-entered literal.
    pub expression: Option<Expression>,
    /// Optional min/max bounds in SI units.
    pub bounds: Option<(f64, f64)>,
    /// Whether this is a built-in read-only parameter (`pi`, `e`, `tau`).
    pub read_only: bool,
    /// Warning from the most recent evaluation.
    #[serde(skip)]
    pub eval_warning: Option<EvaluationWarning>,
    /// Last successfully computed value (for MEASURE staleness recovery).
    #[serde(default)]
    pub last_known: Option<f64>,
}

/// A table of all named parameters used by the feature tree.
///
/// Parameters are keyed by stable [`ParameterId`] handles.  A secondary
/// `name_index` allows O(1) lookup by qualified name.  Global parameters are
/// indexed by their bare name; feature-scoped parameters are indexed by a
/// caller-supplied qualified name such as `"extrude1.distance"`.
///
/// Built-in constants (`pi`, `e`, `tau`) are added at construction time with
/// `read_only = true`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterTable {
    params: SlotMap<ParameterId, Parameter>,
    name_index: HashMap<String, ParameterId>,
}

impl Default for ParameterTable {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterTable {
    /// Create a new table with built-in constants (`pi`, `e`, `tau`).
    pub fn new() -> Self {
        let mut table = Self {
            params: SlotMap::with_key(),
            name_index: HashMap::new(),
        };
        // Insert built-in read-only constants.
        let builtins = [
            ("pi", std::f64::consts::PI),
            ("e", std::f64::consts::E),
            ("tau", std::f64::consts::TAU),
        ];
        for (name, value) in builtins {
            let id = table.params.insert(Parameter {
                name: name.to_string(),
                scope: ParamScope::Global,
                value,
                unit: Unit::Ratio,
                expression: None,
                bounds: None,
                read_only: true,
                eval_warning: None,
                last_known: None,
            });
            table.name_index.insert(name.to_string(), id);
        }
        table
    }

    /// Add a user-defined global parameter.
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if the name is already
    /// in use.
    pub fn add_global(
        &mut self,
        name: &str,
        value: f64,
        unit: Unit,
    ) -> Result<ParameterId, HistoryError> {
        if self.name_index.contains_key(name) {
            return Err(HistoryError::DuplicateParameterName(name.to_string()));
        }
        let id = self.params.insert(Parameter {
            name: name.to_string(),
            scope: ParamScope::Global,
            value,
            unit,
            expression: None,
            bounds: None,
            read_only: false,
            eval_warning: None,
            last_known: None,
        });
        self.name_index.insert(name.to_string(), id);
        Ok(id)
    }

    /// Add a global parameter using a value in the given display unit.
    ///
    /// The value is automatically converted to SI for internal storage.
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if the name is already
    /// in use.
    pub fn add_global_display(
        &mut self,
        name: &str,
        display_value: f64,
        unit: Unit,
    ) -> Result<ParameterId, HistoryError> {
        self.add_global(name, unit.to_si(display_value), unit)
    }

    /// Add a feature-scoped parameter. Called internally during feature
    /// registration.
    ///
    /// The `name` argument is the qualified name (e.g., `"extrude1.distance"`)
    /// and is used as the index key.  The short name stored on the parameter
    /// is everything after the last `.`, or the full string if no `.` is
    /// present.
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if the qualified name
    /// is already in use.
    #[allow(dead_code)]
    pub(crate) fn add_feature_param(
        &mut self,
        name: &str,
        feature: FeatureId,
        value: f64,
        unit: Unit,
    ) -> Result<ParameterId, HistoryError> {
        if self.name_index.contains_key(name) {
            return Err(HistoryError::DuplicateParameterName(name.to_string()));
        }
        let short_name = name.rsplit('.').next().unwrap_or(name).to_string();
        let id = self.params.insert(Parameter {
            name: short_name,
            scope: ParamScope::Feature(feature),
            value,
            unit,
            expression: None,
            bounds: None,
            read_only: false,
            eval_warning: None,
            last_known: None,
        });
        self.name_index.insert(name.to_string(), id);
        Ok(id)
    }

    /// Get a parameter by ID.
    pub fn get(&self, id: ParameterId) -> Option<&Parameter> {
        self.params.get(id)
    }

    /// Get a mutable reference to a parameter (crate-internal only).
    #[allow(dead_code)]
    pub(crate) fn get_mut(&mut self, id: ParameterId) -> Option<&mut Parameter> {
        self.params.get_mut(id)
    }

    /// Resolve a qualified name to a [`ParameterId`].
    ///
    /// Returns `None` if no parameter with that name exists.
    pub fn resolve_name(&self, qualified_name: &str) -> Option<ParameterId> {
        self.name_index.get(qualified_name).copied()
    }

    /// Rename a parameter. Updates the name index.
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if `new_name` is
    /// already in use, or [`HistoryError::FeatureNotFound`] if `id` is invalid.
    pub fn rename(&mut self, id: ParameterId, new_name: &str) -> Result<(), HistoryError> {
        if self.name_index.contains_key(new_name) {
            return Err(HistoryError::DuplicateParameterName(new_name.to_string()));
        }
        // Find the current indexed name for this id.
        let old_key = self
            .name_index
            .iter()
            .find(|(_, v)| **v == id)
            .map(|(k, _)| k.clone())
            .ok_or(HistoryError::FeatureNotFound)?;

        self.name_index.remove(&old_key);
        self.name_index.insert(new_name.to_string(), id);

        let param = self
            .params
            .get_mut(id)
            .ok_or(HistoryError::FeatureNotFound)?;
        param.name = new_name.rsplit('.').next().unwrap_or(new_name).to_string();
        Ok(())
    }

    /// Remove a parameter.
    ///
    /// Returns [`HistoryError::ParameterInUse`] if other parameters reference
    /// it (currently always succeeds since expressions are not yet implemented).
    /// Returns [`HistoryError::ReadOnlyParameter`] if the parameter is
    /// read-only.
    pub fn remove(&mut self, id: ParameterId) -> Result<(), HistoryError> {
        let param = self.params.get(id).ok_or(HistoryError::FeatureNotFound)?;
        if param.read_only {
            return Err(HistoryError::ReadOnlyParameter);
        }
        let dependents = self.dependents(id);
        if !dependents.is_empty() {
            return Err(HistoryError::ParameterInUse(dependents));
        }
        // Remove from name_index.
        self.name_index.retain(|_, &mut v| v != id);
        self.params.remove(id);
        Ok(())
    }

    /// Remove a parameter and clear all expressions referencing it.
    ///
    /// Any parameter whose expression references `id` has its expression
    /// cleared (keeping its last resolved value). Silently ignores invalid
    /// `id`.
    pub fn remove_cascade(&mut self, id: ParameterId) {
        if self.params.get(id).is_none() {
            return;
        }
        // Clear expressions that reference this parameter.
        let dependents = self.dependents(id);
        for dep_id in dependents {
            if let Some(param) = self.params.get_mut(dep_id) {
                param.expression = None;
            }
        }
        self.name_index.retain(|_, &mut v| v != id);
        self.params.remove(id);
    }

    /// Remove all parameters scoped to a feature.
    pub fn remove_feature_params(&mut self, feature: FeatureId) {
        let to_remove: Vec<ParameterId> = self
            .params
            .iter()
            .filter(|(_, p)| p.scope == ParamScope::Feature(feature))
            .map(|(id, _)| id)
            .collect();
        for id in to_remove {
            self.name_index.retain(|_, &mut v| v != id);
            self.params.remove(id);
        }
    }

    /// Iterate all parameters.
    pub fn all_params(&self) -> impl Iterator<Item = (ParameterId, &Parameter)> {
        self.params.iter()
    }

    /// Return parameter IDs scoped to a feature.
    pub fn feature_params(&self, feature: FeatureId) -> Vec<ParameterId> {
        self.params
            .iter()
            .filter(|(_, p)| p.scope == ParamScope::Feature(feature))
            .map(|(id, _)| id)
            .collect()
    }

    /// Return IDs of parameters whose expressions reference this one.
    pub fn dependents(&self, id: ParameterId) -> Vec<ParameterId> {
        self.params
            .iter()
            .filter(|(_, param)| {
                param
                    .expression
                    .as_ref()
                    .map(|expr| expression_param_refs(expr).contains(&id))
                    .unwrap_or(false)
            })
            .map(|(pid, _)| pid)
            .collect()
    }

    /// Set the value of a parameter.
    ///
    /// # Errors
    ///
    /// - [`HistoryError::ReadOnlyParameter`] — parameter is read-only.
    /// - [`HistoryError::DimensionMismatch`] — `unit` has a different
    ///   [`DimensionKind`] than the parameter's existing unit.
    /// - [`HistoryError::FeatureNotFound`] — `id` is invalid.
    pub fn set_value(
        &mut self,
        id: ParameterId,
        value: f64,
        unit: Unit,
    ) -> Result<(), HistoryError> {
        let param = self.params.get(id).ok_or(HistoryError::FeatureNotFound)?;
        if param.read_only {
            return Err(HistoryError::ReadOnlyParameter);
        }
        let existing_kind = param.unit.kind();
        let new_kind = unit.kind();
        if existing_kind != new_kind {
            return Err(HistoryError::DimensionMismatch {
                param: id,
                expected: existing_kind,
                got: new_kind,
            });
        }
        let param = self
            .params
            .get_mut(id)
            .ok_or(HistoryError::FeatureNotFound)?;
        param.value = unit.to_si(value);
        param.unit = unit;
        Ok(())
    }

    /// Set a parameter's value using a value in the given display unit.
    ///
    /// Equivalent to [`set_value`](Self::set_value) — the name makes the
    /// display-unit convention explicit at the call site.
    ///
    /// # Errors
    ///
    /// Same as [`set_value`](Self::set_value).
    pub fn set_value_display(
        &mut self,
        id: ParameterId,
        display_value: f64,
        unit: Unit,
    ) -> Result<(), HistoryError> {
        self.set_value(id, display_value, unit)
    }

    /// Set the expression driving a parameter.
    ///
    /// # Errors
    ///
    /// - [`HistoryError::FeatureNotFound`] if `id` is invalid.
    /// - [`HistoryError::ReadOnlyParameter`] if the parameter is read-only.
    pub fn set_expression(
        &mut self,
        id: ParameterId,
        expr: Expression,
    ) -> Result<(), HistoryError> {
        let param = self
            .params
            .get_mut(id)
            .ok_or(HistoryError::FeatureNotFound)?;
        if param.read_only {
            return Err(HistoryError::ReadOnlyParameter);
        }
        param.expression = Some(expr);
        Ok(())
    }

    /// Clear the expression on a parameter, leaving its current value.
    pub fn clear_expression(&mut self, id: ParameterId) {
        if let Some(param) = self.params.get_mut(id) {
            param.expression = None;
        }
    }

    /// Set optional min/max bounds (in SI units) on a parameter.
    pub fn set_bounds(&mut self, id: ParameterId, min: f64, max: f64) {
        if let Some(param) = self.params.get_mut(id) {
            param.bounds = Some((min, max));
        }
    }

    /// Clear bounds on a parameter.
    pub fn clear_bounds(&mut self, id: ParameterId) {
        if let Some(param) = self.params.get_mut(id) {
            param.bounds = None;
        }
    }

    /// Return the qualified display name, or `None` if `id` is invalid.
    pub fn qualified_name(&self, id: ParameterId) -> Option<String> {
        self.name_index
            .iter()
            .find(|(_, v)| **v == id)
            .map(|(k, _)| k.clone())
    }

    /// Replace the name-index entry for `id` with `qualified_name`.
    pub(crate) fn update_name_index(&mut self, id: ParameterId, qualified_name: &str) {
        // Remove old key(s) pointing to this id.
        self.name_index.retain(|_, &mut v| v != id);
        self.name_index.insert(qualified_name.to_string(), id);
    }

    /// Remove all `name_index` entries pointing to `id`, without touching the
    /// SlotMap.  Used by undo to "hide" a parameter without destroying its key.
    pub(crate) fn deregister_name(&mut self, id: ParameterId) {
        self.name_index.retain(|_, &mut v| v != id);
    }

    /// Re-insert a `name -> id` mapping.  Used by undo/redo to restore a
    /// parameter that was hidden via [`deregister_name`].
    ///
    /// Silently ignores conflicts (the name may already be present if redo
    /// is called repeatedly).
    pub(crate) fn reregister_name(&mut self, id: ParameterId, name: &str) {
        self.name_index.entry(name.to_string()).or_insert(id);
    }

    /// Rename all parameters scoped to a feature when the feature is renamed.
    ///
    /// Updates every key in `name_index` that starts with `"<old_label>."` to
    /// start with `"<new_label>."` instead.
    pub(crate) fn rename_feature_params(
        &mut self,
        feature: FeatureId,
        old_label: &str,
        new_label: &str,
    ) {
        let old_prefix = format!("{old_label}.");
        let new_prefix = format!("{new_label}.");

        let to_rename: Vec<(String, ParameterId)> = self
            .name_index
            .iter()
            .filter(|(k, v)| {
                k.starts_with(&old_prefix)
                    && self
                        .params
                        .get(**v)
                        .map(|p| p.scope == ParamScope::Feature(feature))
                        .unwrap_or(false)
            })
            .map(|(k, v)| (k.clone(), *v))
            .collect();

        for (old_key, param_id) in to_rename {
            let new_key = format!("{new_prefix}{}", &old_key[old_prefix.len()..]);
            self.name_index.remove(&old_key);
            self.name_index.insert(new_key, param_id);
        }
    }
}

#[cfg(test)]
#[path = "parameter_tests.rs"]
mod tests;
