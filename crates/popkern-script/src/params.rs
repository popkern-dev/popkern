//! Parameter verbs beyond declaration and value edits (spec §4.1):
//! expressions, bounds, and renames on any parameter by qualified
//! name.

use popkern_history::unit::Unit;

use crate::error::ScriptError;
use crate::model::Model;

/// A configuration as the surface lists it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConfigInfo {
    /// Configuration name.
    pub name: String,
    /// Description text.
    pub description: String,
    /// Parent configuration name, for derived configurations.
    pub parent: Option<String>,
    /// Parameters overridden (qualified names).
    pub overrides: Vec<String>,
    /// Features suppressed (labels) and their suppression state.
    pub suppressions: Vec<(String, bool)>,
}

impl Model {
    fn param_by_name(
        &self,
        name: &str,
    ) -> Result<popkern_history::parameter::ParameterId, ScriptError> {
        self.tree
            .parameters()
            .resolve_name(name)
            .ok_or_else(|| ScriptError::Argument {
                argument: "parameter".into(),
                reason: format!("unknown parameter '{name}'"),
            })
    }

    /// Drive a parameter by an expression.
    pub fn set_expr(&mut self, name: &str, text: &str) -> Result<(), ScriptError> {
        let id = self.param_by_name(name)?;
        let expression = popkern_history::parse::parse_expression(text, self.tree.parameters())
            .map_err(|e| ScriptError::Argument {
                argument: name.to_string(),
                reason: format!("expression '{text}': {e}"),
            })?;
        self.tree.set_parameter_expression(id, expression)?;
        Ok(())
    }

    /// Remove a parameter's expression (it keeps its last value).
    pub fn clear_expr(&mut self, name: &str) -> Result<(), ScriptError> {
        let id = self.param_by_name(name)?;
        self.tree.clear_parameter_expression(id);
        Ok(())
    }

    /// Declare bounds in the parameter's display unit.
    pub fn set_bounds(&mut self, name: &str, min: f64, max: f64) -> Result<(), ScriptError> {
        let id = self.param_by_name(name)?;
        if !(min.is_finite() && max.is_finite()) || min > max {
            return Err(ScriptError::Argument {
                argument: "bounds".into(),
                reason: format!("bounds must be finite with min <= max, got [{min}, {max}]"),
            });
        }
        let unit = self
            .tree
            .parameters()
            .get(id)
            .map(|p| p.unit)
            .unwrap_or(Unit::Meters);
        self.tree
            .set_parameter_bounds(id, unit.to_si(min), unit.to_si(max));
        Ok(())
    }

    /// Remove a parameter's bounds.
    pub fn clear_bounds(&mut self, name: &str) -> Result<(), ScriptError> {
        let id = self.param_by_name(name)?;
        self.tree.clear_parameter_bounds(id);
        Ok(())
    }

    /// Rename a global parameter (references follow the id).
    pub fn rename_param(&mut self, old: &str, new: &str) -> Result<(), ScriptError> {
        let id = self.param_by_name(old)?;
        self.tree.rename_parameter(id, new)?;
        Ok(())
    }
}
