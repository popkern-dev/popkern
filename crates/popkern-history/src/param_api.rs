//! Undo-safe parameter mutation API for [`crate::FeatureTree`].
//!
//! All methods in this module record a [`crate::undo::TreeCommand`] entry so
//! that changes can be reversed with [`crate::FeatureTree::undo`].

use crate::FeatureTree;
use crate::error::HistoryError;
use crate::expression::Expression;
use crate::node::FeatureId;
use crate::parameter::ParameterId;
use crate::undo::TreeCommand;
use crate::unit::{DimensionKind, Unit};

impl FeatureTree {
    /// Set an expression on a parameter. Records undo.
    ///
    /// Validates dimensional consistency before accepting: the expression's
    /// [`DimensionKind`] must match the parameter's unit, unless the expression
    /// is [`DimensionKind::Dimensionless`] (literal numbers are allowed on any
    /// parameter).
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::FeatureNotFound`] if `id` is invalid,
    /// [`HistoryError::ReadOnlyParameter`] if the parameter is read-only, or
    /// [`HistoryError::DimensionMismatch`] if the expression produces a
    /// dimension incompatible with the parameter's unit.
    pub fn set_parameter_expression(
        &mut self,
        id: ParameterId,
        expr: Expression,
    ) -> Result<(), HistoryError> {
        let param = self
            .parameter_table
            .get(id)
            .ok_or(HistoryError::FeatureNotFound)?;
        let param_kind = param.unit.kind();
        let old = param.expression.clone();

        // Validate dimensional consistency.
        let expr_kind = crate::expression::dimension_kind(&expr, &self.parameter_table)?;
        if param_kind != DimensionKind::Dimensionless
            && expr_kind != param_kind
            && expr_kind != DimensionKind::Dimensionless
        {
            return Err(HistoryError::DimensionMismatch {
                param: id,
                expected: param_kind,
                got: expr_kind,
            });
        }

        self.parameter_table.set_expression(id, expr.clone())?;
        self.record_undo(TreeCommand::SetExpression {
            param: id,
            old,
            new: Some(expr),
        });
        Ok(())
    }

    /// Clear the expression on a parameter, leaving its current value. Records undo.
    pub fn clear_parameter_expression(&mut self, id: ParameterId) {
        let old = self
            .parameter_table
            .get(id)
            .and_then(|p| p.expression.clone());
        self.parameter_table.clear_expression(id);
        self.record_undo(TreeCommand::SetExpression {
            param: id,
            old,
            new: None,
        });
    }

    /// Set a parameter's value and unit. Records undo.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::FeatureNotFound`] if `id` is invalid,
    /// [`HistoryError::ReadOnlyParameter`] if the parameter is read-only, or
    /// [`HistoryError::DimensionMismatch`] if `unit` has the wrong dimension.
    pub fn set_parameter_value(
        &mut self,
        id: ParameterId,
        value: f64,
        unit: Unit,
    ) -> Result<(), HistoryError> {
        let param = self
            .parameter_table
            .get(id)
            .ok_or(HistoryError::FeatureNotFound)?;
        let old_value = param.value;
        let old_unit = param.unit;
        self.parameter_table.set_value(id, value, unit)?;
        self.record_undo(TreeCommand::SetParamValue {
            param: id,
            old_value,
            old_unit,
            new_value: value,
            new_unit: unit,
        });
        Ok(())
    }

    /// Set a parameter's value using a display-unit value. Records undo.
    ///
    /// Equivalent to [`set_parameter_value`](Self::set_parameter_value) — the
    /// name makes the display-unit convention explicit at the call site.
    ///
    /// # Errors
    ///
    /// Same as [`set_parameter_value`](Self::set_parameter_value).
    pub fn set_parameter_value_display(
        &mut self,
        id: ParameterId,
        display_value: f64,
        unit: Unit,
    ) -> Result<(), HistoryError> {
        self.set_parameter_value(id, display_value, unit)
    }

    /// Set min/max bounds on a parameter. Records undo.
    pub fn set_parameter_bounds(&mut self, id: ParameterId, min: f64, max: f64) {
        let old = self.parameter_table.get(id).and_then(|p| p.bounds);
        self.parameter_table.set_bounds(id, min, max);
        self.record_undo(TreeCommand::SetBounds {
            param: id,
            old,
            new: Some((min, max)),
        });
    }

    /// Clear bounds on a parameter. Records undo.
    pub fn clear_parameter_bounds(&mut self, id: ParameterId) {
        let old = self.parameter_table.get(id).and_then(|p| p.bounds);
        self.parameter_table.clear_bounds(id);
        self.record_undo(TreeCommand::SetBounds {
            param: id,
            old,
            new: None,
        });
    }

    /// Add a user-defined global parameter. Records undo.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if the name is already
    /// in use.
    pub fn add_global_parameter(
        &mut self,
        name: &str,
        value: f64,
        unit: Unit,
    ) -> Result<ParameterId, HistoryError> {
        let id = self.parameter_table.add_global(name, value, unit)?;
        // Snapshot is guaranteed to exist since we just inserted it.
        let snapshot = self
            .parameter_table
            .get(id)
            .cloned()
            .ok_or(HistoryError::FeatureNotFound)?;
        self.record_undo(TreeCommand::AddGlobalParam { id, snapshot });
        Ok(id)
    }

    /// Add a global parameter using a display-unit value. Records undo.
    ///
    /// The value is automatically converted to SI for internal storage.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::DuplicateParameterName`] if the name is already
    /// in use.
    pub fn add_global_parameter_display(
        &mut self,
        name: &str,
        display_value: f64,
        unit: Unit,
    ) -> Result<ParameterId, HistoryError> {
        self.add_global_parameter(name, unit.to_si(display_value), unit)
    }

    /// Remove a global parameter. Records undo. Returns an error if in use.
    ///
    /// On undo the parameter's name is restored in the index; the SlotMap entry
    /// is kept intact (same pattern as feature delete) so that the
    /// [`ParameterId`] remains valid for redo.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::FeatureNotFound`] if `id` is invalid,
    /// [`HistoryError::ReadOnlyParameter`] if read-only, or
    /// [`HistoryError::ParameterInUse`] if other parameters reference it.
    pub fn remove_global_parameter(&mut self, id: ParameterId) -> Result<(), HistoryError> {
        let snapshot = self
            .parameter_table
            .get(id)
            .ok_or(HistoryError::FeatureNotFound)?
            .clone();
        // Soft-delete: remove from name_index only so undo can restore by name.
        // First validate (read-only / in-use checks mirror `remove`).
        if snapshot.read_only {
            return Err(HistoryError::ReadOnlyParameter);
        }
        let dependents = self.parameter_table.dependents(id);
        if !dependents.is_empty() {
            return Err(HistoryError::ParameterInUse(dependents));
        }
        self.parameter_table.deregister_name(id);
        self.record_undo(TreeCommand::RemoveGlobalParam { id, snapshot });
        Ok(())
    }

    /// Remove a global parameter, clearing all expressions that reference it.
    /// Records undo.
    pub fn remove_global_parameter_cascade(&mut self, id: ParameterId) {
        if let Some(snapshot) = self.parameter_table.get(id).cloned() {
            // Clear dependent expressions before deregistering the name.
            let dependents = self.parameter_table.dependents(id);
            for dep_id in dependents {
                self.parameter_table.clear_expression(dep_id);
            }
            self.parameter_table.deregister_name(id);
            self.record_undo(TreeCommand::RemoveGlobalParam { id, snapshot });
        }
    }

    /// Rename a global parameter. Records undo.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::FeatureNotFound`] if `id` is invalid, or
    /// [`HistoryError::DuplicateParameterName`] if `new_name` is already in use.
    pub fn rename_parameter(
        &mut self,
        id: ParameterId,
        new_name: &str,
    ) -> Result<(), HistoryError> {
        let old_name = self
            .parameter_table
            .get(id)
            .ok_or(HistoryError::FeatureNotFound)?
            .name
            .clone();
        self.parameter_table.rename(id, new_name)?;
        self.record_undo(TreeCommand::RenameParam {
            id,
            old_name,
            new_name: new_name.to_string(),
        });
        Ok(())
    }

    /// Set or clear the suppress_if expression on a feature. Records undo.
    ///
    /// If an expression is provided, it must produce
    /// [`DimensionKind::Dimensionless`] — suppress_if is a boolean condition.
    ///
    /// # Errors
    ///
    /// Returns [`HistoryError::SuppressIfNotDimensionless`] if the expression
    /// has a non-dimensionless dimension.
    pub fn set_suppress_if(
        &mut self,
        feature: FeatureId,
        expr: Option<Expression>,
    ) -> Result<(), HistoryError> {
        if let Some(ref e) = expr {
            let kind = crate::expression::dimension_kind(e, &self.parameter_table)?;
            if kind != DimensionKind::Dimensionless {
                return Err(HistoryError::SuppressIfNotDimensionless(kind));
            }
        }
        let old = self
            .features
            .get(feature)
            .and_then(|n| n.suppress_if.clone());
        if let Some(node) = self.features.get_mut(feature) {
            node.suppress_if = expr.clone();
        }
        self.record_undo(TreeCommand::SetSuppressIf {
            feature,
            old,
            new: expr,
        });
        Ok(())
    }
}

#[cfg(test)]
#[path = "param_api_tests.rs"]
mod tests;
