//! [`Sketch`] methods for expression-driven dimensions.
//!
//! Expressions let constraint dimension values be computed from formulas
//! that reference other constraints' target values (e.g. "distance A =
//! 2 * distance B").

use crate::Sketch;
use crate::config::get_constraint_value;
use crate::constraint::ConstraintHandle;
use crate::error::ConstraintError;
use crate::expression::{self, Expression};
use serde::Deserialize;

impl Sketch {
    /// Set a dimension expression on a constraint.
    ///
    /// The expression is evaluated before each solve, writing the result
    /// into the constraint's target value. Returns an error if the
    /// constraint does not exist, is not dimensional, or if the expression
    /// would create a circular dependency.
    pub fn set_dimension_expression(
        &mut self,
        handle: crate::ConstraintHandle,
        expr: Expression,
    ) -> Result<(), ConstraintError> {
        // Validate constraint exists.
        let entry = self
            .constraints
            .get(handle)
            .ok_or(ConstraintError::InvalidHandle)?;

        // Validate constraint is dimensional.
        if get_constraint_value(&entry.constraint).is_none() {
            return Err(ConstraintError::NotDimensional);
        }

        // Check for cycles.
        if expression::has_cycle(handle, &expr, &self.expressions) {
            return Err(ConstraintError::CircularExpression);
        }

        if self.recording {
            self.push_undo_snapshot();
        }
        self.expressions.insert(handle, expr);
        Ok(())
    }

    /// Clear any expression on a constraint, reverting to its current
    /// constant value.
    ///
    /// Does nothing if the handle has no expression.
    pub fn clear_dimension_expression(&mut self, handle: crate::ConstraintHandle) {
        if self.expressions.contains_key(&handle) {
            if self.recording {
                self.push_undo_snapshot();
            }
            self.expressions.remove(&handle);
        }
    }

    /// Evaluate an expression without storing it.
    ///
    /// Useful for previewing what a formula would produce given the
    /// current sketch state.
    pub fn evaluate_expression(&self, expr: &Expression) -> Result<f64, ConstraintError> {
        expression::evaluate(expr, self).map_err(ConstraintError::ExpressionError)
    }
}

/// Serialize the expressions `HashMap` as a `Vec` of pairs (JSON-compatible).
pub(crate) fn serialize_expressions<S>(
    map: &std::collections::HashMap<ConstraintHandle, Expression>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    use serde::ser::SerializeSeq;
    let mut seq = serializer.serialize_seq(Some(map.len()))?;
    for (k, v) in map {
        seq.serialize_element(&(*k, v))?;
    }
    seq.end()
}

/// Deserialize the expressions `HashMap` from a `Vec` of pairs.
pub(crate) fn deserialize_expressions<'de, D>(
    deserializer: D,
) -> Result<std::collections::HashMap<ConstraintHandle, Expression>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let pairs: Vec<(ConstraintHandle, Expression)> = Deserialize::deserialize(deserializer)?;
    Ok(pairs.into_iter().collect())
}
