//! Configuration-aware parameter evaluation for replay (split from
//! `tree_replay.rs` under the file-size ratchet).

use crate::FeatureTree;
use crate::config_eval;
use crate::expression::Expression;
use crate::parameter::{EvaluationWarning, ParameterId};

impl FeatureTree {
    /// Apply config overrides to the parameter table and evaluate all expressions.
    ///
    /// Saves the original default expressions, swaps in config overrides for
    /// evaluation, runs [`expression::evaluate_all`], then restores the default
    /// expressions so the table's expression fields always represent the
    /// default model state.  The resolved `value` fields are left in place so
    /// [`replay::execute_feature`] reads the config-overridden values.
    pub(crate) fn apply_config_and_evaluate(&mut self) -> Vec<(ParameterId, EvaluationWarning)> {
        let (literals, expr_overrides) = config_eval::build_effective_expressions(
            self.active_config,
            &self.configs,
            &self.parameter_table,
        );

        // Track what we change so we can restore expressions afterwards.
        // A literal override pins a value that has no expression to
        // re-derive it from, so its BASE value is remembered in
        // `config_pinned_base` (persisted) and restored by the next base
        // replay — leaving the effective value in the table meanwhile,
        // which is what a GUI shows under an active configuration.
        // (Measured before: `size` stayed 2.0 after returning from the
        // `tall` configuration, PK-384.)
        let mut saved: Vec<(ParameterId, Option<Expression>)> = Vec::new();

        // Apply literal overrides — pin value and clear expression.
        for (pid, v) in &literals {
            if let Some(param) = self.parameter_table.get_mut(*pid) {
                if param.expression.is_none() {
                    self.config_pinned_base.entry(*pid).or_insert(param.value);
                }
                saved.push((*pid, param.expression.clone()));
                param.value = *v;
                param.expression = None;
            }
        }

        // Swap in config expression overrides.
        for (pid, expr) in &expr_overrides {
            if let Some(param) = self.parameter_table.get_mut(*pid) {
                if !saved.iter().any(|(id, _)| id == pid) {
                    saved.push((*pid, param.expression.clone()));
                }
                param.expression = Some(expr.clone());
            }
        }

        // Evaluate — literal-overridden params have no expression, so they
        // keep their pinned value; expression-overridden params are evaluated
        // in topological order.
        let warnings = crate::expression::evaluate_all(
            &mut self.parameter_table,
            &crate::measure_query::MeasureContext::empty(),
        );

        // Restore default expressions (leave evaluated values in place).
        for (pid, original_expr) in saved {
            if let Some(param) = self.parameter_table.get_mut(pid) {
                param.expression = original_expr;
            }
        }

        warnings
    }

    /// Put back the base values a configuration's literal overrides
    /// pinned (see `config_pinned_base`); called by every base replay.
    pub(crate) fn restore_config_pinned_base(&mut self) {
        for (pid, base) in std::mem::take(&mut self.config_pinned_base) {
            if let Some(param) = self.parameter_table.get_mut(pid)
                && param.expression.is_none()
            {
                param.value = base;
            }
        }
    }
}
