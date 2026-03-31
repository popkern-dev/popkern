//! Private replay helpers for [`FeatureTree`].
//!
//! Extracted from `lib.rs` to stay within the 500-line structural limit.

use std::collections::HashMap;

#[cfg(test)]
#[path = "replay_integration_tests.rs"]
mod integration_tests;

#[cfg(test)]
#[path = "config_replay_tests.rs"]
mod config_replay_tests;

#[cfg(test)]
#[path = "reference_replay_tests.rs"]
mod reference_replay_tests;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::config::ConfigId;
use crate::config_eval;
use crate::error::{FeatureResult, ReplayError};
use crate::expression::Expression;
use crate::feature::Feature;
use crate::measure_query::MeasureContext;
use crate::measure_resolve;
use crate::node::{BodyId, FeatureId};
use crate::parameter::{EvaluationWarning, ParameterId};
use crate::unit::{DimensionKind, Unit};
use crate::{FeatureTree, ReplayRangeOutput, replay};

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
        let mut saved: Vec<(ParameterId, Option<Expression>)> = Vec::new();

        // Apply literal overrides — pin value and clear expression.
        for (pid, v) in &literals {
            if let Some(param) = self.parameter_table.get_mut(*pid) {
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

    /// Persist replay results: update body_manager and current_brep.
    pub(crate) fn store_body_breps(&mut self, body_breps: &HashMap<BodyId, BRep>) {
        for (body_id, brep) in body_breps {
            self.body_manager.set_brep(*body_id, brep.clone());
        }
        self.current_brep = body_breps
            .get(&BodyId(0))
            .cloned()
            .or_else(|| body_breps.values().next().cloned());
    }

    /// Execute features `from..to` using per-body BRep maps, caching results
    /// and returning the final body BReps together with per-feature outcomes.
    ///
    /// This is the shared core used by both [`replay`](FeatureTree::replay) and
    /// [`replay_from`](FeatureTree::replay_from).
    ///
    /// Each feature's body BRep is temporarily removed from the map before
    /// calling [`replay::execute_feature`] so that Boolean features can freely
    /// access (and consume) the tool body's BRep without borrow conflicts.
    pub(crate) fn replay_range(
        &mut self,
        body_breps: HashMap<BodyId, BRep>,
        from: usize,
        to: usize,
        tolerance: &Tolerance,
    ) -> Result<ReplayRangeOutput, ReplayError> {
        let mut body_breps = body_breps;
        let mut results = HashMap::new();

        for i in from..to {
            let id = self.timeline[i];
            let node = &self.features[id];

            // Config suppression overrides take priority over the node default.
            let config_suppressed =
                config_eval::resolve_suppression(id, self.active_config, &self.configs);
            let base_suppressed = config_suppressed.unwrap_or(node.suppressed);

            if base_suppressed {
                results.insert(id, FeatureResult::Suppressed);
                continue;
            }

            // Evaluate conditional suppression. Clone suppress_if first to
            // release the immutable borrow on `node` before we need `self`
            // mutably below.
            let suppress_if = node.suppress_if.clone();
            let auto_suppressed = suppress_if
                .as_ref()
                .and_then(|expr| {
                    crate::expression::evaluate(
                        expr,
                        &self.parameter_table,
                        &crate::measure_query::MeasureContext::empty(),
                    )
                    .ok()
                })
                .map(|v| v != 0.0)
                .unwrap_or(false);

            if auto_suppressed {
                results.insert(id, FeatureResult::Suppressed);
                continue;
            }

            // Reference features are evaluated at the tree level — they
            // don't modify the BRep, so we skip execute_feature entirely.
            if let Feature::Reference { name: _, query } = &node.feature.clone() {
                let body = node.body;
                let measure_ctx = MeasureContext {
                    bodies: &body_breps,
                    active_body: body,
                    tolerance,
                };
                let eval_result = measure_resolve::resolve_and_evaluate(query, &measure_ctx);

                let label = node.label.clone();
                let qualified_name = format!("{}.value", label.to_lowercase().replace(' ', "_"));
                let dimension_unit = dimension_kind_to_unit(query.dimension_kind());

                match eval_result {
                    Ok(value) => {
                        self.upsert_reference_param(id, &qualified_name, value, dimension_unit);
                        let history = popkern_core::history::HistoryMap::default();
                        results.insert(
                            id,
                            FeatureResult::Success {
                                history: history.clone(),
                            },
                        );
                        self.features[id].last_result = Some(FeatureResult::Success { history });
                    }
                    Err(e) => {
                        // Check for a last_known stale value.
                        let stale_value = self
                            .parameter_table
                            .resolve_name(&qualified_name)
                            .and_then(|pid| {
                                self.parameter_table.get(pid).and_then(|p| p.last_known)
                            });

                        if let Some(stale) = stale_value {
                            // Use stale value — feature succeeds with a warning.
                            self.upsert_reference_param(id, &qualified_name, stale, dimension_unit);
                            let history = popkern_core::history::HistoryMap::default();
                            results.insert(
                                id,
                                FeatureResult::Success {
                                    history: history.clone(),
                                },
                            );
                            self.features[id].last_result =
                                Some(FeatureResult::Success { history });
                        } else {
                            let error_msg = e.to_string();
                            results.insert(
                                id,
                                FeatureResult::Failed {
                                    error: error_msg.clone(),
                                },
                            );
                            self.features[id].last_result =
                                Some(FeatureResult::Failed { error: error_msg });
                        }
                    }
                }
                continue;
            }

            let mut ctx = TrackingContext::new(node.operation_id);
            let feature = node.feature.clone();
            let body = node.body;

            // Take this body's BRep out of the map so we can pass both
            // `&mut brep` and `&mut body_breps` to execute_feature.
            let mut brep = body_breps
                .remove(&body)
                .unwrap_or_else(|| BRep::new(*tolerance));

            // Snapshot BRep before execution for atomic rollback on failure.
            let brep_backup = brep.clone();

            let exec_result = replay::execute_feature(
                &mut brep,
                &feature,
                &mut ctx,
                tolerance,
                &mut body_breps,
                &self.parameter_table,
            );

            match exec_result {
                Ok(side_effects) => {
                    // Apply deferred side effects only on success.
                    for consumed in &side_effects.consumed_bodies {
                        body_breps.remove(consumed);
                    }
                    let history = ctx.finish();
                    results.insert(
                        id,
                        FeatureResult::Success {
                            history: history.clone(),
                        },
                    );
                    self.features[id].last_result = Some(FeatureResult::Success { history });
                    self.cache.store(self.active_config, id, brep.clone(), i);
                }
                Err(e) => {
                    // Atomically restore the pre-execution BRep so the body
                    // remains in a consistent state for subsequent features.
                    // Side effects NOT applied — tool body remains in the map.
                    brep = brep_backup;
                    let error_msg = e.to_string();
                    results.insert(
                        id,
                        FeatureResult::Failed {
                            error: error_msg.clone(),
                        },
                    );
                    self.features[id].last_result =
                        Some(FeatureResult::Failed { error: error_msg });
                }
            }

            // Put the BRep back (restored to backup on failure, or updated on success).
            body_breps.insert(body, brep);
        }

        Ok((body_breps, results))
    }
}

// -- Reference parameter helpers ----------------------------------------------

/// Map a [`DimensionKind`] to its default SI [`Unit`].
fn dimension_kind_to_unit(kind: DimensionKind) -> Unit {
    match kind {
        DimensionKind::Length => Unit::Meters,
        DimensionKind::Area => Unit::SquareMeters,
        DimensionKind::Volume => Unit::CubicMeters,
        DimensionKind::Angle => Unit::Radians,
        DimensionKind::InverseLength => Unit::PerMeter,
        DimensionKind::MomentOfArea => Unit::MetersToFourth,
        DimensionKind::Dimensionless => Unit::Ratio,
    }
}

impl FeatureTree {
    /// Create or update the read-only parameter for a Reference feature.
    ///
    /// If the parameter already exists (by qualified name), its value and
    /// `last_known` are updated in place.  Otherwise a new feature-scoped
    /// parameter is created.
    fn upsert_reference_param(
        &mut self,
        feature_id: FeatureId,
        qualified_name: &str,
        value: f64,
        unit: Unit,
    ) {
        if let Some(pid) = self.parameter_table.resolve_name(qualified_name) {
            // Update existing parameter.
            if let Some(param) = self.parameter_table.get_mut(pid) {
                param.value = value;
                param.last_known = Some(value);
            }
        } else if let Ok(pid) =
            self.parameter_table
                .add_feature_param(qualified_name, feature_id, value, unit)
            && let Some(param) = self.parameter_table.get_mut(pid)
        {
            // Create new feature-scoped, read-only parameter.
            param.read_only = true;
            param.last_known = Some(value);
        }
    }
}

// -- Config-aware parameter accessors -----------------------------------------
//
// These methods live here rather than in lib.rs to keep lib.rs under the
// 500-line structural limit.

impl FeatureTree {
    /// Effective parameter value for the active config.
    ///
    /// Reads directly from the parameter table, which always reflects the most
    /// recent replay evaluation.  For no active config this returns the default
    /// evaluated value.
    pub fn effective_value(&self, param: ParameterId) -> f64 {
        self.parameter_table
            .get(param)
            .map(|p| p.value)
            .unwrap_or(0.0)
    }

    /// Effective suppression state for a feature considering the active config.
    ///
    /// Config suppression overrides take priority over the node's `suppressed`
    /// field.  Returns `false` if the feature does not exist.
    pub fn effective_suppressed(&self, feature: FeatureId) -> bool {
        if let Some(overridden) =
            config_eval::resolve_suppression(feature, self.active_config, &self.configs)
        {
            return overridden;
        }
        self.features
            .get(feature)
            .map(|n| n.suppressed)
            .unwrap_or(false)
    }

    /// Evaluate all parameters for a specific config (expensive — full evaluation).
    ///
    /// Clones the parameter table, applies overrides for `config`, evaluates
    /// all expressions, and returns the resulting value map.  Does not mutate
    /// `self`.
    pub fn evaluate_config(&self, config: Option<ConfigId>) -> HashMap<ParameterId, f64> {
        config_eval::evaluate_config(config, &self.configs, &self.parameter_table)
    }
}
