//! Private replay helpers for [`FeatureTree`].
//!
//! Extracted from `lib.rs` to stay within the 500-line structural limit.

use std::collections::{HashMap, HashSet};

#[cfg(test)]
#[path = "replay_integration_tests.rs"]
mod integration_tests;

#[cfg(test)]
#[path = "config_replay_tests.rs"]
mod config_replay_tests;

#[cfg(test)]
#[path = "reference_replay_tests.rs"]
mod reference_replay_tests;

#[path = "tree_replay_helpers.rs"]
mod helpers;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::auto_merge;
use crate::config_eval;
use crate::error::{FeatureResult, ReplayError};
use crate::expression::Expression;
use crate::feature::{Feature, OperationMode};
use crate::measure_cache::MeasureCache;
use crate::measure_query::MeasureContext;
use crate::node::{BodyId, FeatureId};
use crate::parameter::{EvaluationWarning, ParamScope, ParameterId};
use crate::replay_report::{
    FeatureEvidence, FeatureReport, ReplayOutcome, ReplayReport, SuppressionReason,
};
use crate::{FeatureTree, ReplayRangeOutput, replay};
use helpers::{dimension_kind_to_unit, is_pattern_feature};

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
        let mut report = ReplayReport::default();
        let mut measure_cache = MeasureCache::new();

        // Collect seed features of all pattern features in the replay range.
        // Seed features should NOT be auto-merged — the pattern handles that.
        let seed_feature_ids: HashSet<FeatureId> = self.timeline[from..to]
            .iter()
            .filter_map(|&fid| {
                if let Some(node) = self.features.get(fid) {
                    match &node.feature {
                        Feature::LinearPattern { seed_features, .. }
                        | Feature::CircularPattern { seed_features, .. }
                        | Feature::CurvePattern { seed_features, .. }
                        | Feature::Mirror { seed_features, .. } => {
                            Some(seed_features.iter().copied())
                        }
                        _ => None,
                    }
                } else {
                    None
                }
            })
            .flatten()
            .collect();

        // Track stale reference parameter IDs for StaleInput propagation.
        let mut stale_param_ids: HashSet<ParameterId> = HashSet::new();
        // Track failed (never-valid) reference feature IDs for suppression cascade.
        let mut failed_ref_ids: HashSet<FeatureId> = HashSet::new();

        for i in from..to {
            let id = self.timeline[i];
            let node = &self.features[id];
            let label = node.label.clone();
            let operation = node.feature.operation_label();
            let suppressed_report = |reason: SuppressionReason| {
                FeatureReport::new(id, &label, operation, ReplayOutcome::Suppressed { reason })
            };

            // Config suppression overrides take priority over the node default.
            let config_suppressed =
                config_eval::resolve_suppression(id, self.active_config, &self.configs);
            let base_suppressed = config_suppressed.unwrap_or(node.suppressed);

            if base_suppressed {
                results.insert(id, FeatureResult::Suppressed);
                report
                    .features
                    .push(suppressed_report(if config_suppressed.is_some() {
                        SuppressionReason::Configuration
                    } else {
                        SuppressionReason::Explicit
                    }));
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
                report
                    .features
                    .push(suppressed_report(SuppressionReason::Conditional));
                continue;
            }

            // Suppress features that depend on a failed (never-valid) Reference.
            // Check both entity-level deps and parameter expression refs.
            if !failed_ref_ids.is_empty() {
                let depends_on_failed = node
                    .depends_on
                    .iter()
                    .any(|dep| failed_ref_ids.contains(dep));
                let params_depend_on_failed = self.parameter_table.feature_params(id).iter().any(|&pid| {
                    self.parameter_table.get(pid)
                        .and_then(|p| p.expression.as_ref())
                        .map(|expr| {
                            crate::expression::expression_param_refs(expr).iter().any(|ref_id| {
                                self.parameter_table.get(*ref_id)
                                    .map(|rp| matches!(rp.scope, ParamScope::Feature(fid) if failed_ref_ids.contains(&fid)))
                                    .unwrap_or(false)
                            })
                        })
                        .unwrap_or(false)
                });

                if depends_on_failed || params_depend_on_failed {
                    results.insert(id, FeatureResult::Suppressed);
                    failed_ref_ids.insert(id); // cascade
                    report.features.push(suppressed_report(
                        SuppressionReason::DependsOnFailedReference,
                    ));
                    continue;
                }
            }

            // Propagate StaleInput warnings to parameters that reference stale measures.
            if !stale_param_ids.is_empty() {
                for param_id in self.parameter_table.feature_params(id) {
                    let refs = self
                        .parameter_table
                        .get(param_id)
                        .and_then(|p| p.expression.as_ref())
                        .map(crate::expression::expression_param_refs)
                        .unwrap_or_default();
                    for ref_id in refs {
                        if stale_param_ids.contains(&ref_id) {
                            let source = self
                                .parameter_table
                                .qualified_name(ref_id)
                                .unwrap_or_else(|| "unknown".to_string());
                            if let Some(p) = self.parameter_table.get_mut(param_id) {
                                p.eval_warning = Some(EvaluationWarning::StaleInput { source });
                            }
                        }
                    }
                }
            }

            // Reference features are evaluated at the tree level — they
            // don't modify the BRep, so we skip execute_feature entirely.
            if let Feature::Reference { name: _, query } = &node.feature.clone() {
                let timer = popkern_core::time::Stopwatch::start();
                let body = node.body;
                let ctx = MeasureContext {
                    bodies: &body_breps,
                    active_body: body,
                    tolerance,
                };
                let eval_result = measure_cache.get_or_evaluate(query, &ctx);
                let label = node.label.clone();
                let qname = format!("{}.value", crate::labels::qualified_prefix(&label));
                let dim_unit = dimension_kind_to_unit(query.dimension_kind());
                let fr = match eval_result {
                    Ok(value) => {
                        self.upsert_reference_param(id, &qname, value, dim_unit);
                        FeatureResult::Success {
                            history: popkern_core::history::HistoryMap::default(),
                        }
                    }
                    Err(e) => {
                        let stale = self.parameter_table.resolve_name(&qname).and_then(|pid| {
                            self.parameter_table.get(pid).and_then(|p| p.last_known)
                        });
                        if let Some(sv) = stale {
                            self.upsert_reference_param(id, &qname, sv, dim_unit);
                            if let Some(pid) = self.parameter_table.resolve_name(&qname) {
                                if let Some(p) = self.parameter_table.get_mut(pid) {
                                    p.eval_warning = Some(EvaluationWarning::MeasureStale {
                                        message: e.to_string(),
                                    });
                                }
                                stale_param_ids.insert(pid);
                            }
                            FeatureResult::Success {
                                history: popkern_core::history::HistoryMap::default(),
                            }
                        } else {
                            failed_ref_ids.insert(id);
                            FeatureResult::Failed {
                                error: e.to_string(),
                            }
                        }
                    }
                };
                let mut fr_report = FeatureReport::new(id, &label, operation, outcome_of(&fr));
                fr_report.elapsed = timer.elapsed();
                report.features.push(fr_report);
                results.insert(id, fr.clone());
                self.features[id].last_result = Some(fr);
                continue;
            }

            // Pattern features are executed by the pattern_execute module.
            if is_pattern_feature(&node.feature) {
                let timer = popkern_core::time::Stopwatch::start();
                let (fr, evidence) =
                    match self.execute_pattern_feature(id, &mut body_breps, tolerance) {
                        Ok(pattern_result) => {
                            measure_cache.invalidate();
                            let h = popkern_core::history::HistoryMap::default();
                            (
                                FeatureResult::Success { history: h },
                                FeatureEvidence::pattern(&pattern_result),
                            )
                        }
                        Err(e) => (
                            FeatureResult::Failed {
                                error: e.to_string(),
                            },
                            FeatureEvidence::empty(),
                        ),
                    };
                let mut fr_report = FeatureReport::new(id, &label, operation, outcome_of(&fr))
                    .with_evidence(evidence);
                fr_report.elapsed = timer.elapsed();
                report.features.push(fr_report);
                results.insert(id, fr.clone());
                self.features[id].last_result = Some(fr);
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

            let timer = popkern_core::time::Stopwatch::start();
            let exact_profile = match &feature {
                Feature::Extrude { face, .. } => crate::replay_exact_extrude::exact_profile_for(
                    self.features.values().map(|n| (&n.feature, n.operation_id)),
                    face,
                    &brep,
                    &self.parameter_table,
                ),
                _ => None,
            };
            let exec_result = replay::execute_feature(
                &mut brep,
                &feature,
                &mut ctx,
                tolerance,
                &mut body_breps,
                &self.parameter_table,
                exact_profile.as_ref(),
            );

            match exec_result {
                Ok(side_effects) => {
                    let mut evidence = side_effects.evidence;
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
                    measure_cache.invalidate();

                    // Put the BRep back before auto-merge so detect_host_body can
                    // search all body BReps including this one.
                    body_breps.insert(body, brep);

                    // Auto-merge: Boss/Cut modes merge the tool body into the host.
                    let mode = auto_merge::feature_mode(&feature);
                    if !matches!(mode, OperationMode::NewBody)
                        && !seed_feature_ids.contains(&id)
                        && let Some(host) =
                            auto_merge::detect_host_body(id, &self.features, &body_breps)
                        && host != body
                    {
                        let mut merge_ctx = TrackingContext::new(self.features[id].operation_id);
                        let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
                        auto_merge::execute_auto_merge(
                            mode,
                            body,
                            host,
                            &mut body_breps,
                            &mut merge_ctx,
                            tolerance,
                        )?;
                        evidence.record_boolean_fallbacks(
                            popkern_bool::merge_fallback::take_merge_fallback_events(),
                        );
                    }

                    let mut fr_report =
                        FeatureReport::new(id, &label, operation, ReplayOutcome::Served)
                            .with_evidence(evidence);
                    fr_report.elapsed = timer.elapsed();
                    report.features.push(fr_report);

                    // Skip the normal BRep re-insert below — already handled above.
                    continue;
                }
                Err(e) => {
                    // Atomically restore the pre-execution BRep so the body
                    // remains in a consistent state for subsequent features.
                    // Side effects NOT applied — tool body remains in the map.
                    brep = brep_backup;
                    let error_msg = e.to_string();
                    let mut fr_report = FeatureReport::new(
                        id,
                        &label,
                        operation,
                        ReplayOutcome::Failed {
                            error: error_msg.clone(),
                        },
                    );
                    fr_report.elapsed = timer.elapsed();
                    report.features.push(fr_report);
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

        Ok(ReplayRangeOutput {
            bodies: body_breps,
            results,
            report,
        })
    }
}

/// The report outcome matching a feature result.
fn outcome_of(result: &FeatureResult) -> ReplayOutcome {
    match result {
        FeatureResult::Success { .. } => ReplayOutcome::Served,
        FeatureResult::Failed { error } => ReplayOutcome::Failed {
            error: error.clone(),
        },
        FeatureResult::Suppressed => ReplayOutcome::Suppressed {
            reason: SuppressionReason::Explicit,
        },
    }
}
