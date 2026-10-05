//! One timeline step of a replay: execute (or skip) a single feature
//! against the working body map and say what the step published.
//!
//! Split from `tree_replay.rs` (500-line ratchet) when replay began
//! recording every step for the step cache (PK-659): the loop in
//! `replay_range` turns each [`StepPublished`] into a report row, a
//! result, a [`crate::document_state::DocumentState`] and a cache record.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::helpers::{dimension_kind_to_unit, is_pattern_feature};
use crate::auto_merge;
use crate::config_eval;
use crate::error::{FeatureResult, ReplayError};
use crate::feature::{Feature, OperationMode};
use crate::measure_cache::MeasureCache;
use crate::measure_query::MeasureContext;
use crate::node::{BodyId, FeatureId};
use crate::parameter::{EvaluationWarning, ParamScope, ParameterId};
use crate::replay::{self};
use crate::replay_report::{FeatureEvidence, FeatureReport, ReplayOutcome, SuppressionReason};
use crate::{FeatureTree, failure::FeatureFailure};

/// State a replay range carries from step to step besides the bodies.
pub(super) struct RangeCarry {
    /// Seed features of the range's patterns (never auto-merged: the
    /// pattern owns their merge).
    pub(super) seed_feature_ids: HashSet<FeatureId>,
    /// Reference parameters whose value went stale (StaleInput).
    pub(super) stale_param_ids: HashSet<ParameterId>,
    /// Never-valid Reference features (suppression cascade).
    pub(super) failed_ref_ids: HashSet<FeatureId>,
    /// Measurement memo, invalidated whenever a body changes.
    pub(super) measure_cache: MeasureCache,
}

/// Which bodies a step may have changed.
pub(super) enum StepTouched {
    /// None (suppressed, reference, failed).
    Nothing,
    /// Exactly these.
    Bodies(HashSet<BodyId>),
    /// Untracked: every body (patterns).
    All,
}

/// What one step published.
pub(super) struct StepPublished {
    pub(super) result: FeatureResult,
    pub(super) last_result: Option<FeatureResult>,
    pub(super) row: FeatureReport,
    pub(super) touched: StepTouched,
    pub(super) merged_into: Option<BodyId>,
    pub(super) reference: bool,
}

impl StepPublished {
    fn skipped(row: FeatureReport) -> Self {
        StepPublished {
            result: FeatureResult::Suppressed,
            last_result: None,
            row,
            touched: StepTouched::Nothing,
            merged_into: None,
            reference: false,
        }
    }
}

impl FeatureTree {
    /// Execute timeline step `id` on `body_breps` (see the module docs).
    pub(super) fn replay_step(
        &mut self,
        id: FeatureId,
        body_breps: &mut HashMap<BodyId, BRep>,
        carry: &mut RangeCarry,
        tolerance: &Tolerance,
    ) -> Result<StepPublished, ReplayError> {
        // This execution's evidence replaces the feature's report
        // row, so its kept deferred-certification BRep goes too.
        self.forget_deferred_body(id);
        let node = &self.features[id];
        let label = node.label.clone();
        let operation = node.feature.operation_label();
        let suppressed_report = |reason: SuppressionReason| {
            FeatureReport::new(id, &label, operation, ReplayOutcome::Suppressed { reason })
        };

        // Config suppression overrides take priority over the node default.
        let config_suppressed =
            config_eval::resolve_suppression(id, self.active_config, &self.configs);
        if config_suppressed.unwrap_or(node.suppressed) {
            return Ok(StepPublished::skipped(suppressed_report(
                if config_suppressed.is_some() {
                    SuppressionReason::Configuration
                } else {
                    SuppressionReason::Explicit
                },
            )));
        }

        // Conditional suppression.
        let auto_suppressed = node
            .suppress_if
            .as_ref()
            .and_then(|expr| {
                crate::expression::evaluate(expr, &self.parameter_table, &MeasureContext::empty())
                    .ok()
            })
            .is_some_and(|v| v != 0.0);
        if auto_suppressed {
            return Ok(StepPublished::skipped(suppressed_report(
                SuppressionReason::Conditional,
            )));
        }

        // Suppress features that depend on a failed (never-valid) Reference,
        // through entity-level deps or parameter expression refs.
        if !carry.failed_ref_ids.is_empty() {
            let depends_on_failed = node
                .depends_on
                .iter()
                .any(|dep| carry.failed_ref_ids.contains(dep));
            let params_depend_on_failed =
                self.parameter_table.feature_params(id).iter().any(|&pid| {
                    self.parameter_table
                        .get(pid)
                        .and_then(|p| p.expression.as_ref())
                        .is_some_and(|expr| {
                            crate::expression::expression_param_refs(expr)
                                .iter()
                                .any(|ref_id| {
                                    self.parameter_table.get(*ref_id).is_some_and(|rp| {
                                        matches!(rp.scope, ParamScope::Feature(fid) if carry.failed_ref_ids.contains(&fid))
                                    })
                                })
                        })
                });
            if depends_on_failed || params_depend_on_failed {
                carry.failed_ref_ids.insert(id); // cascade
                return Ok(StepPublished::skipped(suppressed_report(
                    SuppressionReason::DependsOnFailedReference,
                )));
            }
        }

        // Propagate StaleInput warnings to parameters that reference stale measures.
        if !carry.stale_param_ids.is_empty() {
            for param_id in self.parameter_table.feature_params(id) {
                let refs = self
                    .parameter_table
                    .get(param_id)
                    .and_then(|p| p.expression.as_ref())
                    .map(crate::expression::expression_param_refs)
                    .unwrap_or_default();
                for ref_id in refs {
                    if carry.stale_param_ids.contains(&ref_id) {
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

        let node = &self.features[id];
        if let Feature::Reference { query, .. } = &node.feature.clone() {
            return Ok(self.reference_step(id, query, body_breps, carry, tolerance));
        }

        // Pattern features are executed by the pattern_execute module.
        if is_pattern_feature(&node.feature) {
            let timer = popkern_core::time::Stopwatch::start();
            let (fr, evidence) = match self.execute_pattern_feature(id, body_breps, tolerance) {
                Ok(pattern_result) => {
                    carry.measure_cache.invalidate();
                    (
                        FeatureResult::Success {
                            history: popkern_core::history::HistoryMap::default(),
                        },
                        FeatureEvidence::pattern(&pattern_result),
                    )
                }
                Err(e) => (
                    FeatureResult::Failed {
                        error: FeatureFailure::from(e),
                    },
                    FeatureEvidence::empty(),
                ),
            };
            let mut row =
                FeatureReport::new(id, &label, operation, outcome_of(&fr)).with_evidence(evidence);
            row.elapsed = timer.elapsed();
            self.features[id].last_result = Some(fr.clone());
            return Ok(StepPublished {
                result: fr.clone(),
                last_result: Some(fr),
                row,
                touched: StepTouched::All,
                merged_into: None,
                reference: false,
            });
        }

        let mut ctx =
            TrackingContext::new(node.operation_id).with_naming_scheme(self.naming_scheme);
        let body = node.body;

        // Take this body's BRep out of the map so we can pass both
        // `&mut brep` and `&mut body_breps` to execute_feature.
        let mut brep = body_breps
            .remove(&body)
            .unwrap_or_else(|| BRep::new(*tolerance));
        // Snapshot BRep before execution for atomic rollback on failure.
        let brep_backup = brep.clone();
        // References an upstream edit renumbered rebind by their
        // recorded provenance; an ambiguous one declines (PK-677).
        let node = &self.features[id];
        let rebound = crate::reference_rebind::rebind(&node.feature, &node.provenance, &brep);
        let (feature, rebind_fallbacks) = match &rebound {
            Ok((feature, fallbacks)) => (feature.clone(), fallbacks.clone()),
            Err(_) => (node.feature.clone(), Vec::new()),
        };

        let timer = popkern_core::time::Stopwatch::start();
        // The body carrying the face this feature's sketch was placed
        // on (detected through the on-face reference): the Boss/Cut
        // merge host, and the BRep an on-face sketch plane or an
        // up-to / offset target resolves against — a tool body's own
        // BRep never carries them. (Pattern seeds detect their host
        // too — their planes and targets live there — but the pattern
        // owns their merge.)
        let merge_host = auto_merge::detect_host_body(id, &self.features, body_breps)
            .or_else(|| auto_merge::host_for_extrude_target(&feature, body, body_breps));
        let exact_profiles = crate::replay_exact_extrude::exact_profile_set(
            &feature,
            self.features.values().map(|n| (&n.feature, n.operation_id)),
            merge_host.and_then(|h| body_breps.get(&h)).unwrap_or(&brep),
            &self.parameter_table,
        );
        let exec_result = match rebound {
            Err(error) => Err(error),
            Ok(_) => replay::execute_feature(
                &mut brep,
                &feature,
                &mut ctx,
                tolerance,
                body_breps,
                &self.parameter_table,
                &exact_profiles,
                merge_host,
                self.deferred.policy,
            ),
        };

        match exec_result {
            Ok(side_effects) => {
                let mut evidence = side_effects.evidence;
                evidence.fallbacks.extend(rebind_fallbacks);
                let recorded = self.features[id].feature.clone();
                crate::reference_rebind::refresh(
                    &mut self.features[id].provenance,
                    &recorded,
                    &brep_backup,
                );
                let mut touched: HashSet<BodyId> = HashSet::from([body]);
                // Apply deferred side effects only on success.
                for consumed in &side_effects.consumed_bodies {
                    body_breps.remove(consumed);
                }
                let history = ctx.finish();
                let fr = FeatureResult::Success { history };
                self.features[id].last_result = Some(fr.clone());
                self.keep_deferred_body(id, &evidence, &brep);
                carry.measure_cache.invalidate();

                // Put the BRep back before auto-merge so detect_host_body can
                // search all body BReps including this one.
                body_breps.insert(body, brep);

                // Auto-merge: Boss/Cut modes merge the tool body into the host.
                let mut merged_into = None;
                let mode = auto_merge::feature_mode(&feature);
                if !matches!(mode, OperationMode::NewBody)
                    && !carry.seed_feature_ids.contains(&id)
                    && let Some(host) = merge_host
                    && host != body
                {
                    let mut merge_ctx = TrackingContext::new(self.features[id].operation_id)
                        .with_naming_scheme(self.naming_scheme);
                    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
                    auto_merge::execute_auto_merge(
                        mode,
                        body,
                        host,
                        body_breps,
                        &mut merge_ctx,
                        tolerance,
                    )?;
                    self.merged_into.insert(id, host);
                    merged_into = Some(host);
                    touched.insert(host);
                    evidence.record_boolean_fallbacks(
                        popkern_bool::merge_fallback::take_merge_fallback_events(),
                    );
                }

                let mut row = FeatureReport::new(id, &label, operation, ReplayOutcome::Served)
                    .with_evidence(evidence);
                row.elapsed = timer.elapsed();
                Ok(StepPublished {
                    result: fr.clone(),
                    last_result: Some(fr),
                    row,
                    touched: StepTouched::Bodies(touched),
                    merged_into,
                    reference: false,
                })
            }
            Err(e) => {
                // Atomically restore the pre-execution BRep so the body
                // remains in a consistent state for subsequent features.
                // Side effects NOT applied — tool body remains in the map.
                let error = FeatureFailure::from(e);
                let mut row = FeatureReport::new(
                    id,
                    &label,
                    operation,
                    ReplayOutcome::Failed {
                        error: error.clone(),
                    },
                );
                row.elapsed = timer.elapsed();
                let fr = FeatureResult::Failed { error };
                self.features[id].last_result = Some(fr.clone());
                body_breps.insert(body, brep_backup);
                Ok(StepPublished {
                    result: fr.clone(),
                    last_result: Some(fr),
                    row,
                    touched: StepTouched::Nothing,
                    merged_into: None,
                    reference: false,
                })
            }
        }
    }

    /// Evaluate a Reference feature at the tree level: it changes no
    /// body, only its read-only parameter.
    fn reference_step(
        &mut self,
        id: FeatureId,
        query: &crate::measure_query::MeasureQuery,
        body_breps: &HashMap<BodyId, BRep>,
        carry: &mut RangeCarry,
        tolerance: &Tolerance,
    ) -> StepPublished {
        let timer = popkern_core::time::Stopwatch::start();
        let node = &self.features[id];
        let operation = node.feature.operation_label();
        let ctx = MeasureContext {
            bodies: body_breps,
            active_body: node.body,
            tolerance,
        };
        let eval_result = carry.measure_cache.get_or_evaluate(query, &ctx);
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
                let stale = self
                    .parameter_table
                    .resolve_name(&qname)
                    .and_then(|pid| self.parameter_table.get(pid).and_then(|p| p.last_known));
                if let Some(sv) = stale {
                    self.upsert_reference_param(id, &qname, sv, dim_unit);
                    if let Some(pid) = self.parameter_table.resolve_name(&qname) {
                        if let Some(p) = self.parameter_table.get_mut(pid) {
                            p.eval_warning = Some(EvaluationWarning::MeasureStale {
                                message: e.to_string(),
                            });
                        }
                        carry.stale_param_ids.insert(pid);
                    }
                    FeatureResult::Success {
                        history: popkern_core::history::HistoryMap::default(),
                    }
                } else {
                    carry.failed_ref_ids.insert(id);
                    FeatureResult::Failed {
                        error: FeatureFailure::from(e),
                    }
                }
            }
        };
        let mut row = FeatureReport::new(id, &label, operation, outcome_of(&fr));
        row.elapsed = timer.elapsed();
        self.features[id].last_result = Some(fr.clone());
        StepPublished {
            result: fr.clone(),
            last_result: Some(fr),
            row,
            touched: StepTouched::Nothing,
            merged_into: None,
            reference: true,
        }
    }

    /// The deferred-certification BRep the last step kept for `id`.
    pub(super) fn kept_deferred_body(&self, id: FeatureId) -> Option<Arc<BRep>> {
        self.deferred.bodies.get(&id).cloned()
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
