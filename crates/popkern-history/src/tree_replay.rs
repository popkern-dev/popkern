//! The replay core of [`FeatureTree`]: execute a range of timeline
//! steps, record each in the step cache, and resume a replay after the
//! longest cached prefix (see `cache` and `step_key`).

use std::collections::{HashMap, HashSet};

use crate::replay_report::ReplayResult;
use crate::{expression, measure_query};

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

#[path = "tree_replay_step.rs"]
mod step;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::FeatureTree;
use crate::ReplayRangeOutput;
use crate::cache::StepRecord;
use crate::document_state::{DocumentState, Touched};
use crate::error::{FeatureResult, ReplayError};
use crate::feature::Feature;
use crate::measure_cache::MeasureCache;
use crate::node::{BodyId, FeatureId};
use crate::parameter::{EvaluationWarning, ParameterId};
use crate::replay_report::ReplayReport;
use crate::step_key::StepKey;
use step::{RangeCarry, StepTouched};

impl FeatureTree {
    /// Persist replay results: update body_manager and current_brep.
    pub(crate) fn store_body_breps(&mut self, body_breps: &HashMap<BodyId, BRep>) {
        // A body absent from this replay's result no longer exists at
        // the end of the timeline (a boolean consumed it, or the
        // rollback bar sits before it); a stale snapshot from an
        // earlier replay must not keep serving it (measured: consumed
        // tool bodies reported alive, PK-385).
        self.body_manager.clear_breps();
        for (body_id, brep) in body_breps {
            self.body_manager.set_brep(*body_id, brep.clone());
        }
        // Fallback: the lowest surviving id — never hash order (PK-525).
        self.current_brep = body_breps.get(&BodyId(0)).cloned().or_else(|| {
            body_breps
                .iter()
                .min_by_key(|(id, _)| id.0)
                .map(|(_, brep)| brep.clone())
        });
    }

    /// Execute timeline steps `from..to` on `start`, recording every
    /// step under its key (`keys` is indexed by timeline position), and
    /// return the final bodies and state with the per-feature outcomes.
    ///
    /// The shared core of every replay entry point. Each feature's body
    /// BRep is temporarily removed from the working map while it
    /// executes so Boolean features can freely access (and consume)
    /// the tool body's BRep without borrow conflicts.
    pub(crate) fn replay_range(
        &mut self,
        start: DocumentState,
        from: usize,
        to: usize,
        keys: &[Option<StepKey>],
        tolerance: &Tolerance,
    ) -> Result<ReplayRangeOutput, ReplayError> {
        let mut body_breps = start.to_working();
        let mut state = start;
        let mut results = HashMap::new();
        let mut report = ReplayReport::default();
        let mut carry = RangeCarry {
            seed_feature_ids: self.pattern_seed_ids(to),
            stale_param_ids: HashSet::new(),
            failed_ref_ids: HashSet::new(),
            measure_cache: MeasureCache::new(),
        };
        let none = HashSet::new();
        for i in from..to {
            let id = self.timeline[i];
            let published = self.replay_step(id, &mut body_breps, &mut carry, tolerance)?;
            let touched = match &published.touched {
                StepTouched::Nothing => Touched::Bodies(&none),
                StepTouched::Bodies(set) => Touched::Bodies(set),
                StepTouched::All => Touched::All,
            };
            state = state.advanced(&body_breps, touched, id);
            if self.frozen_index == Some(i + 1) {
                self.cache.set_frozen(i + 1, state.clone());
            }
            if let Some(Some(key)) = keys.get(i) {
                let mut record =
                    StepRecord::executed(id, published.result.clone(), published.row.clone());
                record.state = Some(state.clone());
                record.last_result = published.last_result;
                record.deferred = self.kept_deferred_body(id);
                record.merged_into = published.merged_into;
                record.reference = published.reference;
                self.cache.insert(*key, i, record);
            }
            results.insert(id, published.result);
            report.features.push(published.row);
        }
        Ok(ReplayRangeOutput {
            bodies: body_breps,
            state,
            results,
            report,
        })
    }

    /// Seed features of the patterns in steps `0..to`: a seed does not
    /// auto-merge (its pattern owns the merge), so whether a step is a
    /// seed depends on LATER steps — the step key hashes this same set.
    pub(crate) fn pattern_seed_ids(&self, to: usize) -> HashSet<FeatureId> {
        self.timeline[..to.min(self.timeline.len())]
            .iter()
            .filter_map(|&fid| match &self.features.get(fid)?.feature {
                Feature::LinearPattern { seed_features, .. }
                | Feature::CircularPattern { seed_features, .. }
                | Feature::CurvePattern { seed_features, .. }
                | Feature::Mirror { seed_features, .. } => Some(seed_features.iter().copied()),
                _ => None,
            })
            .flatten()
            .collect()
    }

    /// Evaluate every parameter the way a replay sees them: under the
    /// active configuration's overrides, or with pinned base values
    /// restored when none is active.
    pub(crate) fn evaluate_parameters_for_replay(
        &mut self,
    ) -> Vec<(ParameterId, EvaluationWarning)> {
        if self.active_config.is_some() {
            self.apply_config_and_evaluate()
        } else {
            self.restore_config_pinned_base();
            expression::evaluate_all(
                &mut self.parameter_table,
                &measure_query::MeasureContext::empty(),
            )
        }
    }

    /// The longest cached prefix a replay to `end` may skip, at most
    /// `limit` steps: the position to execute from and the state there.
    /// Every skipped step must have a resumable record (see
    /// [`StepRecord::resumable`]) and the last one a kept state.
    pub(crate) fn resume_point(
        &mut self,
        keys: &[Option<StepKey>],
        limit: usize,
    ) -> (usize, DocumentState) {
        let mut best = (0, DocumentState::empty());
        for (i, key) in keys.iter().enumerate().take(limit) {
            let Some(record) = key.and_then(|k| self.cache.get(k)) else {
                break;
            };
            if !record.resumable() || record.feature != self.timeline[i] {
                break;
            }
            if let Some(state) = &record.state {
                best = (i + 1, state.clone());
            }
        }
        best
    }

    /// Republish skipped steps `0..from` from their records exactly as
    /// executing them would have: results, node results, deferred-
    /// certification bodies, auto-merge hosts, and their report rows
    /// (returned in timeline order).
    pub(crate) fn republish_prefix(
        &mut self,
        keys: &[Option<StepKey>],
        from: usize,
        results: &mut HashMap<FeatureId, FeatureResult>,
    ) -> ReplayReport {
        let mut rows = ReplayReport::default();
        for key in keys.iter().take(from).flatten() {
            let Some(record) = self.cache.get(*key).cloned() else {
                continue;
            };
            let id = record.feature;
            results.insert(id, record.result);
            if let Some(last) = record.last_result {
                self.features[id].last_result = Some(last);
            }
            self.forget_deferred_body(id);
            if let Some(body) = record.deferred {
                self.deferred.bodies.insert(id, body);
            }
            if let Some(host) = record.merged_into {
                self.merged_into.insert(id, host);
            }
            rows.features.extend(record.row);
        }
        rows
    }

    /// Replay steps `0..end`, resuming after the longest cached prefix
    /// of at most `limit` steps (`limit = 0` executes everything), and
    /// publish the result: bodies, report, and the cache chain. The
    /// returned report carries the republished rows too when
    /// `report_cached` is set (a whole-document replay reports every
    /// feature); otherwise only the executed steps' rows.
    pub(crate) fn replay_resuming(
        &mut self,
        end: usize,
        limit: usize,
        report_cached: bool,
        tolerance: Tolerance,
    ) -> Result<ReplayResult, ReplayError> {
        let start = popkern_core::time::Stopwatch::start();
        let eval_warnings = self.evaluate_parameters_for_replay();
        let keys = self.step_keys(end, &tolerance);
        let (from, state) = match self.frozen_index.map(|f| f.min(end)) {
            Some(f) if f > 0 && self.cache.frozen(f).is_some() => {
                (f, self.cache.frozen(f).cloned().unwrap_or_default())
            }
            _ => self.resume_point(&keys, limit.min(end)),
        };
        let mut out = self.replay_range(state, from, end, &keys, &tolerance)?;
        let mut cached = self.republish_prefix(&keys, from, &mut out.results);
        self.store_body_breps(&out.bodies);
        self.last_report.absorb(cached.clone(), &self.timeline);
        self.last_report.absorb(out.report.clone(), &self.timeline);
        if report_cached {
            cached.features.append(&mut out.report.features);
            out.report = cached;
        }
        let chain = self.timeline[..end].iter().copied().zip(keys).collect();
        self.cache.set_chain(chain);
        Ok(ReplayResult {
            feature_results: out.results,
            elapsed: start.elapsed(),
            eval_warnings,
            report: out.report,
            tolerance,
            cached_steps: from,
        })
    }
}
