//! The step cache: what every timeline step of a replay published,
//! content-addressed by [`StepKey`].
//!
//! Each [`StepRecord`] holds a step's [`DocumentState`] (every alive
//! body, shared with neighbouring steps) and everything the step
//! published besides geometry — its feature result, report row,
//! deferred-certification body and auto-merge host — so a replay that
//! resumes after the step republishes it exactly as a fresh replay
//! would. Keys chain every input a step reads (see `step_key`), so a
//! record can never serve stale geometry: an edit changes the keys
//! from the first affected step on, and those steps miss.
//!
//! The cache also remembers the CHAIN of the last published replay
//! (one key per timeline position), which is what the history scrubber
//! and the in-place parameter route read.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use crate::document_state::DocumentState;
use crate::error::FeatureResult;
use crate::node::{BodyId, FeatureId};
use crate::replay_report::FeatureReport;
use crate::step_key::StepKey;

// ---------------------------------------------------------------------------
// CacheStrategy
// ---------------------------------------------------------------------------

/// Controls which intermediate document states are retained in memory.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum CacheStrategy {
    /// Keep the document state after every feature — instant scrubbing at
    /// the cost of higher memory usage.
    EveryStep,
    /// Keep the state after every *N*th feature — balanced trade-off between
    /// memory and replay latency.
    EveryNth(usize),
    /// Never keep states automatically; one can still be pinned at a
    /// step of the last replay (the in-place parameter route does).
    OnDemand,
    /// Keep states until their estimated footprint exceeds `max_bytes`,
    /// then drop the least recently used (shared bodies count once).
    MemoryBudget {
        /// Maximum estimated size of the kept states in bytes.
        max_bytes: usize,
    },
}

/// What one executed step published (see the module docs).
#[derive(Debug, Clone)]
pub(crate) struct StepRecord {
    /// The timeline feature the step executed.
    pub(crate) feature: FeatureId,
    /// The document after the step, when the strategy keeps it.
    pub(crate) state: Option<DocumentState>,
    /// The step's entry in the replay's per-feature results.
    pub(crate) result: FeatureResult,
    /// What the step wrote to the node's `last_result` (`None`: the
    /// step left it untouched, as suppression does).
    pub(crate) last_result: Option<FeatureResult>,
    /// The step's report row; `None` for a state pinned without an
    /// execution (the in-place route), which a replay cannot resume
    /// past because it has no evidence to republish.
    pub(crate) row: Option<FeatureReport>,
    /// The BRep a deferred-certification feature keeps.
    pub(crate) deferred: Option<Arc<BRep>>,
    /// The body a Boss/Cut feature merged into.
    pub(crate) merged_into: Option<BodyId>,
    /// Whether the step evaluated a Reference feature. Reference
    /// values and the stale-input warnings they propagate live in the
    /// parameter table, not here, so a replay never resumes past one.
    pub(crate) reference: bool,
    used: u64,
}

impl StepRecord {
    /// A record of an executed step.
    pub(crate) fn executed(feature: FeatureId, result: FeatureResult, row: FeatureReport) -> Self {
        StepRecord {
            feature,
            state: None,
            result,
            last_result: None,
            row: Some(row),
            deferred: None,
            merged_into: None,
            reference: false,
            used: 0,
        }
    }

    /// A state pinned without an execution (no row to republish).
    pub(crate) fn state_only(feature: FeatureId, state: DocumentState) -> Self {
        StepRecord {
            feature,
            state: Some(state),
            result: FeatureResult::Suppressed,
            last_result: None,
            row: None,
            deferred: None,
            merged_into: None,
            reference: false,
            used: 0,
        }
    }

    /// Whether a replay may resume after this step: it executed (has a
    /// row to republish) and is not a Reference evaluation.
    pub(crate) fn resumable(&self) -> bool {
        self.row.is_some() && !self.reference
    }
}

/// The step cache (see the module docs).
#[derive(Debug, Clone)]
pub(crate) struct StepCache {
    records: HashMap<StepKey, StepRecord>,
    strategy: CacheStrategy,
    tick: u64,
    /// The last published replay's keys, one per timeline position.
    chain: Vec<(FeatureId, Option<StepKey>)>,
    /// The state at the freeze point, pinned when a replay passes it.
    frozen: Option<(usize, DocumentState)>,
    /// Per feature, the step whose record last wrote the node's
    /// `last_result` (executed under that key or republished from it).
    /// A republish from the same key would rewrite an equal value, so
    /// it skips the copy (PK-805: every recorded feature re-copied the
    /// naming history of every earlier one). Any other write clears
    /// the entry; an empty map only costs copies.
    node_results: HashMap<FeatureId, StepKey>,
}

impl Default for StepCache {
    fn default() -> Self {
        Self::new(CacheStrategy::EveryStep)
    }
}

impl StepCache {
    /// A new, empty cache with the given strategy.
    pub(crate) fn new(strategy: CacheStrategy) -> Self {
        Self {
            records: HashMap::new(),
            strategy,
            tick: 0,
            chain: Vec::new(),
            frozen: None,
            node_results: HashMap::new(),
        }
    }

    /// The record stored under `key`, marked as used.
    pub(crate) fn get(&mut self, key: StepKey) -> Option<&StepRecord> {
        self.tick += 1;
        let tick = self.tick;
        self.records.get_mut(&key).map(|r| {
            r.used = tick;
            &*r
        })
    }

    /// Store `record` for the step at timeline `index`; the strategy
    /// decides whether its state is kept.
    pub(crate) fn insert(&mut self, key: StepKey, index: usize, mut record: StepRecord) {
        let keep = match self.strategy {
            CacheStrategy::EveryStep | CacheStrategy::MemoryBudget { .. } => true,
            CacheStrategy::EveryNth(n) => n > 0 && index.is_multiple_of(n),
            CacheStrategy::OnDemand => false,
        };
        if !keep {
            record.state = None;
        }
        self.tick += 1;
        record.used = self.tick;
        self.records.insert(key, record);
    }

    /// Publish the chain of the replay that just completed, then trim
    /// records the strategy no longer keeps.
    pub(crate) fn set_chain(&mut self, chain: Vec<(FeatureId, Option<StepKey>)>) {
        self.chain = chain;
        self.trim();
    }

    /// The last published chain.
    pub(crate) fn chain(&self) -> &[(FeatureId, Option<StepKey>)] {
        &self.chain
    }

    /// The record of `feature` in the last published chain.
    pub(crate) fn chain_record_mut(&mut self, feature: FeatureId) -> Option<&mut StepRecord> {
        let key = self
            .chain
            .iter()
            .find(|(id, _)| *id == feature)
            .and_then(|(_, key)| *key)?;
        self.records.get_mut(&key)
    }

    /// The state after timeline position `index` of the last published
    /// chain, when it is kept.
    pub(crate) fn chain_state(&self, index: usize) -> Option<&DocumentState> {
        let key = self.chain.get(index)?.1?;
        self.records.get(&key)?.state.as_ref()
    }

    /// Pin `state` at position `index` of the last published chain
    /// (the in-place parameter route's refreshed state, which has no
    /// report row and so ends any resume there).
    pub(crate) fn pin_chain_state(&mut self, index: usize, key: StepKey, state: DocumentState) {
        let Some(&(feature, _)) = self.chain.get(index) else {
            return;
        };
        self.chain[index].1 = Some(key);
        self.chain.truncate(index + 1);
        let mut record = StepRecord::state_only(feature, state);
        self.tick += 1;
        record.used = self.tick;
        self.records.insert(key, record);
    }

    /// The state pinned at the freeze point, when its position is `index`.
    pub(crate) fn frozen(&self, index: usize) -> Option<&DocumentState> {
        self.frozen
            .as_ref()
            .filter(|(at, _)| *at == index)
            .map(|(_, state)| state)
    }

    /// Pin the freeze point's state.
    pub(crate) fn set_frozen(&mut self, index: usize, state: DocumentState) {
        self.frozen = Some((index, state));
    }

    /// Forget the pinned freeze-point state.
    pub(crate) fn clear_frozen(&mut self) {
        self.frozen = None;
    }

    /// Remove every record, the chain and the frozen state.
    /// Record which step's record the node of `id` now holds as its
    /// `last_result` (`None`: a write no record matches).
    pub(crate) fn note_node_result(&mut self, id: FeatureId, key: Option<StepKey>) {
        match key {
            Some(key) => {
                self.node_results.insert(id, key);
            }
            None => {
                self.node_results.remove(&id);
            }
        }
    }

    /// Whether the node of `id` holds the `last_result` the record of
    /// `key` wrote.
    pub(crate) fn node_holds(&self, id: FeatureId, key: StepKey) -> bool {
        self.node_results.get(&id) == Some(&key)
    }

    pub(crate) fn clear(&mut self) {
        self.records.clear();
        self.chain.clear();
        self.frozen = None;
        self.node_results.clear();
    }

    /// Estimated bytes of every kept state, shared bodies counted once.
    pub(crate) fn memory_bytes(&self) -> usize {
        let mut seen = HashSet::new();
        self.records
            .values()
            .filter_map(|r| r.state.as_ref())
            .map(|s| s.unique_bytes(&mut seen))
            .sum()
    }

    /// The record under `key` without marking it used.
    #[cfg(test)]
    pub(crate) fn peek(&self, key: StepKey) -> Option<&StepRecord> {
        self.records.get(&key)
    }

    /// Replace the chain with hand-authored state-only steps.
    #[cfg(test)]
    pub(crate) fn plant_chain(&mut self, steps: Vec<(FeatureId, StepKey, DocumentState)>) {
        self.chain = steps.iter().map(|(id, key, _)| (*id, Some(*key))).collect();
        for (id, key, state) in steps {
            self.records.insert(key, StepRecord::state_only(id, state));
        }
    }

    /// Replace the strategy and clear the cache.
    pub(crate) fn set_strategy(&mut self, strategy: CacheStrategy) {
        self.strategy = strategy;
        self.clear();
    }

    /// Keep the published chain's records plus, off the chain, at most
    /// as many recently used records as the chain is long (an earlier
    /// variant a user toggles back to stays warm). Under a memory
    /// budget, drop least-recently-used states until the kept states
    /// fit; their records stay (rows and results are small).
    fn trim(&mut self) {
        let on_chain: HashSet<StepKey> = self.chain.iter().filter_map(|(_, k)| *k).collect();
        let mut off_chain: Vec<(u64, StepKey)> = self
            .records
            .iter()
            .filter(|(k, _)| !on_chain.contains(k))
            .map(|(k, r)| (r.used, *k))
            .collect();
        off_chain.sort_unstable_by_key(|(used, _)| std::cmp::Reverse(*used));
        for (_, key) in off_chain.iter().skip(self.chain.len().max(1)) {
            self.records.remove(key);
        }
        if let CacheStrategy::MemoryBudget { max_bytes } = self.strategy {
            let mut by_age: Vec<(u64, StepKey)> = self
                .records
                .iter()
                .filter(|(_, r)| r.state.is_some())
                .map(|(k, r)| (r.used, *k))
                .collect();
            by_age.sort_unstable_by_key(|(used, _)| *used);
            let mut oldest = by_age.into_iter();
            while self.memory_bytes() > max_bytes {
                let Some((_, key)) = oldest.next() else { break };
                if let Some(r) = self.records.get_mut(&key) {
                    r.state = None;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "cache_tests.rs"]
mod tests;
