//! BRep intermediate-state cache for the parametric feature tree.
//!
//! [`BRepCache`] stores cloned [`BRep`] snapshots keyed by [`FeatureId`],
//! allowing [`crate::FeatureTree::brep_at`] and
//! [`crate::FeatureTree::replay_from`] to resume replays from the nearest
//! already-computed ancestor rather than from scratch.

use std::collections::{HashMap, VecDeque};

use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use crate::node::FeatureId;

// ---------------------------------------------------------------------------
// CacheStrategy
// ---------------------------------------------------------------------------

/// Controls which intermediate BRep states are retained in memory.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum CacheStrategy {
    /// Cache the BRep after every feature — instant scrubbing at the cost of
    /// higher memory usage.
    EveryStep,
    /// Cache the BRep after every *N*th feature — balanced trade-off between
    /// memory and replay latency.
    EveryNth(usize),
    /// Never cache automatically; snapshots can still be stored explicitly via
    /// [`BRepCache::store_explicit`].
    OnDemand,
    /// Cache until the estimated memory footprint exceeds `max_bytes`, then
    /// evict the oldest entry (LRU eviction by insertion order).
    MemoryBudget {
        /// Maximum estimated cache size in bytes before eviction begins.
        max_bytes: usize,
    },
}

// ---------------------------------------------------------------------------
// BRepCache
// ---------------------------------------------------------------------------

/// Cache of intermediate BRep states, keyed by the [`FeatureId`] of the last
/// feature that was applied to produce each snapshot.
#[derive(Debug, Clone)]
pub(crate) struct BRepCache {
    cache: HashMap<FeatureId, BRep>,
    strategy: CacheStrategy,
    /// Tracks insertion order, oldest entry first, for LRU eviction.
    insertion_order: VecDeque<FeatureId>,
}

impl Default for BRepCache {
    fn default() -> Self {
        Self::new(CacheStrategy::EveryStep)
    }
}

#[allow(dead_code)] // Several methods are used only in tests or future tasks.
impl BRepCache {
    /// Create a new, empty cache with the given strategy.
    pub(crate) fn new(strategy: CacheStrategy) -> Self {
        Self {
            cache: HashMap::new(),
            strategy,
            insertion_order: VecDeque::new(),
        }
    }

    /// Return a reference to the cached BRep produced after `id`, or `None`
    /// if no snapshot exists for that feature.
    pub(crate) fn get(&self, id: FeatureId) -> Option<&BRep> {
        self.cache.get(&id)
    }

    /// Store a snapshot according to the current [`CacheStrategy`].
    ///
    /// `timeline_index` is the zero-based position of this feature in the
    /// timeline; it is used to decide whether to cache for `EveryNth`.
    pub(crate) fn store(&mut self, id: FeatureId, brep: BRep, timeline_index: usize) {
        let should_store = match self.strategy {
            CacheStrategy::EveryStep => true,
            CacheStrategy::EveryNth(n) => n > 0 && timeline_index.is_multiple_of(n),
            CacheStrategy::OnDemand => false,
            CacheStrategy::MemoryBudget { .. } => true,
        };

        if !should_store {
            return;
        }

        self.insert(id, brep);

        // For MemoryBudget, evict oldest entries until we are under budget.
        if let CacheStrategy::MemoryBudget { max_bytes } = self.strategy {
            while self.memory_bytes() > max_bytes && self.insertion_order.len() > 1 {
                if let Some(oldest) = self.insertion_order.pop_front() {
                    self.cache.remove(&oldest);
                }
            }
        }
    }

    /// Store a snapshot unconditionally, regardless of the current strategy.
    ///
    /// Use this to manually pin a snapshot when using [`CacheStrategy::OnDemand`].
    pub(crate) fn store_explicit(&mut self, id: FeatureId, brep: BRep) {
        self.insert(id, brep);
    }

    /// Remove cached BReps for every feature at or after `index` in
    /// `timeline`.
    ///
    /// This is called after editing a feature so that stale downstream
    /// snapshots are evicted, while snapshots for unaffected ancestors are
    /// retained.
    pub(crate) fn invalidate_from(&mut self, index: usize, timeline: &[FeatureId]) {
        let to_remove: Vec<FeatureId> = timeline[index..].to_vec();
        for id in &to_remove {
            self.cache.remove(id);
        }
        self.insertion_order.retain(|id| !to_remove.contains(id));
    }

    /// Remove all cached snapshots.
    pub(crate) fn clear(&mut self) {
        self.cache.clear();
        self.insertion_order.clear();
    }

    /// Estimate the memory occupied by all cached BReps, in bytes.
    ///
    /// Uses a simple heuristic: each entity contributes a fixed number of
    /// bytes (vertices: 64, edges: 96, faces: 256).  The result is suitable
    /// for budget enforcement but is not exact.
    pub(crate) fn memory_bytes(&self) -> usize {
        self.cache.values().map(estimate_brep_bytes).sum()
    }

    /// Return the current caching strategy.
    pub(crate) fn strategy(&self) -> CacheStrategy {
        self.strategy
    }

    /// Replace the current caching strategy and clear all cached snapshots.
    pub(crate) fn set_strategy(&mut self, strategy: CacheStrategy) {
        self.strategy = strategy;
        self.clear();
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Insert an entry, maintaining `insertion_order`.
    fn insert(&mut self, id: FeatureId, brep: BRep) {
        if !self.cache.contains_key(&id) {
            self.insertion_order.push_back(id);
        }
        self.cache.insert(id, brep);
    }
}

/// Estimate the heap memory consumed by a single BRep snapshot, in bytes.
///
/// The heuristic is: vertices × 64 + edges × 96 + faces × 256.
fn estimate_brep_bytes(brep: &BRep) -> usize {
    brep.vertices.len() * 64 + brep.edges.len() * 96 + brep.faces.len() * 256
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "cache_tests.rs"]
mod tests;
