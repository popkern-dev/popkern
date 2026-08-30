//! Session-state controls on [`FeatureTree`]: the replay cache and the
//! feature freeze point.

use crate::FeatureTree;
use crate::cache::CacheStrategy;

impl FeatureTree {
    // -- Cache management -----------------------------------------------------

    /// Replace the caching strategy.  All existing cached snapshots are evicted.
    pub fn set_cache_strategy(&mut self, strategy: CacheStrategy) {
        self.cache.set_strategy(strategy);
    }

    /// Evict all cached BRep snapshots.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }

    // -- Feature freeze -------------------------------------------------------

    /// Freeze all features up to (but not including) `index`.
    ///
    /// Frozen features are never re-evaluated during replay — the cached
    /// BRep at the freeze point is used as the starting state.
    pub fn freeze_to(&mut self, index: usize) {
        self.frozen_index = Some(index);
    }

    /// Remove the freeze point, allowing all features to be rebuilt.
    pub fn unfreeze(&mut self) {
        self.frozen_index = None;
    }

    /// Return the current freeze index, or `None` if nothing is frozen.
    pub fn frozen_index(&self) -> Option<usize> {
        self.frozen_index
    }
}
