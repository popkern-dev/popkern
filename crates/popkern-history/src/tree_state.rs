//! Session-state controls on [`FeatureTree`]: the replay cache and the
//! feature freeze point.

use crate::FeatureTree;
use crate::cache::CacheStrategy;

impl FeatureTree {
    // -- Cache management -----------------------------------------------------

    /// Replace the caching strategy.  All existing cached states are evicted.
    pub fn set_cache_strategy(&mut self, strategy: CacheStrategy) {
        self.cache.set_strategy(strategy);
    }

    /// Estimated bytes of every cached step state (bodies shared
    /// between steps count once).
    pub fn cache_memory_bytes(&self) -> usize {
        self.cache.memory_bytes()
    }

    /// Evict all cached step states.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }

    // -- Feature freeze -------------------------------------------------------

    /// Freeze all features up to (but not including) `index`.
    ///
    /// Frozen features are never re-evaluated during replay — the
    /// document state at the freeze point is used as the starting state,
    /// even after upstream edits. The state is pinned from the last
    /// replay when it is cached, else by the next replay that reaches it.
    pub fn freeze_to(&mut self, index: usize) {
        self.frozen_index = Some(index);
        match self.cached_state_at(index).filter(|_| index > 0).cloned() {
            Some(state) => self.cache.set_frozen(index, state),
            None => self.cache.clear_frozen(),
        }
    }

    /// Remove the freeze point, allowing all features to be rebuilt.
    pub fn unfreeze(&mut self) {
        self.frozen_index = None;
        self.cache.clear_frozen();
    }

    /// Return the current freeze index, or `None` if nothing is frozen.
    pub fn frozen_index(&self) -> Option<usize> {
        self.frozen_index
    }
}
