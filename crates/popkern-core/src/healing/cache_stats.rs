//! §18 cache hit/miss telemetry.
//!
//! [`CacheStats`] is a small counter block that any in-engine cache
//! embeds to satisfy §18's "cache hit rates — all queryable"
//! requirement. `HealingCache` is the reference consumer; other
//! caches (current or future) should embed the same type so the
//! metric stays consistent across popkern.

use std::cell::Cell;

/// §18 cache hit/miss telemetry for a single cache instance.
///
/// Records how many `get_*` calls found an entry in the cache (hits)
/// versus how many missed (misses). The hit rate is queryable via
/// [`hit_rate`](Self::hit_rate) and returns `None` before any queries
/// have been recorded — callers should treat an empty cache and an
/// unused cache the same way in UIs.
///
/// Uses `Cell<u64>` for interior mutability, so recording a hit or
/// miss does not require `&mut self` — compatible with `get_*`
/// accessors that take `&self`. The counters are not thread-safe;
/// if a cache ever moves to a `Sync` context, switch to
/// `AtomicU64` — the public API does not change.
#[derive(Debug, Default, Clone)]
pub struct CacheStats {
    hits: Cell<u64>,
    misses: Cell<u64>,
}

impl CacheStats {
    /// Create a fresh stats block with zero hits and zero misses.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of successful `get_*` lookups so far.
    pub fn hits(&self) -> u64 {
        self.hits.get()
    }

    /// Number of failed `get_*` lookups so far.
    pub fn misses(&self) -> u64 {
        self.misses.get()
    }

    /// Total queries (hits + misses).
    pub fn total_queries(&self) -> u64 {
        self.hits.get() + self.misses.get()
    }

    /// Hit rate in the range `0.0..=1.0`, or `None` if no queries
    /// have been recorded yet.
    pub fn hit_rate(&self) -> Option<f64> {
        let total = self.total_queries();
        if total == 0 {
            None
        } else {
            Some(self.hits.get() as f64 / total as f64)
        }
    }

    /// Reset both counters to zero. Useful when reusing a cache
    /// instance across operations and you want to track each
    /// operation's hit rate independently.
    pub fn reset(&self) {
        self.hits.set(0);
        self.misses.set(0);
    }

    pub(crate) fn record_hit(&self) {
        self.hits.set(self.hits.get() + 1);
    }

    pub(crate) fn record_miss(&self) {
        self.misses.set(self.misses.get() + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_stats_are_zero() {
        let stats = CacheStats::new();
        assert_eq!(stats.hits(), 0);
        assert_eq!(stats.misses(), 0);
        assert_eq!(stats.total_queries(), 0);
        // None (not 0.0) when nothing has been queried — lets a UI
        // distinguish "cache is cold" from "cache is broken".
        assert_eq!(stats.hit_rate(), None);
    }

    #[test]
    fn record_hits_and_misses() {
        let stats = CacheStats::new();
        stats.record_hit();
        stats.record_hit();
        stats.record_miss();
        assert_eq!(stats.hits(), 2);
        assert_eq!(stats.misses(), 1);
        assert!((stats.hit_rate().unwrap() - (2.0 / 3.0)).abs() < 1e-12);
    }

    #[test]
    fn reset_zeroes_counters() {
        let stats = CacheStats::new();
        stats.record_hit();
        stats.record_miss();
        stats.reset();
        assert_eq!(stats.total_queries(), 0);
        assert_eq!(stats.hit_rate(), None);
    }

    #[test]
    fn all_hits_reaches_one() {
        let stats = CacheStats::new();
        for _ in 0..10 {
            stats.record_hit();
        }
        assert_eq!(stats.hit_rate(), Some(1.0));
    }

    #[test]
    fn all_misses_reaches_zero() {
        let stats = CacheStats::new();
        for _ in 0..10 {
            stats.record_miss();
        }
        assert_eq!(stats.hit_rate(), Some(0.0));
    }
}
