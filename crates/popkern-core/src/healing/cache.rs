//! Computation cache for healing analysis.
//!
//! Caches expensive computations (curvature fields, deviation measurements)
//! with entity-level invalidation: when an entity is modified, only its
//! cache and dependent caches are invalidated.
//!
//! # §18 hit/miss telemetry
//!
//! Every `get_*` call on [`HealingCache`] records a hit or a miss in a
//! [`CacheStats`] counter owned by the cache. Call [`HealingCache::stats`]
//! to retrieve the current hit rate. §18 of `fillets_blends.md` requires
//! cache hit rates to be queryable; this module is the reference
//! implementation of that pattern for popkern. Additional caches added
//! to the engine should reuse [`CacheStats`] (it's re-exported at
//! `healing::CacheStats`) so callers have a single type to reason about.

use std::cell::Cell;
use std::collections::HashMap;

use crate::topology::{EdgeHandle, FaceHandle};

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

/// Cached curvature data for a face.
#[derive(Debug, Clone)]
pub struct CachedCurvature {
    /// Gaussian curvature at sampled points.
    pub gaussian_samples: Vec<(f64, f64, f64)>,
    /// Mean curvature at sampled points.
    pub mean_samples: Vec<(f64, f64, f64)>,
}

/// Cached deviation data for an edge.
#[derive(Debug, Clone)]
pub struct CachedDeviation {
    /// Maximum deviation from adjacent surface.
    pub max_deviation: f64,
    /// Mean deviation.
    pub mean_deviation: f64,
}

/// Cache for healing computations with entity-level invalidation.
///
/// Stores precomputed curvature and deviation data keyed by entity handle.
/// When an entity is modified, call [`invalidate_face`] or [`invalidate_edge`]
/// to remove stale entries and their dependents.
///
/// # Invalidation Cascade
///
/// Face invalidation removes the face's curvature entry and all deviation
/// entries for edges adjacent to that face (registered via
/// [`register_face_edges`]). Edge invalidation removes only the edge's
/// deviation entry.
#[derive(Debug, Clone, Default)]
pub struct HealingCache {
    /// Curvature data keyed by face handle.
    face_curvature: HashMap<FaceHandle, CachedCurvature>,
    /// Deviation data keyed by edge handle.
    edge_deviation: HashMap<EdgeHandle, CachedDeviation>,
    /// Tracks which edges are adjacent to which faces for invalidation.
    face_edges: HashMap<FaceHandle, Vec<EdgeHandle>>,
    /// §18 hit/miss telemetry — one counter block for all entries
    /// in this cache, regardless of whether they came from
    /// `get_face_curvature` or `get_edge_deviation`. Callers that
    /// need per-kind breakdown can split the cache into two.
    stats: CacheStats,
}

impl HealingCache {
    /// Create a new empty cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// Get cached curvature for a face, or `None` if not cached.
    /// Records a hit or miss in [`stats`](Self::stats).
    pub fn get_face_curvature(&self, face: FaceHandle) -> Option<&CachedCurvature> {
        let result = self.face_curvature.get(&face);
        if result.is_some() {
            self.stats.record_hit();
        } else {
            self.stats.record_miss();
        }
        result
    }

    /// Store curvature data for a face.
    pub fn set_face_curvature(&mut self, face: FaceHandle, data: CachedCurvature) {
        self.face_curvature.insert(face, data);
    }

    /// Get cached deviation for an edge, or `None` if not cached.
    /// Records a hit or miss in [`stats`](Self::stats).
    pub fn get_edge_deviation(&self, edge: EdgeHandle) -> Option<&CachedDeviation> {
        let result = self.edge_deviation.get(&edge);
        if result.is_some() {
            self.stats.record_hit();
        } else {
            self.stats.record_miss();
        }
        result
    }

    /// Store deviation data for an edge.
    pub fn set_edge_deviation(&mut self, edge: EdgeHandle, data: CachedDeviation) {
        self.edge_deviation.insert(edge, data);
    }

    /// §18 telemetry: access the hit/miss counters for this cache.
    /// Returned reference is valid for the cache's lifetime.
    pub fn stats(&self) -> &CacheStats {
        &self.stats
    }

    /// Register which edges are adjacent to a face (for cascade invalidation).
    ///
    /// When [`invalidate_face`] is called, all registered edge deviation
    /// entries for that face are also removed.
    pub fn register_face_edges(&mut self, face: FaceHandle, edges: Vec<EdgeHandle>) {
        self.face_edges.insert(face, edges);
    }

    /// Invalidate all cached data for a face and its adjacent edges.
    ///
    /// Removes the face's curvature entry. Then looks up the registered
    /// adjacent edges and removes their deviation entries as well.
    pub fn invalidate_face(&mut self, face: FaceHandle) {
        self.face_curvature.remove(&face);

        if let Some(edges) = self.face_edges.remove(&face) {
            for edge in &edges {
                self.edge_deviation.remove(edge);
            }
        }
    }

    /// Invalidate cached data for an edge.
    ///
    /// Removes only the edge's deviation entry. Does not affect
    /// adjacent face curvature data.
    pub fn invalidate_edge(&mut self, edge: EdgeHandle) {
        self.edge_deviation.remove(&edge);
    }

    /// Clear all cached data.
    pub fn clear(&mut self) {
        self.face_curvature.clear();
        self.edge_deviation.clear();
        self.face_edges.clear();
    }

    /// Number of cached entries (faces + edges).
    pub fn len(&self) -> usize {
        self.face_curvature.len() + self.edge_deviation.len()
    }

    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.face_curvature.is_empty() && self.edge_deviation.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use slotmap::SlotMap;

    /// Shared handle factory that produces distinct handles per call.
    struct HandleFactory {
        faces: SlotMap<FaceHandle, ()>,
        edges: SlotMap<EdgeHandle, ()>,
    }

    impl HandleFactory {
        fn new() -> Self {
            Self {
                faces: SlotMap::with_key(),
                edges: SlotMap::with_key(),
            }
        }

        fn face(&mut self) -> FaceHandle {
            self.faces.insert(())
        }

        fn edge(&mut self) -> EdgeHandle {
            self.edges.insert(())
        }
    }

    fn sample_curvature() -> CachedCurvature {
        CachedCurvature {
            gaussian_samples: vec![(0.0, 0.0, 1.0)],
            mean_samples: vec![(0.0, 0.0, 0.5)],
        }
    }

    fn sample_deviation() -> CachedDeviation {
        CachedDeviation {
            max_deviation: 0.001,
            mean_deviation: 0.0005,
        }
    }

    #[test]
    fn new_cache_is_empty() {
        let cache = HealingCache::new();
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);
    }

    #[test]
    fn store_and_retrieve_face_curvature() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let mut cache = HealingCache::new();
        assert!(cache.get_face_curvature(face).is_none());

        cache.set_face_curvature(face, sample_curvature());
        let cached = cache.get_face_curvature(face);
        assert!(cached.is_some());
        assert_eq!(cached.map(|c| c.gaussian_samples.len()), Some(1));
        assert_eq!(cache.len(), 1);
    }

    #[test]
    fn store_and_retrieve_edge_deviation() {
        let mut hf = HandleFactory::new();
        let edge = hf.edge();
        let mut cache = HealingCache::new();
        assert!(cache.get_edge_deviation(edge).is_none());

        cache.set_edge_deviation(edge, sample_deviation());
        let cached = cache.get_edge_deviation(edge);
        assert!(cached.is_some());
        assert!((cached.map(|c| c.max_deviation).unwrap_or(0.0) - 0.001).abs() < 1e-12);
        assert_eq!(cache.len(), 1);
    }

    #[test]
    fn invalidate_face_removes_curvature_and_dependent_edges() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let edge1 = hf.edge();
        let edge2 = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(face, sample_curvature());
        cache.set_edge_deviation(edge1, sample_deviation());
        cache.set_edge_deviation(edge2, sample_deviation());
        cache.register_face_edges(face, vec![edge1, edge2]);

        assert_eq!(cache.len(), 3);

        cache.invalidate_face(face);

        assert!(cache.get_face_curvature(face).is_none());
        assert!(cache.get_edge_deviation(edge1).is_none());
        assert!(cache.get_edge_deviation(edge2).is_none());
        assert!(cache.is_empty());
    }

    #[test]
    fn invalidate_edge_only_removes_edge_entry() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let edge = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(face, sample_curvature());
        cache.set_edge_deviation(edge, sample_deviation());

        cache.invalidate_edge(edge);

        assert!(cache.get_edge_deviation(edge).is_none());
        assert!(cache.get_face_curvature(face).is_some());
        assert_eq!(cache.len(), 1);
    }

    #[test]
    fn clear_empties_everything() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let edge = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(face, sample_curvature());
        cache.set_edge_deviation(edge, sample_deviation());
        cache.register_face_edges(face, vec![edge]);

        assert!(!cache.is_empty());

        cache.clear();
        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);
        assert!(cache.get_face_curvature(face).is_none());
        assert!(cache.get_edge_deviation(edge).is_none());
    }

    #[test]
    fn invalidate_face_with_no_registered_edges_only_removes_curvature() {
        let mut hf = HandleFactory::new();
        let face = hf.face();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(face, sample_curvature());

        cache.invalidate_face(face);

        assert!(cache.get_face_curvature(face).is_none());
        assert!(cache.is_empty());
    }

    #[test]
    fn invalidate_nonexistent_entries_is_no_op() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let edge = hf.edge();

        let mut cache = HealingCache::new();
        // Should not panic or error.
        cache.invalidate_face(face);
        cache.invalidate_edge(edge);
        assert!(cache.is_empty());
    }

    #[test]
    fn len_counts_faces_and_edges() {
        let mut hf = HandleFactory::new();
        let f1 = hf.face();
        let f2 = hf.face();
        let e1 = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(f1, sample_curvature());
        cache.set_face_curvature(f2, sample_curvature());
        cache.set_edge_deviation(e1, sample_deviation());
        assert_eq!(cache.len(), 3);
    }

    #[test]
    fn overwrite_face_curvature() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let mut cache = HealingCache::new();

        cache.set_face_curvature(face, sample_curvature());
        assert_eq!(cache.len(), 1);

        // Overwrite with new data.
        let updated = CachedCurvature {
            gaussian_samples: vec![(1.0, 1.0, 2.0), (2.0, 2.0, 3.0)],
            mean_samples: vec![(1.0, 1.0, 1.5)],
        };
        cache.set_face_curvature(face, updated);
        assert_eq!(cache.len(), 1); // Still one entry.
        assert_eq!(
            cache
                .get_face_curvature(face)
                .map(|c| c.gaussian_samples.len()),
            Some(2)
        );
    }

    #[test]
    fn overwrite_edge_deviation() {
        let mut hf = HandleFactory::new();
        let edge = hf.edge();
        let mut cache = HealingCache::new();

        cache.set_edge_deviation(edge, sample_deviation());
        assert_eq!(cache.len(), 1);

        let updated = CachedDeviation {
            max_deviation: 0.05,
            mean_deviation: 0.025,
        };
        cache.set_edge_deviation(edge, updated);
        assert_eq!(cache.len(), 1);
        assert!(
            (cache
                .get_edge_deviation(edge)
                .map(|c| c.max_deviation)
                .unwrap_or(0.0)
                - 0.05)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn register_face_edges_overwrites_previous() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let e1 = hf.edge();
        let e2 = hf.edge();
        let e3 = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(face, sample_curvature());
        cache.set_edge_deviation(e1, sample_deviation());
        cache.set_edge_deviation(e2, sample_deviation());
        cache.set_edge_deviation(e3, sample_deviation());

        // Register e1 and e2 as adjacent.
        cache.register_face_edges(face, vec![e1, e2]);
        // Overwrite: now only e3 is adjacent.
        cache.register_face_edges(face, vec![e3]);

        cache.invalidate_face(face);

        // e1 and e2 should still be cached (no longer registered).
        assert!(cache.get_edge_deviation(e1).is_some());
        assert!(cache.get_edge_deviation(e2).is_some());
        // e3 should be invalidated along with the face.
        assert!(cache.get_edge_deviation(e3).is_none());
    }

    #[test]
    fn multiple_faces_share_edge() {
        let mut hf = HandleFactory::new();
        let f1 = hf.face();
        let f2 = hf.face();
        let shared_edge = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(f1, sample_curvature());
        cache.set_face_curvature(f2, sample_curvature());
        cache.set_edge_deviation(shared_edge, sample_deviation());

        cache.register_face_edges(f1, vec![shared_edge]);
        cache.register_face_edges(f2, vec![shared_edge]);

        // Invalidate f1 — removes shared edge deviation.
        cache.invalidate_face(f1);
        assert!(cache.get_edge_deviation(shared_edge).is_none());
        // f2 curvature still cached.
        assert!(cache.get_face_curvature(f2).is_some());
    }

    #[test]
    fn clear_after_registrations() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let edge = hf.edge();

        let mut cache = HealingCache::new();
        cache.set_face_curvature(face, sample_curvature());
        cache.set_edge_deviation(edge, sample_deviation());
        cache.register_face_edges(face, vec![edge]);

        cache.clear();
        assert!(cache.is_empty());

        // Invalidating after clear should not panic.
        cache.invalidate_face(face);
        cache.invalidate_edge(edge);
        assert!(cache.is_empty());
    }

    #[test]
    fn default_is_same_as_new() {
        let d = HealingCache::default();
        let n = HealingCache::new();
        assert_eq!(d.len(), n.len());
        assert!(d.is_empty());
        assert!(n.is_empty());
    }

    // -----------------------------------------------------------------
    // §18 hit/miss telemetry
    // -----------------------------------------------------------------

    #[test]
    fn fresh_cache_stats_are_zero() {
        let cache = HealingCache::new();
        let stats = cache.stats();
        assert_eq!(stats.hits(), 0);
        assert_eq!(stats.misses(), 0);
        assert_eq!(stats.total_queries(), 0);
        // Hit rate is None (not 0.0) when nothing has been queried —
        // lets a UI distinguish "cache is cold" from "cache is broken".
        assert_eq!(stats.hit_rate(), None);
    }

    #[test]
    fn face_curvature_miss_then_hit_are_counted() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let mut cache = HealingCache::new();

        // Cold: miss.
        assert!(cache.get_face_curvature(face).is_none());
        assert_eq!(cache.stats().misses(), 1);
        assert_eq!(cache.stats().hits(), 0);
        assert_eq!(cache.stats().hit_rate(), Some(0.0));

        // Warm: populate and re-query.
        cache.set_face_curvature(face, sample_curvature());
        assert!(cache.get_face_curvature(face).is_some());
        assert_eq!(cache.stats().misses(), 1);
        assert_eq!(cache.stats().hits(), 1);
        assert_eq!(cache.stats().hit_rate(), Some(0.5));
        assert_eq!(cache.stats().total_queries(), 2);
    }

    #[test]
    fn edge_deviation_miss_then_hit_are_counted() {
        let mut hf = HandleFactory::new();
        let edge = hf.edge();
        let mut cache = HealingCache::new();

        // Cold.
        assert!(cache.get_edge_deviation(edge).is_none());
        // Warm.
        cache.set_edge_deviation(edge, sample_deviation());
        assert!(cache.get_edge_deviation(edge).is_some());
        // Third call — another hit.
        assert!(cache.get_edge_deviation(edge).is_some());

        assert_eq!(cache.stats().hits(), 2);
        assert_eq!(cache.stats().misses(), 1);
        assert!((cache.stats().hit_rate().unwrap() - (2.0 / 3.0)).abs() < 1e-12);
    }

    #[test]
    fn stats_aggregate_hits_across_both_cache_kinds() {
        // Both face-curvature and edge-deviation lookups increment
        // the same shared stats block — §18's "cache hit rate" is a
        // per-cache metric, not a per-bucket one.
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let edge = hf.edge();
        let mut cache = HealingCache::new();

        cache.set_face_curvature(face, sample_curvature());
        cache.set_edge_deviation(edge, sample_deviation());

        let _ = cache.get_face_curvature(face); // hit
        let _ = cache.get_edge_deviation(edge); // hit
        let _ = cache.get_face_curvature(hf.face()); // miss (different handle)
        let _ = cache.get_edge_deviation(hf.edge()); // miss

        assert_eq!(cache.stats().hits(), 2);
        assert_eq!(cache.stats().misses(), 2);
        assert_eq!(cache.stats().hit_rate(), Some(0.5));
    }

    #[test]
    fn stats_reset_zeroes_counters() {
        let mut hf = HandleFactory::new();
        let face = hf.face();
        let mut cache = HealingCache::new();

        cache.set_face_curvature(face, sample_curvature());
        let _ = cache.get_face_curvature(face); // hit
        let _ = cache.get_face_curvature(hf.face()); // miss

        assert_eq!(cache.stats().total_queries(), 2);
        cache.stats().reset();
        assert_eq!(cache.stats().hits(), 0);
        assert_eq!(cache.stats().misses(), 0);
        assert_eq!(cache.stats().hit_rate(), None);
    }
}
