//! Computation cache for healing analysis.
//!
//! Caches expensive computations (curvature fields, deviation measurements)
//! with entity-level invalidation: when an entity is modified, only its
//! cache and dependent caches are invalidated.
//!
//! # §18 hit/miss telemetry
//!
//! Every `get_*` call on [`HealingCache`] records a hit or a miss in a
//! [`CacheStats`](super::cache_stats::CacheStats) counter owned by the
//! cache. Call [`HealingCache::stats`] to retrieve the current hit rate.
//! §18 of `fillets_blends.md` requires cache hit rates to be queryable;
//! this module is the reference consumer of `CacheStats`. Additional
//! caches added to the engine should embed the same type.

use std::collections::HashMap;

use super::cache_stats::CacheStats;
use crate::topology::{EdgeHandle, FaceHandle};

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
#[path = "cache_tests.rs"]
mod tests;
