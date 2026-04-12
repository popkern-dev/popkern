//! Computation cache for healing analysis.
//!
//! Caches expensive computations (curvature fields, deviation measurements)
//! with entity-level invalidation: when an entity is modified, only its
//! cache and dependent caches are invalidated.

use std::collections::HashMap;

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
}

impl HealingCache {
    /// Create a new empty cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// Get cached curvature for a face, or `None` if not cached.
    pub fn get_face_curvature(&self, face: FaceHandle) -> Option<&CachedCurvature> {
        self.face_curvature.get(&face)
    }

    /// Store curvature data for a face.
    pub fn set_face_curvature(&mut self, face: FaceHandle, data: CachedCurvature) {
        self.face_curvature.insert(face, data);
    }

    /// Get cached deviation for an edge, or `None` if not cached.
    pub fn get_edge_deviation(&self, edge: EdgeHandle) -> Option<&CachedDeviation> {
        self.edge_deviation.get(&edge)
    }

    /// Store deviation data for an edge.
    pub fn set_edge_deviation(&mut self, edge: EdgeHandle, data: CachedDeviation) {
        self.edge_deviation.insert(edge, data);
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
}
