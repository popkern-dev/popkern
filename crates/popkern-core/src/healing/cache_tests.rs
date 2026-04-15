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
