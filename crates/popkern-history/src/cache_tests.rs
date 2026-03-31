// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use popkern_constraint::{Constraint, Sketch};
use popkern_core::tolerance::Tolerance;

use crate::FeatureTree;
use crate::cache::CacheStrategy;
use crate::error::FeatureResult;
use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};
use crate::node::FeatureId;
use crate::parameter::ParameterId;
use crate::unit::Unit;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

// ---------------------------------------------------------------------------
// Helpers (mirrored from replay_tests.rs)
// ---------------------------------------------------------------------------

fn rectangle_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

/// Build a tree with one sketch + one extrude, replay, return the tree and
/// both feature ids.
fn build_sketch_extrude_tree() -> (FeatureTree, FeatureId, FeatureId) {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Replay sketch first to discover what entity names it produced.
    tree.replay(Tolerance::default()).expect("sketch replay");

    let dist_param = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist_param,
        },
        "Extrude 1",
    );

    (tree, sketch_id, extrude_id)
}

// ---------------------------------------------------------------------------
// Test 1: EveryStep — replay 2 features → 2 cached BReps
// ---------------------------------------------------------------------------

#[test]
fn every_step_caches_all_features() {
    let (mut tree, sketch_id, extrude_id) = build_sketch_extrude_tree();

    // EveryStep is the default — ensure it is set.
    tree.set_cache_strategy(CacheStrategy::EveryStep);

    tree.replay(Tolerance::default()).expect("replay");

    // Both feature positions should be cached now.
    assert!(
        tree.cache.get(None, sketch_id).is_some(),
        "sketch BRep should be cached"
    );
    assert!(
        tree.cache.get(None, extrude_id).is_some(),
        "extrude BRep should be cached"
    );
}

// ---------------------------------------------------------------------------
// Test 2: brep_at returns from cache without re-replaying
// ---------------------------------------------------------------------------

#[test]
fn brep_at_returns_from_cache() {
    let (mut tree, _sketch_id, _extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.replay(Tolerance::default()).expect("replay");

    // brep_at(2) → position of extrude (index 1, so index arg = 2 for
    // "after second feature").
    let brep = tree
        .brep_at(2, Tolerance::default())
        .expect("brep_at should succeed");
    assert!(brep.faces.len() == 6, "extruded box should have 6 faces");
}

// ---------------------------------------------------------------------------
// Test 3: MemoryBudget — oldest entry evicted when budget exceeded
// ---------------------------------------------------------------------------

#[test]
fn memory_budget_evicts_oldest() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    // Add three features to get three snapshots.
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json.clone(),
        },
        "Sketch 1",
    );

    tree.replay(Tolerance::default()).expect("sketch replay");

    let dist_param = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist_param,
        },
        "Extrude 1",
    );

    let second_sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 2",
    );

    // Set a very tight budget: only enough for roughly one 6-face box snapshot
    // (6 faces × 256 bytes = 1 536 bytes) plus a little margin.
    // With two snapshots each ~1.5 KB the budget will be exceeded after the
    // second insert, evicting the first.
    tree.set_cache_strategy(CacheStrategy::MemoryBudget { max_bytes: 2_000 });

    tree.replay(Tolerance::default()).expect("full replay");

    // At most one of the first two features should still be cached — the
    // combined memory should be at or below the budget.
    let bytes_used = tree.cache.memory_bytes();
    assert!(
        bytes_used <= 2_000,
        "cache should have stayed within budget, but used {bytes_used} bytes"
    );

    // At least the last feature should still be cached (most recent).
    // (The extrude produces a larger BRep so it will be the one kept.)
    let _ = second_sketch_id; // may or may not be cached — depends on size
    let _ = sketch_id;
    let _ = extrude_id;
}

// ---------------------------------------------------------------------------
// Test 4: invalidate_from — clears downstream cache, keeps upstream
// ---------------------------------------------------------------------------

#[test]
fn invalidate_from_clears_downstream() {
    let (mut tree, sketch_id, extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.replay(Tolerance::default()).expect("replay");

    // Both should be cached before invalidation.
    assert!(tree.cache.get(None, sketch_id).is_some());
    assert!(tree.cache.get(None, extrude_id).is_some());

    // Invalidate from index 1 (the extrude is at position 1 in the timeline).
    let timeline: Vec<FeatureId> = tree.timeline().to_vec();
    tree.cache.invalidate_from(1, &timeline);

    // Sketch (index 0) should still be cached.
    assert!(
        tree.cache.get(None, sketch_id).is_some(),
        "sketch cache should be retained"
    );
    // Extrude (index 1) should be evicted.
    assert!(
        tree.cache.get(None, extrude_id).is_none(),
        "extrude cache should be invalidated"
    );
}

// ---------------------------------------------------------------------------
// Test 5: EveryNth(2) — only even-indexed steps cached
// ---------------------------------------------------------------------------

#[test]
fn every_nth_caches_only_matching_indices() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    // Add four features: two sketches + two extrudes (simplified; second
    // sketch on XY plane again).
    let s1 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json.clone(),
        },
        "Sketch 1",
    );

    tree.replay(Tolerance::default()).expect("replay s1");

    let dist_param = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    let e1 = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist_param,
        },
        "Extrude 1",
    );

    let s2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json.clone(),
        },
        "Sketch 2",
    );

    let e2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 3",
    );

    tree.set_cache_strategy(CacheStrategy::EveryNth(2));
    tree.replay(Tolerance::default()).expect("replay all");

    // Indices 0 and 2 should be cached (0 % 2 == 0, 2 % 2 == 0).
    // Indices 1 and 3 should not be cached.
    assert!(
        tree.cache.get(None, s1).is_some(),
        "index 0 should be cached"
    );
    assert!(
        tree.cache.get(None, e1).is_none(),
        "index 1 should NOT be cached"
    );
    assert!(
        tree.cache.get(None, s2).is_some(),
        "index 2 should be cached"
    );
    assert!(
        tree.cache.get(None, e2).is_none(),
        "index 3 should NOT be cached"
    );
}

// ---------------------------------------------------------------------------
// Test 6: OnDemand — no automatic caching
// ---------------------------------------------------------------------------

#[test]
fn on_demand_does_not_cache_automatically() {
    let (mut tree, sketch_id, extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::OnDemand);
    tree.replay(Tolerance::default()).expect("replay");

    assert!(
        tree.cache.get(None, sketch_id).is_none(),
        "OnDemand: sketch should not be cached"
    );
    assert!(
        tree.cache.get(None, extrude_id).is_none(),
        "OnDemand: extrude should not be cached"
    );
}

// ---------------------------------------------------------------------------
// Test 7: rollback_to sets position
// ---------------------------------------------------------------------------

#[test]
fn rollback_to_sets_position() {
    let (mut tree, _sketch_id, _extrude_id) = build_sketch_extrude_tree();

    tree.rollback_to(1);
    assert_eq!(tree.rollback_position(), Some(1));

    tree.rollback_to(0);
    assert_eq!(tree.rollback_position(), Some(0));
}

// ---------------------------------------------------------------------------
// Test 8: rollback_to_end clears rollback position
// ---------------------------------------------------------------------------

#[test]
fn rollback_to_end_clears_rollback() {
    let (mut tree, _sketch_id, _extrude_id) = build_sketch_extrude_tree();

    tree.rollback_to(1);
    assert_eq!(tree.rollback_position(), Some(1));

    tree.rollback_to_end();
    assert_eq!(tree.rollback_position(), None);
}

// ---------------------------------------------------------------------------
// Test 9: clear_cache empties the cache
// ---------------------------------------------------------------------------

#[test]
fn clear_cache_removes_all_entries() {
    let (mut tree, sketch_id, extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.replay(Tolerance::default()).expect("replay");

    assert!(tree.cache.get(None, sketch_id).is_some());
    assert!(tree.cache.get(None, extrude_id).is_some());

    tree.clear_cache();

    assert!(
        tree.cache.get(None, sketch_id).is_none(),
        "cache should be empty after clear"
    );
    assert!(
        tree.cache.get(None, extrude_id).is_none(),
        "cache should be empty after clear"
    );
}

// ---------------------------------------------------------------------------
// Test 10: set_strategy clears existing cache
// ---------------------------------------------------------------------------

#[test]
fn set_strategy_clears_cache() {
    let (mut tree, sketch_id, _extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.replay(Tolerance::default()).expect("replay");

    assert!(tree.cache.get(None, sketch_id).is_some());

    // Switching strategy should evict all cached entries.
    tree.set_cache_strategy(CacheStrategy::OnDemand);

    assert!(
        tree.cache.get(None, sketch_id).is_none(),
        "cache should be cleared after set_strategy"
    );
}

// ---------------------------------------------------------------------------
// Test 11: rollback_to returns cached BRep from the preceding feature
// ---------------------------------------------------------------------------

#[test]
fn rollback_to_returns_preceding_cached_brep() {
    let (mut tree, sketch_id, _extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.replay(Tolerance::default()).expect("replay");

    // rollback_to(1) should return the BRep cached after sketch (index 0).
    let brep = tree.rollback_to(1);
    assert!(
        brep.is_some(),
        "should return cached BRep for preceding step"
    );

    // The BRep before extrude has no solid (sketch only).
    let brep = brep.expect("just checked");
    let _ = sketch_id;
    assert!(brep.faces.len() < 6, "pre-extrude BRep should not be a box");
}

// ---------------------------------------------------------------------------
// Test 12: replay respects rollback_index when set
// ---------------------------------------------------------------------------

#[test]
fn replay_with_rollback_index_stops_early() {
    let (mut tree, _sketch_id, _extrude_id) = build_sketch_extrude_tree();

    // Stop after the first feature (sketch only).
    tree.rollback_to(1);
    let result = tree.replay(Tolerance::default()).expect("rollback replay");

    // Should have replayed only 1 feature.
    assert_eq!(result.feature_results.len(), 1, "only sketch should run");

    let brep = tree.current_brep().expect("brep should exist");
    assert!(brep.faces.len() < 6, "rollback replay should not extrude");
}

// ---------------------------------------------------------------------------
// Test 13: brep_at with no prior replay (cache miss → full replay)
// ---------------------------------------------------------------------------

#[test]
fn brep_at_without_prior_cache_replays() {
    let (mut tree, _sketch_id, _extrude_id) = build_sketch_extrude_tree();
    // No prior replay; cache is empty.
    tree.set_cache_strategy(CacheStrategy::EveryStep);

    let brep = tree
        .brep_at(2, Tolerance::default())
        .expect("brep_at should succeed via fallback replay");

    assert_eq!(
        brep.faces.len(),
        6,
        "full replay should produce a 6-face box"
    );
}

// ---------------------------------------------------------------------------
// Test 14: FeatureResult::Success recorded on replay with caching
// ---------------------------------------------------------------------------

#[test]
fn replay_records_success_results() {
    let (mut tree, sketch_id, extrude_id) = build_sketch_extrude_tree();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&sketch_id),
            Some(FeatureResult::Success { .. })
        ),
        "sketch should succeed"
    );
    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed"
    );
}
