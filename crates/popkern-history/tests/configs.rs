//! Integration tests for the configuration / design-table system.
//!
//! Exercises the full pipeline through the public API: config lifecycle,
//! derived inheritance, suppression, undo/redo, CSV round-trip, save/load,
//! delete cleanup, and config-aware caching.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

// ---------------------------------------------------------------------------
// Helpers
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
    sk.solve().expect("sketch solve");
    serde_json::to_string(&sk).expect("sketch serialize")
}

fn sketch_face_name(op_id: u64) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    }
}

/// Sketch + extrude tree. Returns `(tree, sketch_id, extrude_id, dist_param)`.
fn build_sketch_extrude(distance_m: f64) -> (FeatureTree, FeatureId, FeatureId, ParameterId) {
    let mut tree = FeatureTree::new();
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let dist = tree
        .parameters_mut()
        .add_global("dist", distance_m, Unit::Meters)
        .expect("dist");
    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );
    (tree, sketch_id, extrude_id, dist)
}

/// Maximum Z across all vertices of the current BRep.
fn max_z(tree: &FeatureTree) -> f64 {
    tree.current_brep()
        .expect("brep")
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Find a config by name, panicking if not found.
fn find_config(tree: &FeatureTree, name: &str) -> ConfigId {
    tree.configs()
        .into_iter()
        .find(|&id| tree.config(id).map(|c| c.name == name).unwrap_or(false))
        .unwrap_or_else(|| panic!("config '{name}' not found"))
}

// ---------------------------------------------------------------------------
// Test 1: Full config pipeline — switch configs, replay, verify different BReps
// ---------------------------------------------------------------------------

#[test]
fn full_config_pipeline_different_brep_states() {
    let (mut tree, _sk, _ex, dist) = build_sketch_extrude(1.0);

    let small = tree.add_config("Small", "Small variant").expect("Small");
    tree.set_config_override(small, dist, ParamOverride::Literal(0.5))
        .expect("Small override");

    let large = tree.add_config("Large", "Large variant").expect("Large");
    tree.set_config_override(large, dist, ParamOverride::Literal(2.0))
        .expect("Large override");

    tree.set_active_config(None);
    tree.replay(Tolerance::default()).expect("replay default");
    let z_default = max_z(&tree);
    assert!(
        (z_default - 1.0).abs() < 0.01,
        "default height: {z_default}"
    );
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);

    tree.set_active_config(Some(small));
    tree.replay(Tolerance::default()).expect("replay Small");
    let z_small = max_z(&tree);
    assert!((z_small - 0.5).abs() < 0.01, "Small height: {z_small}");

    tree.set_active_config(Some(large));
    tree.replay(Tolerance::default()).expect("replay Large");
    let z_large = max_z(&tree);
    assert!((z_large - 2.0).abs() < 0.01, "Large height: {z_large}");

    assert!(
        z_small < z_default && z_default < z_large,
        "ordering: {z_small} < {z_default} < {z_large}"
    );
}

// ---------------------------------------------------------------------------
// Test 2: Derived configurations — inheritance chain
// ---------------------------------------------------------------------------

#[test]
fn derived_config_inherits_parent_overrides() {
    let mut tree = FeatureTree::new();
    let width = tree
        .parameters_mut()
        .add_global("width", 1.0, Unit::Meters)
        .expect("width");
    let depth = tree
        .parameters_mut()
        .add_global("depth", 1.0, Unit::Meters)
        .expect("depth");

    let base = tree.add_config("Base", "").expect("Base");
    tree.set_config_override(base, width, ParamOverride::Literal(2.0))
        .expect("width override");

    let derived = tree
        .add_derived_config("Base-v1", "", base)
        .expect("Base-v1");
    tree.set_config_override(derived, depth, ParamOverride::Literal(3.0))
        .expect("depth override");

    // Own overrides.
    let base_cfg = tree.config(base).expect("base");
    assert!(base_cfg.overrides.contains_key(&width), "Base: width");
    assert!(!base_cfg.overrides.contains_key(&depth), "Base: no depth");

    let derived_cfg = tree.config(derived).expect("derived");
    assert_eq!(derived_cfg.parent, Some(base), "parent link");
    assert!(derived_cfg.overrides.contains_key(&depth), "Base-v1: depth");
    assert!(
        !derived_cfg.overrides.contains_key(&width),
        "Base-v1: width inherited"
    );

    // Effective evaluation — Base-v1 inherits width from Base.
    let values = tree.evaluate_config(Some(derived));
    let w = values.get(&width).copied().unwrap_or(0.0);
    let d = values.get(&depth).copied().unwrap_or(0.0);
    assert!((w - 2.0).abs() < 1e-10, "inherited width: {w}");
    assert!((d - 3.0).abs() < 1e-10, "own depth: {d}");
}

// ---------------------------------------------------------------------------
// Test 3: Config suppression — extrude suppressed in one config
// ---------------------------------------------------------------------------

#[test]
fn config_suppression_skips_extrude() {
    let (mut tree, _sk, extrude_id, _dist) = build_sketch_extrude(1.0);

    tree.replay(Tolerance::default()).expect("replay default");
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);

    let flat = tree.add_config("FlatOnly", "").expect("FlatOnly");
    tree.set_config_suppression(flat, extrude_id, true);

    tree.set_active_config(Some(flat));
    tree.replay(Tolerance::default()).expect("replay FlatOnly");

    let faces = tree.current_brep().unwrap().faces.len();
    assert!(
        faces < 6,
        "suppressed extrude → fewer than 6 faces, got {faces}"
    );
    assert!(
        tree.effective_suppressed(extrude_id),
        "effective_suppressed = true"
    );

    tree.set_active_config(None);
    assert!(
        !tree.effective_suppressed(extrude_id),
        "effective_suppressed = false for default"
    );
}

// ---------------------------------------------------------------------------
// Test 4: Config undo / redo
// ---------------------------------------------------------------------------

#[test]
fn config_add_undo_redo() {
    let mut tree = FeatureTree::new();
    let id = tree.add_config("Undo-me", "").expect("add");
    assert_eq!(tree.configs().len(), 1);

    tree.undo();
    assert!(tree.configs().is_empty(), "gone after undo");
    assert!(tree.config(id).is_none());

    tree.redo();
    assert_eq!(tree.configs().len(), 1, "back after redo");
    assert!(tree.config(id).is_some());
}

#[test]
fn config_override_undo_redo() {
    let mut tree = FeatureTree::new();
    let param = tree
        .parameters_mut()
        .add_global("h", 1.0, Unit::Meters)
        .expect("h");
    let cfg = tree.add_config("Cfg", "").expect("Cfg");

    tree.set_config_override(cfg, param, ParamOverride::Literal(5.0))
        .expect("override");
    assert_eq!(tree.config_overrides(cfg).unwrap().len(), 1);

    tree.undo();
    assert!(
        tree.config_overrides(cfg).unwrap().is_empty(),
        "cleared by undo"
    );

    tree.redo();
    assert_eq!(
        tree.config_overrides(cfg).unwrap().len(),
        1,
        "restored by redo"
    );
}

// ---------------------------------------------------------------------------
// Test 5: CSV round-trip
// ---------------------------------------------------------------------------

#[test]
fn csv_round_trip() {
    let mut tree = FeatureTree::new();
    let width = tree
        .parameters_mut()
        .add_global("width", 0.05, Unit::Meters)
        .expect("width");

    let small = tree.add_config("Small", "Small variant").expect("Small");
    tree.set_config_override(small, width, ParamOverride::Literal(0.025))
        .expect("Small override");

    let large = tree.add_config("Large", "Large variant").expect("Large");
    tree.set_config_override(large, width, ParamOverride::Literal(0.1))
        .expect("Large override");

    let csv = tree.export_configs_csv();
    assert!(csv.contains("Small") && csv.contains("Large"), "CSV: {csv}");

    let mut tree2 = FeatureTree::new();
    let width2 = tree2
        .parameters_mut()
        .add_global("width", 0.05, Unit::Meters)
        .expect("width2");

    let result = tree2.import_configs_csv(&csv).expect("import");
    assert!(
        result.warnings.is_empty(),
        "warnings: {:?}",
        result.warnings
    );
    assert_eq!(result.created.len(), 2);

    let small2 = find_config(&tree2, "Small");
    match tree2
        .config_overrides(small2)
        .expect("overrides")
        .get(&width2)
    {
        Some(ParamOverride::Literal(v)) => assert!((*v - 0.025).abs() < 1e-10, "Small: {v}"),
        other => panic!("expected Literal, got {other:?}"),
    }

    let large2 = find_config(&tree2, "Large");
    match tree2
        .config_overrides(large2)
        .expect("overrides")
        .get(&width2)
    {
        Some(ParamOverride::Literal(v)) => assert!((*v - 0.1).abs() < 1e-10, "Large: {v}"),
        other => panic!("expected Literal, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Test 6: Save / load preserves configs and active_config
// ---------------------------------------------------------------------------

#[test]
fn save_load_preserves_configs_and_active_config() {
    let (mut tree, _sk, _ex, dist) = build_sketch_extrude(1.0);

    let small = tree.add_config("Small", "").expect("Small");
    tree.set_config_override(small, dist, ParamOverride::Literal(0.5))
        .expect("Small override");

    let large = tree.add_config("Large", "").expect("Large");
    tree.set_active_config(Some(large));

    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");

    assert_eq!(loaded.configs().len(), 2, "two configs preserved");

    let small_l = find_config(&loaded, "Small");
    let large_l = find_config(&loaded, "Large");

    let overrides_l = loaded.config_overrides(small_l).expect("Small overrides");
    assert_eq!(overrides_l.len(), 1, "Small: one override preserved");

    assert_eq!(
        loaded.active_config(),
        Some(large_l),
        "active_config preserved"
    );
}

// ---------------------------------------------------------------------------
// Test 7: Delete cascade cleans up config suppression overrides
// ---------------------------------------------------------------------------

#[test]
fn delete_cascade_cleans_config_suppression_override() {
    let (mut tree, _sk, extrude_id, _dist) = build_sketch_extrude(1.0);

    let cfg = tree.add_config("Cfg", "").expect("Cfg");
    tree.set_config_suppression(cfg, extrude_id, true);
    assert_eq!(tree.config_suppressions(cfg).unwrap().len(), 1);

    tree.delete_feature(extrude_id, DeleteBehavior::Cascade)
        .expect("delete");

    assert!(
        tree.config_suppressions(cfg).unwrap().is_empty(),
        "suppression cleaned up after cascade delete"
    );
}

// ---------------------------------------------------------------------------
// Test 8: Config-aware caching — brep_at uses config-specific cache slots
// ---------------------------------------------------------------------------
//
// The cache is keyed by `(active_config, feature_id)`.  Switching configs
// never serves stale BReps from a different set of overrides.  We verify
// through `brep_at`, which returns the cached BRep when available.

#[test]
fn config_aware_caching_brep_at_uses_config_specific_slots() {
    let (mut tree, _sk, _ex, dist) = build_sketch_extrude(1.0);
    let tol = Tolerance::default();

    // Replay default — populates cache under `None`.
    tree.set_active_config(None);
    tree.replay(tol).expect("replay default");

    let (z_default, faces_default) = {
        let brep = tree.brep_at(2, tol).expect("brep_at(2) default");
        let z = brep
            .vertices
            .iter()
            .map(|(_, v)| v.position.z)
            .fold(f64::NEG_INFINITY, f64::max);
        (z, brep.faces.len())
    };
    assert!(
        (z_default - 1.0).abs() < 0.01,
        "default height: {z_default}"
    );

    // Replay under a named config — populates cache under `Some(cfg)`.
    let cfg = tree.add_config("Alt", "").expect("Alt");
    tree.set_config_override(cfg, dist, ParamOverride::Literal(0.5))
        .expect("Alt override");
    tree.set_active_config(Some(cfg));
    tree.replay(tol).expect("replay Alt");

    let (z_alt, faces_alt) = {
        let brep = tree.brep_at(2, tol).expect("brep_at(2) Alt");
        let z = brep
            .vertices
            .iter()
            .map(|(_, v)| v.position.z)
            .fold(f64::NEG_INFINITY, f64::max);
        (z, brep.faces.len())
    };
    assert!((z_alt - 0.5).abs() < 0.01, "Alt height: {z_alt}");

    // Different heights confirm separate cache slots per config.
    assert!(
        (z_default - z_alt).abs() > 0.1,
        "cache slots differ: {z_default} vs {z_alt}"
    );
    assert_eq!(faces_default, faces_alt, "both configs: 6-face box");
}

/// CSV import is a single undo step — one undo reverses the entire import.
#[test]
fn csv_import_is_single_undo_step() {
    let (mut tree, _, _, _) = build_sketch_extrude(0.1);

    // Export, then import 2 configs via CSV.
    let _cfg1 = tree.add_config("Alpha", "first").expect("cfg1");
    let _cfg2 = tree.add_config("Beta", "second").expect("cfg2");
    let csv = tree.export_configs_csv();

    // Start fresh tree, import the CSV.
    let (mut tree2, _, _, _) = build_sketch_extrude(0.1);
    assert!(tree2.configs().is_empty());

    let import_result = tree2.import_configs_csv(&csv).expect("import");
    assert_eq!(import_result.created.len(), 2);
    assert_eq!(tree2.configs().len(), 2);

    // Single undo should reverse the entire import.
    assert!(tree2.undo());
    assert_eq!(
        tree2.configs().len(),
        0,
        "one undo cleared both imported configs"
    );

    // Single redo restores both.
    assert!(tree2.redo());
    assert_eq!(tree2.configs().len(), 2, "one redo restored both configs");
}
