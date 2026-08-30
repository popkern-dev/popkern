//! Integration tests for feature patterns in `popkern-history` (part 2 of 2).
//!
//! Tests 5–9: OperationId persistence, instance_op_ids serde, CircularPattern,
//! Mirror, and geometry_pattern flag.
//!
//! See also: `pattern_integration.rs` for tests 1–4.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_history::*;

// ---------------------------------------------------------------------------
// Helpers (duplicated from pattern_integration.rs — each integration test
// binary is a standalone crate and cannot share non-pub items across files)
// ---------------------------------------------------------------------------

/// Solved 1×1 rectangle sketch serialised to JSON.
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

/// Construct a generated EntityName referencing operation `op_id`, index 0.
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

/// Add a global parameter and return its id.
fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add_param should succeed")
}

/// Build a Sketch + Extrude(NewBody) tree.  Returns (tree, sketch_id, extrude_id).
fn build_sketch_extrude(distance: f64) -> (FeatureTree, FeatureId, FeatureId) {
    let mut tree = FeatureTree::new();
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let dist_param = add_param(&mut tree, "_extrude1_distance", distance, Unit::Meters);
    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind {
                distance: dist_param,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );
    (tree, sketch_id, extrude_id)
}

// ---------------------------------------------------------------------------
// Test 5: instance_op_ids default to empty and survive a save/load round-trip
//
// The `instance_op_ids` field uses `#[serde(skip_serializing_if = "is_empty")]`
// so a freshly created pattern (before any replay) serialises without that
// key.  On load the field defaults to an empty HashMap.
// ---------------------------------------------------------------------------

#[test]
fn pattern_instance_op_ids_default_empty_and_survive_round_trip() {
    let (mut tree, _sketch_id, extrude_id) = build_sketch_extrude(1.0);

    let count_param = add_param(&mut tree, "_pat_count", 3.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_pat_total", 20.0, Unit::Meters);

    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![extrude_id],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPatternOps",
    );

    // Before replay: instance_op_ids must be empty.
    let node = tree.feature(pattern_id).expect("pattern node must exist");
    assert!(
        node.instance_op_ids.is_empty(),
        "instance_op_ids must be empty before replay"
    );

    // Round-trip: empty map survives save/load.
    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    let loaded_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::LinearPattern { .. })
            )
        })
        .expect("LinearPattern must survive round-trip");

    let loaded_node = loaded.feature(loaded_id).expect("loaded node must exist");
    assert!(
        loaded_node.instance_op_ids.is_empty(),
        "instance_op_ids must remain empty after round-trip (no replay performed)"
    );
}

// ---------------------------------------------------------------------------
// Test 6: OperationId counter advances and does not collide across patterns
//
// Verifies that two pattern features added to the same tree receive distinct
// OperationIds, and that those ids survive a save/load cycle.
//
// TODO: once the full replay pipeline is available, extend this to verify that
//   replay allocates instance_op_ids and that re-replay reuses the same ids.
// ---------------------------------------------------------------------------

#[test]
fn pattern_op_ids_distinct_and_stable_across_save_load() {
    let (mut tree, _sk, ex_id) = build_sketch_extrude(1.0);

    let count_param = add_param(&mut tree, "_count", 2.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_total", 10.0, Unit::Meters);

    let pat_a = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![ex_id],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "PatA",
    );

    let count2 = add_param(&mut tree, "_count2", 3.0, Unit::Ratio);
    let total2 = add_param(&mut tree, "_total2", 15.0, Unit::Meters);

    let pat_b = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![ex_id],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([0.0, 1.0, 0.0]),
                count: count2,
                total_distance: total2,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "PatB",
    );

    // The two patterns must have distinct OperationIds.
    let op_a = tree.feature(pat_a).expect("PatA must exist").operation_id;
    let op_b = tree.feature(pat_b).expect("PatB must exist").operation_id;
    assert_ne!(op_a, op_b, "PatA and PatB must have distinct OperationIds");

    // Round-trip: OperationIds survive save/load.
    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    // Collect the loaded patterns in timeline order.
    let loaded_patterns: Vec<_> = loaded
        .timeline()
        .iter()
        .filter(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::LinearPattern { .. })
            )
        })
        .collect();

    assert_eq!(
        loaded_patterns.len(),
        2,
        "both patterns must survive round-trip"
    );

    let loaded_op_a = loaded
        .feature(*loaded_patterns[0])
        .expect("first loaded pattern")
        .operation_id;
    let loaded_op_b = loaded
        .feature(*loaded_patterns[1])
        .expect("second loaded pattern")
        .operation_id;

    assert_ne!(
        loaded_op_a, loaded_op_b,
        "loaded patterns must still have distinct OperationIds"
    );

    // The loaded ops must match the originals (order preserved by timeline).
    assert_eq!(
        loaded_op_a, op_a,
        "PatA OperationId must survive round-trip"
    );
    assert_eq!(
        loaded_op_b, op_b,
        "PatB OperationId must survive round-trip"
    );
}

// ---------------------------------------------------------------------------
// Test 7: CircularPattern round-trips with axis and angle parameters
// ---------------------------------------------------------------------------

#[test]
fn circular_pattern_serializes_and_loads_correctly() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(1.0);

    let count_param = add_param(&mut tree, "_circ_count", 4.0, Unit::Ratio);
    let angle_param = add_param(
        &mut tree,
        "_circ_angle",
        std::f64::consts::TAU,
        Unit::Radians,
    );

    let pattern_id = tree.add_feature(
        Feature::CircularPattern {
            seed_features: vec![extrude_id],
            axis: AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            count: count_param,
            total_angle: angle_param,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "CircPattern1",
    );

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    let loaded_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::CircularPattern { .. })
            )
        })
        .expect("CircularPattern must survive round-trip");

    let loaded_node = loaded
        .feature(loaded_id)
        .expect("circular pattern node must exist");
    assert_eq!(
        loaded_node.label, "CircPattern1",
        "label must survive round-trip"
    );

    let Feature::CircularPattern {
        count, total_angle, ..
    } = &loaded_node.feature
    else {
        panic!("expected CircularPattern after load");
    };

    let count_val = loaded
        .parameters()
        .get(*count)
        .expect("count param must exist after load")
        .value;
    assert!(
        (count_val - 4.0).abs() < 1e-12,
        "count must be 4.0 after load, got {count_val}"
    );

    let angle_val = loaded
        .parameters()
        .get(*total_angle)
        .expect("total_angle param must exist after load")
        .value;
    assert!(
        (angle_val - std::f64::consts::TAU).abs() < 1e-12,
        "total_angle must be TAU after load, got {angle_val}"
    );

    // The loaded pattern must have exactly one dependency
    // (on the extrude feature's corresponding loaded node).
    assert_eq!(
        loaded_node.depends_on.len(),
        1,
        "circular pattern must have exactly one dependency"
    );

    // original pattern_id is dead after round-trip (new FeatureId assigned),
    // but we can verify the dependency points to the loaded extrude.
    let loaded_extrude_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::Extrude { .. })
            )
        })
        .expect("Extrude must survive round-trip");

    assert!(
        loaded_node.depends_on.contains(&loaded_extrude_id),
        "circular pattern must depend on the extrude feature"
    );

    // TODO: replay and assert 4 instances appear in current_brep once
    //   the full geometric pipeline is available.
    let _ = pattern_id; // suppress unused warning
}

// ---------------------------------------------------------------------------
// Test 8: Mirror feature round-trips correctly
// ---------------------------------------------------------------------------

#[test]
fn mirror_feature_serializes_and_loads_correctly() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(0.5);

    let pattern_id = tree.add_feature(
        Feature::Mirror {
            seed_features: vec![extrude_id],
            plane: pattern::MirrorPlaneDesc::XY,
            merge_result: true,
        },
        "Mirror1",
    );

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    let loaded_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::Mirror { .. })
            )
        })
        .expect("Mirror feature must survive round-trip");

    let loaded_node = loaded
        .feature(loaded_id)
        .expect("mirror node must exist after load");
    assert_eq!(
        loaded_node.label, "Mirror1",
        "label must survive round-trip"
    );

    let Feature::Mirror {
        merge_result,
        plane,
        ..
    } = &loaded_node.feature
    else {
        panic!("expected Mirror after load");
    };

    assert!(*merge_result, "merge_result must be true after load");
    assert!(
        matches!(plane, pattern::MirrorPlaneDesc::XY),
        "plane must be XY after load"
    );

    // TODO: replay and verify the mirrored geometry once the full BRep
    //   pipeline is available.
    let _ = pattern_id; // suppress unused warning
}

// ---------------------------------------------------------------------------
// Test 9: LinearPattern geometry_pattern flag persists
// ---------------------------------------------------------------------------

#[test]
fn linear_pattern_geometry_pattern_flag_persists() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(1.0);

    let count_param = add_param(&mut tree, "_gp_count", 2.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_gp_total", 5.0, Unit::Meters);

    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![extrude_id],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([0.0, 1.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: true, // explicitly enabled
            instance_overrides: Default::default(),
        },
        "GeomPattern",
    );

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    let loaded_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::LinearPattern { .. })
            )
        })
        .expect("LinearPattern must survive round-trip");

    let loaded_node = loaded.feature(loaded_id).expect("loaded node must exist");
    let Feature::LinearPattern {
        geometry_pattern, ..
    } = &loaded_node.feature
    else {
        panic!("expected LinearPattern after load");
    };

    assert!(
        *geometry_pattern,
        "geometry_pattern flag must be true after round-trip"
    );

    let _ = pattern_id; // suppress unused warning
}
