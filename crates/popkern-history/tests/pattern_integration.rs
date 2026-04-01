//! Integration tests for feature patterns in `popkern-history`.
//!
//! Tests cover tree-level operations: creation, serialization round-trips,
//! seed validation, skip instances, backward compatibility, and OperationId
//! stability.  Full geometric validation (BRep::validate) requires the
//! sketch+extrude pipeline; those cases are marked with TODO comments.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_history::*;

// ---------------------------------------------------------------------------
// Helpers
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
            sketch_data: rectangle_sketch_json(),
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
// Test 1: FeatureTree with LinearPattern serialises and loads correctly
// ---------------------------------------------------------------------------

#[test]
fn feature_tree_with_linear_pattern_serializes() {
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
        "LinPattern1",
    );

    // Tree should have 3 features: Sketch, Extrude, LinearPattern.
    assert_eq!(tree.timeline().len(), 3, "expected 3 features in timeline");

    // Pattern node must exist.
    let node = tree.feature(pattern_id).expect("pattern node must exist");
    assert_eq!(node.label, "LinPattern1");
    assert!(
        matches!(node.feature, Feature::LinearPattern { .. }),
        "feature variant must be LinearPattern"
    );

    // Save and load round-trip.
    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    assert_eq!(
        loaded.timeline().len(),
        3,
        "timeline length must survive round-trip"
    );

    // Find the pattern feature in the loaded tree.
    let loaded_pattern_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::LinearPattern { .. })
            )
        })
        .expect("LinearPattern must be present after load");

    let loaded_node = loaded
        .feature(loaded_pattern_id)
        .expect("loaded pattern node must exist");
    assert_eq!(
        loaded_node.label, "LinPattern1",
        "pattern label must survive round-trip"
    );

    // The count and total_distance parameters must still resolve.
    let Feature::LinearPattern { primary, .. } = &loaded_node.feature else {
        panic!("expected LinearPattern after load");
    };
    let count_val = loaded
        .parameters()
        .get(primary.count)
        .expect("count parameter must exist after load")
        .value;
    assert!(
        (count_val - 3.0).abs() < 1e-12,
        "count parameter value must be 3.0, got {count_val}"
    );

    let total_val = loaded
        .parameters()
        .get(primary.total_distance)
        .expect("total_distance parameter must exist after load")
        .value;
    assert!(
        (total_val - 20.0).abs() < 1e-12,
        "total_distance parameter value must be 20.0, got {total_val}"
    );

    // TODO: once the full sketch+extrude pipeline is wired, add:
    //   loaded.replay(Tolerance::default()).expect("replay must succeed");
    //   assert pattern instances appear in loaded.current_brep()
}

// ---------------------------------------------------------------------------
// Test 2: Pattern rejects empty seed features during tree-level validation
// ---------------------------------------------------------------------------

#[test]
fn pattern_rejects_empty_seeds_in_serialized_tree() {
    // We can build a LinearPattern with no seeds and add it to the tree
    // (add_feature does not call validate_seed_features).
    // The pattern can be stored and serialised; the error surfaces during replay.
    // This test verifies that the empty-seed pattern survives the serialisation
    // round-trip so that the replay error is deterministic.
    let mut tree = FeatureTree::new();

    let count_param = add_param(&mut tree, "_pat_count", 3.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_pat_total", 10.0, Unit::Meters);

    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![], // intentionally empty
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "EmptySeedPattern",
    );

    // The feature is present in the tree.
    let node = tree.feature(pattern_id).expect("pattern node must exist");
    let Feature::LinearPattern { seed_features, .. } = &node.feature else {
        panic!("expected LinearPattern");
    };
    assert!(
        seed_features.is_empty(),
        "seed_features must be empty as supplied"
    );

    // Round-trip: the empty seed list survives serialisation.
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

    let loaded_node = loaded.feature(loaded_id).expect("node must exist");
    let Feature::LinearPattern {
        seed_features: loaded_seeds,
        ..
    } = &loaded_node.feature
    else {
        panic!("expected LinearPattern after load");
    };
    assert!(
        loaded_seeds.is_empty(),
        "empty seed list must survive round-trip"
    );

    // TODO: once replay is wired end-to-end, assert that replaying this tree
    //   produces a FeatureResult::Failed for the pattern (EmptySeedFeatures).
}

// ---------------------------------------------------------------------------
// Test 3: Skip instance persists across save/load
// ---------------------------------------------------------------------------

#[test]
fn pattern_skip_instance_persists_across_round_trip() {
    let (mut tree, _sketch_id, extrude_id) = build_sketch_extrude(1.0);

    let count_param = add_param(&mut tree, "_pat_count", 4.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_pat_total", 30.0, Unit::Meters);

    // Build instance_overrides with instance index 3 suppressed.
    let mut overrides: std::collections::HashMap<usize, pattern::InstanceOverride> =
        Default::default();
    overrides.insert(
        3,
        pattern::InstanceOverride {
            suppressed: true,
            param_overrides: Default::default(),
        },
    );

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
            instance_overrides: overrides,
        },
        "PatternSkip",
    );

    // Verify instance 3 is suppressed before serialisation.
    {
        let node = tree.feature(pattern_id).expect("pattern node");
        let Feature::LinearPattern {
            instance_overrides, ..
        } = &node.feature
        else {
            panic!("expected LinearPattern");
        };
        let override_3 = instance_overrides
            .get(&3)
            .expect("instance 3 must have an override before save");
        assert!(
            override_3.suppressed,
            "instance 3 must be suppressed before save"
        );
    }

    // Round-trip.
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
        .expect("LinearPattern must exist after load");

    let loaded_node = loaded.feature(loaded_id).expect("loaded pattern node");
    let Feature::LinearPattern {
        instance_overrides: loaded_overrides,
        ..
    } = &loaded_node.feature
    else {
        panic!("expected LinearPattern after load");
    };

    let loaded_override_3 = loaded_overrides
        .get(&3)
        .expect("instance 3 override must survive round-trip");
    assert!(
        loaded_override_3.suppressed,
        "instance 3 must still be suppressed after load"
    );

    // Instances without overrides must not appear.
    assert!(
        loaded_overrides.get(&0).is_none() || !loaded_overrides.get(&0).unwrap().suppressed,
        "instance 0 must not be suppressed"
    );
    assert!(
        loaded_overrides.get(&1).is_none() || !loaded_overrides.get(&1).unwrap().suppressed,
        "instance 1 must not be suppressed"
    );
    assert!(
        loaded_overrides.get(&2).is_none() || !loaded_overrides.get(&2).unwrap().suppressed,
        "instance 2 must not be suppressed"
    );

    // TODO: replay and assert that only 3 instances appear in current_brep()
    //   once the full geometric pipeline is available.
}

// ---------------------------------------------------------------------------
// Test 4: v2 file loads with OperationMode::NewBody on Extrude features
// ---------------------------------------------------------------------------

#[test]
fn v2_extrude_loads_with_default_operation_mode() {
    // Construct a real v3 tree with an Extrude(NewBody), serialise it, downgrade
    // the version field to 2, strip the `mode` field, and re-load.  The v2→v3
    // migration must synthesise `mode: NewBody`.

    let (tree, _sk, _ex) = build_sketch_extrude(0.1);

    let bytes_v3 = tree.save().expect("save v3 should succeed");
    let mut doc: serde_json::Value =
        serde_json::from_slice(&bytes_v3).expect("v3 doc must be valid JSON");

    // Downgrade version.
    doc["version"] = serde_json::json!(2);

    // Walk feature slots and remove `mode` from every Extrude object.
    if let Some(slots) = doc["tree"]["features"].as_array_mut() {
        for slot in slots.iter_mut() {
            let extrude_obj = slot
                .as_object_mut()
                .and_then(|s| s.get_mut("value"))
                .and_then(serde_json::Value::as_object_mut)
                .and_then(|v| v.get_mut("feature"))
                .and_then(serde_json::Value::as_object_mut)
                .and_then(|f| f.get_mut("Extrude"))
                .and_then(serde_json::Value::as_object_mut);

            if let Some(obj) = extrude_obj {
                // v2 layout: flat `distance` field instead of `end_condition`,
                // and no `mode` field.
                let distance_param = obj
                    .get("end_condition")
                    .and_then(|ec| ec.get("Blind"))
                    .and_then(|b| b.get("distance"))
                    .cloned()
                    .expect("v3 must have end_condition.Blind.distance");

                obj.remove("end_condition");
                obj.remove("mode");
                obj.remove("draft_angle");
                obj.remove("second_direction");
                obj.insert("distance".to_string(), distance_param);
            }
        }
    }

    let bytes_v2 = serde_json::to_vec(&doc).expect("re-serialised v2 doc must be valid JSON");
    let loaded = FeatureTree::load(&bytes_v2).expect("v2→v3 migration must succeed");

    // Find the Extrude feature.
    let extrude_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::Extrude { .. })
            )
        })
        .expect("Extrude feature must survive v2 migration");

    let node = loaded.feature(extrude_id).expect("extrude node must exist");
    let Feature::Extrude { mode, .. } = &node.feature else {
        panic!("expected Extrude feature");
    };

    // The migration must default to NewBody.
    assert!(
        matches!(mode, OperationMode::NewBody),
        "v2 migration must default to OperationMode::NewBody, got {mode:?}"
    );
}

// ---------------------------------------------------------------------------
// Test 5: instance_op_ids default to empty and survive a save/load round-trip
//
// The `instance_op_ids` field uses `#[serde(skip_serializing_if = "is_empty")]`
// so a freshly created pattern (before any replay) serialises without that
// key.  On load the field defaults to an empty HashMap.
//
// TODO: once `FeatureTree::features` is accessible (or a public setter is
//   added), extend this test to verify that replay-populated op_ids survive
//   the round-trip.  Currently `HashMap<(u32,u32), OperationId>` cannot be
//   serialised to JSON via serde_json (tuple keys are unsupported); a custom
//   serde impl is required for non-empty maps.
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
