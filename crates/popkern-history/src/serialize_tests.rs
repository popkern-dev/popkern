// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

use crate::FeatureTree;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::parameter::ParameterId;
use crate::unit::Unit;

/// Create a solved rectangle sketch and serialize it to JSON.
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

/// Build an `EntityName` referencing the generated face at `op_id`.
fn face_name(op_id: u64) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    }
}

fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

// ---------------------------------------------------------------------------
// Test 1: Build tree → save → load → same timeline and features
// ---------------------------------------------------------------------------

#[test]
fn save_load_preserves_structure() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let dist_param = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind {
                distance: dist_param,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude 1",
    );

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    // Timeline length.
    assert_eq!(
        loaded.timeline().len(),
        tree.timeline().len(),
        "timeline length should match"
    );

    // next_op_id.
    assert_eq!(
        loaded.next_op_id, tree.next_op_id,
        "next_op_id should match"
    );

    // Feature labels.
    let original_labels: Vec<&str> = tree
        .timeline()
        .iter()
        .map(|&id| tree.feature(id).expect("feature exists").label.as_str())
        .collect();
    let loaded_labels: Vec<&str> = loaded
        .timeline()
        .iter()
        .map(|&id| loaded.feature(id).expect("feature exists").label.as_str())
        .collect();
    assert_eq!(
        loaded_labels, original_labels,
        "feature labels should match"
    );
}

// ---------------------------------------------------------------------------
// Test 2: Session-local fields are NOT serialized
// ---------------------------------------------------------------------------

#[test]
fn session_local_fields_not_serialized() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Replay to populate current_brep and cache.
    tree.replay(Tolerance::default())
        .expect("replay should succeed");
    assert!(
        tree.current_brep().is_some(),
        "current_brep should be Some after replay"
    );

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    // current_brep must be None (session-local, not serialized).
    assert!(
        loaded.current_brep().is_none(),
        "loaded tree should have no current_brep"
    );

    // cache must be empty (session-local, not serialized).
    assert_eq!(
        loaded.cache.memory_bytes(),
        0,
        "loaded tree cache should be empty"
    );

    // undo_log must be empty (session-local, not serialized).
    assert_eq!(
        loaded.undo_log.len(),
        0,
        "loaded tree undo_log should be empty"
    );

    // undo_cursor must be 0.
    assert_eq!(loaded.undo_cursor, 0, "loaded tree undo_cursor should be 0");
}

// ---------------------------------------------------------------------------
// Test 3: Save → load → replay produces the same BRep
// ---------------------------------------------------------------------------

#[test]
fn save_load_replay_produces_same_brep() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let dist_param = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind {
                distance: dist_param,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude 1",
    );

    tree.replay(Tolerance::default())
        .expect("original replay should succeed");
    let original_face_count = tree.current_brep().expect("original brep").faces.len();

    let bytes = tree.save().expect("save should succeed");
    let mut loaded = FeatureTree::load(&bytes).expect("load should succeed");
    loaded
        .replay(Tolerance::default())
        .expect("loaded replay should succeed");

    let loaded_face_count = loaded.current_brep().expect("loaded brep").faces.len();

    assert_eq!(
        loaded_face_count, original_face_count,
        "face count should match after save → load → replay"
    );
}

// ---------------------------------------------------------------------------
// Test 4: Empty tree round-trip
// ---------------------------------------------------------------------------

#[test]
fn empty_tree_round_trip() {
    let tree = FeatureTree::new();

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    assert_eq!(
        loaded.timeline().len(),
        0,
        "empty tree should load with empty timeline"
    );
    assert_eq!(
        loaded.next_op_id, 1,
        "empty tree should start next_op_id at 1"
    );
    assert!(loaded.rollback_index.is_none(), "no rollback bar");
}

// ---------------------------------------------------------------------------
// Test 5: Load from invalid bytes returns error
// ---------------------------------------------------------------------------

#[test]
fn load_invalid_bytes_returns_error() {
    let result = FeatureTree::load(b"not valid json at all !!!!");
    assert!(result.is_err(), "loading invalid bytes should fail");

    let err = result.unwrap_err();
    let msg = err.to_string();
    assert!(
        msg.contains("serialization error"),
        "error message should mention serialization: {msg}"
    );
}

// ---------------------------------------------------------------------------
// Test 6: Rollback index is preserved across round-trip
// ---------------------------------------------------------------------------

#[test]
fn rollback_index_preserved() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let dist_param = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind {
                distance: dist_param,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude 1",
    );

    // Roll back to just before the extrude.
    tree.rollback_to(1);
    assert_eq!(tree.rollback_position(), Some(1));

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    assert_eq!(
        loaded.rollback_position(),
        Some(1),
        "rollback index should survive round-trip"
    );
}

// ---------------------------------------------------------------------------
// Test 7: Round-trip for all feature types
// ---------------------------------------------------------------------------

#[test]
fn save_load_all_feature_types() {
    use crate::feature::{AxisDesc, BoolOpKind, RevolveEnd, ShellDirectionKind, Vec3Desc};
    use crate::node::BodyId;

    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: "base".into(),
        },
        "Sketch",
    );

    let dist = add_param(&mut tree, "_d", 10.0, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );

    let r = add_param(&mut tree, "_r", 2.0, Unit::Millimeters);
    tree.add_feature(
        Feature::Fillet {
            edges: vec![face_name(2)],
            radius: r,
        },
        "Fillet",
    );

    let cd = add_param(&mut tree, "_cd", 1.5, Unit::Millimeters);
    tree.add_feature(
        Feature::Chamfer {
            edges: vec![face_name(2)],
            distance: cd,
        },
        "Chamfer",
    );

    let st = add_param(&mut tree, "_st", 0.5, Unit::Millimeters);
    tree.add_feature(
        Feature::Shell {
            faces_to_remove: vec![face_name(2)],
            thickness: st,
            direction: ShellDirectionKind::Outward,
        },
        "Shell",
    );

    let da = add_param(&mut tree, "_da", 0.1, Unit::Radians);
    tree.add_feature(
        Feature::Draft {
            faces: vec![face_name(2)],
            angle: da,
            neutral_plane: face_name(2),
        },
        "Draft",
    );

    tree.add_feature(
        Feature::Boolean {
            operation: BoolOpKind::Subtraction,
            tool_body: BodyId(1),
        },
        "Boolean",
    );

    let ra = add_param(&mut tree, "_ra", 3.14, Unit::Radians);
    tree.add_feature(
        Feature::Revolve {
            face: face_name(2),
            axis: AxisDesc::Literal {
                origin: [0.0; 3],
                direction: [0.0, 1.0, 0.0],
            },
            end_condition: RevolveEnd::Blind { angle: ra },
            mode: OperationMode::Boss,
        },
        "Revolve",
    );

    let tw = add_param(&mut tree, "_tw", 0.0, Unit::Radians);
    tree.add_feature(
        Feature::Sweep {
            face: face_name(2),
            path_edge: face_name(2),
            twist: tw,
            mode: OperationMode::Cut,
        },
        "Sweep",
    );

    tree.add_feature(
        Feature::Loft {
            sections: vec![face_name(2), face_name(3)],
            closed: true,
            mode: OperationMode::NewBody,
        },
        "Loft",
    );

    let pr = add_param(&mut tree, "_pr", 3.0, Unit::Millimeters);
    tree.add_feature(
        Feature::Pipe {
            path_edge: face_name(2),
            radius: pr,
            mode: OperationMode::NewBody,
        },
        "Pipe",
    );

    let bytes = tree.save().expect("save all features");
    let loaded = FeatureTree::load(&bytes).expect("load all features");

    assert_eq!(loaded.timeline().len(), tree.timeline().len());
    let labels = [
        "Sketch", "Extrude", "Fillet", "Chamfer", "Shell", "Draft", "Boolean", "Revolve", "Sweep",
        "Loft", "Pipe",
    ];
    for (i, &fid) in loaded.timeline().iter().enumerate() {
        assert_eq!(loaded.feature(fid).unwrap().label, labels[i]);
    }
    let ids = loaded.timeline().to_vec();
    let f = |i: usize| &loaded.feature(ids[i]).unwrap().feature;
    assert!(matches!(f(2), Feature::Fillet { .. }));
    assert!(matches!(f(4), Feature::Shell { .. }));
    assert!(matches!(f(6), Feature::Boolean { .. }));
    assert!(matches!(f(7), Feature::Revolve { .. }));
    assert!(matches!(f(9), Feature::Loft { .. }));
    assert!(matches!(f(10), Feature::Pipe { .. }));
}

// ---------------------------------------------------------------------------
// Test 8: Suppressed + config state round-trip
// ---------------------------------------------------------------------------

#[test]
fn save_load_preserves_suppressed_and_config() {
    let mut tree = FeatureTree::new();
    let s0 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        "S0",
    );
    let _s1 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: String::new(),
        },
        "S1",
    );
    tree.suppress(s0);

    let cfg_id = tree.add_config("Big", "big").expect("add config");
    let pid = tree
        .add_global_parameter("width", 5.0, Unit::Meters)
        .expect("param");
    tree.set_config_override(cfg_id, pid, crate::config::ParamOverride::Literal(20.0))
        .expect("override");
    tree.set_active_config(Some(cfg_id));

    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");

    let ids = loaded.timeline().to_vec();
    assert!(loaded.feature(ids[0]).unwrap().suppressed);
    assert!(!loaded.feature(ids[1]).unwrap().suppressed);
    assert_eq!(loaded.configs().len(), 1);
    let cfg = loaded.config(loaded.configs()[0]).unwrap();
    assert_eq!(cfg.name, "Big");
    assert_eq!(
        loaded.config_overrides(loaded.configs()[0]).unwrap().len(),
        1
    );
    assert!(loaded.active_config().is_some());
}
