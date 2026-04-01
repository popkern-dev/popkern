// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
// Tests 7–12 for the v2 serialization format and v1→v2 migration.

use popkern_core::tolerance::Tolerance;

use crate::FeatureTree;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc};
use crate::unit::Unit;

use super::{add_param, face_name, rectangle_sketch_json};

// ---------------------------------------------------------------------------
// Test 7: version number in saved bytes is 3
// ---------------------------------------------------------------------------

#[test]
fn saved_document_version_is_3() {
    let tree = FeatureTree::new();
    let bytes = tree.save().expect("save should succeed");
    let doc: serde_json::Value = serde_json::from_slice(&bytes).expect("should be valid JSON");
    let version = doc["version"].as_u64().expect("version must be u64");
    assert_eq!(version, 3, "saved document should have version 3");
}

// ---------------------------------------------------------------------------
// Test 8: Save v2 → load → ParameterIds still resolve
// ---------------------------------------------------------------------------

#[test]
fn parameters_preserved_across_round_trip() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let dist_param = add_param(&mut tree, "ext_distance", 2.5, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(crate::feature::Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind {
                distance: dist_param,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude 1",
    );
    // After add_feature, `ext_distance` is re-indexed as "extrude_1.distance"
    // because register_feature_params moves feature-scoped params to qualified names.
    // Confirm the qualified name exists pre-save so the test expectation is clear.
    assert!(
        tree.parameters()
            .resolve_name("extrude_1.distance")
            .is_some(),
        "parameter should be indexed as extrude_1.distance before save"
    );

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    // The parameter should resolve by its qualified feature-scoped name.
    let resolved = loaded.parameters().resolve_name("extrude_1.distance");
    assert!(
        resolved.is_some(),
        "extrude_1.distance must resolve after load"
    );

    let param = loaded
        .parameters()
        .get(resolved.expect("already checked"))
        .expect("parameter must exist");
    assert!(
        (param.value - 2.5).abs() < 1e-12,
        "parameter value must be preserved: got {}",
        param.value
    );
}

// ---------------------------------------------------------------------------
// Test 9: Save v2 → load → expression preserved
// ---------------------------------------------------------------------------

#[test]
fn expression_preserved_across_round_trip() {
    use crate::expression::{BinOp, Expression};

    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Add two parameters: base and derived (= base * 2).
    let base_id = add_param(&mut tree, "base_len", 1.0, Unit::Meters);
    let derived_id = add_param(&mut tree, "derived_len", 2.0, Unit::Meters);
    tree.parameters_mut()
        .set_expression(
            derived_id,
            Expression::BinOp(
                Box::new(Expression::Param(base_id)),
                BinOp::Mul,
                Box::new(Expression::Literal(2.0)),
            ),
        )
        .expect("set expression should succeed");

    let bytes = tree.save().expect("save should succeed");
    let loaded = FeatureTree::load(&bytes).expect("load should succeed");

    let derived_loaded = loaded
        .parameters()
        .resolve_name("derived_len")
        .expect("derived_len must resolve");
    let param = loaded
        .parameters()
        .get(derived_loaded)
        .expect("parameter must exist");

    assert!(
        param.expression.is_some(),
        "expression must survive round-trip"
    );
}

// ---------------------------------------------------------------------------
// Test 10: Load version 99 data returns UnsupportedVersion
// ---------------------------------------------------------------------------

#[test]
fn load_unsupported_version_returns_error() {
    use crate::error::HistoryError;

    let future_doc = serde_json::json!({
        "version": 99,
        "tree": {}
    });
    let bytes = serde_json::to_vec(&future_doc).expect("serialize");

    let err = FeatureTree::load(&bytes).expect_err("should fail with unsupported version");
    assert!(
        matches!(
            err,
            HistoryError::UnsupportedVersion {
                file_version: 99,
                max_supported: 3,
            }
        ),
        "expected UnsupportedVersion {{99, 3}}, got: {err:?}"
    );
}

// ---------------------------------------------------------------------------
// Test 11: Load v1 JSON → parameters auto-created from f64 values
// ---------------------------------------------------------------------------

#[test]
fn v1_migration_creates_parameters_from_f64() {
    // Construct a minimal v1 JSON document manually.
    // V1 format: features is a SlotMap array of {version, value} slots.
    // Slot 0 is the vacant sentinel; slot 1 is the Sketch; slot 2 is Extrude.
    // Timeline references slot idx 1 and 2.
    // Extrude had distance as a bare f64 in v1.
    let v1_doc = serde_json::json!({
        "version": 1,
        "tree": {
            "features": [
                {"version": 0, "value": null},
                {
                    "version": 1,
                    "value": {
                        "body": 0,
                        "depends_on": [],
                        "feature": {"Sketch": {"plane": "XY", "sketch_data": "{}"}},
                        "label": "Sketch 1",
                        "operation_id": 1,
                        "suppressed": false,
                        "suppress_if": null
                    }
                },
                {
                    "version": 1,
                    "value": {
                        "body": 0,
                        "depends_on": [{"idx": 1, "version": 1}],
                        "feature": {
                            "Extrude": {
                                "face": {"ancestors": [], "operation": 1,
                                         "role": {"category": "Generated", "index": 0}},
                                "direction": {"Literal": [0.0, 0.0, 1.0]},
                                "distance": 0.05
                            }
                        },
                        "label": "Extrude 1",
                        "operation_id": 2,
                        "suppressed": false,
                        "suppress_if": null
                    }
                }
            ],
            "timeline": [
                {"idx": 1, "version": 1},
                {"idx": 2, "version": 1}
            ],
            "rollback_index": null,
            "next_op_id": 3,
            "frozen_index": null
        }
    });

    let bytes = serde_json::to_vec(&v1_doc).expect("serialize v1 fixture");
    let loaded = FeatureTree::load(&bytes).expect("v1 migration should succeed");

    // Timeline should have 2 features.
    assert_eq!(
        loaded.timeline().len(),
        2,
        "v1 migration should produce 2 features"
    );

    // The Extrude feature's distance parameter must have been auto-created.
    let all_params: Vec<_> = loaded.parameters().all_params().collect();
    let user_params: Vec<_> = all_params.iter().filter(|(_, p)| !p.read_only).collect();
    assert!(
        !user_params.is_empty(),
        "v1 migration should have auto-created at least one parameter for the extrude distance"
    );

    // The auto-created parameter should have value 0.05 (the original f64).
    let has_correct_value = user_params
        .iter()
        .any(|(_, p)| (p.value - 0.05).abs() < 1e-12);
    assert!(
        has_correct_value,
        "migrated parameter should have value 0.05 (the original f64 distance)"
    );
}

// ---------------------------------------------------------------------------
// Test 12: Round-trip with expressions → replay → same face count
// ---------------------------------------------------------------------------

#[test]
fn round_trip_with_expressions_produces_same_brep() {
    use crate::expression::{BinOp, Expression};

    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let base_id = add_param(&mut tree, "base_len", 1.0, Unit::Meters);
    let dist_id = add_param(&mut tree, "dist_len", 1.0, Unit::Meters);
    // dist_len = base_len * 1 (simple expression keeping value at 1.0 m)
    tree.parameters_mut()
        .set_expression(
            dist_id,
            Expression::BinOp(
                Box::new(Expression::Param(base_id)),
                BinOp::Mul,
                Box::new(Expression::Literal(1.0)),
            ),
        )
        .expect("set expression should succeed");

    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(crate::feature::Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_id },
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
        "face count must match after round-trip with expression-driven parameter"
    );
}
