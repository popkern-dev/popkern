// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
// Tests for v2 → v3 migration (old Extrude/Revolve field layout).

use serde_json::Value;

use crate::FeatureTree;
use crate::feature::{AxisDesc, ExtrudeEnd, Feature, OperationMode, RevolveEnd, Vec3Desc};
use crate::unit::Unit;

use super::{add_param, face_name, rectangle_sketch_json};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Navigate a slot JSON value to the inner Extrude fields object if present.
///
/// Uses `as_object` / `get` (read-only) for navigation to avoid the
/// `serde_json::IndexMut` behaviour of inserting `null` entries for missing
/// keys when using the `value["key"]` syntax with a mutable receiver.
fn slot_has_extrude(slot: &Value) -> bool {
    slot.as_object()
        .and_then(|s| s.get("value"))
        .and_then(Value::as_object)
        .and_then(|v| v.get("feature"))
        .and_then(Value::as_object)
        .and_then(|f| f.get("Extrude"))
        .and_then(Value::as_object)
        .is_some()
}

/// Navigate a slot JSON value to the inner Revolve fields object if present.
fn slot_has_revolve(slot: &Value) -> bool {
    slot.as_object()
        .and_then(|s| s.get("value"))
        .and_then(Value::as_object)
        .and_then(|v| v.get("feature"))
        .and_then(Value::as_object)
        .and_then(|f| f.get("Revolve"))
        .and_then(Value::as_object)
        .is_some()
}

/// Build a synthetic v2 document with a valid `parameter_table` (taken from a
/// real v3 save) and an Extrude feature rewritten to the old v2 flat field
/// layout (`distance: ParameterId` instead of `end_condition: ExtrudeEnd`).
fn make_v2_extrude_bytes() -> Vec<u8> {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: crate::feature::SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let dist_param = add_param(&mut tree, "ext_distance", 0.08, Unit::Meters);
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

    // Capture a real (v3) saved document and downgrade it.
    let bytes = tree.save().expect("save should succeed");
    let mut doc: serde_json::Value =
        serde_json::from_slice(&bytes).expect("saved doc should be valid JSON");

    // Downgrade the version number so `load` dispatches to migrate_v2.
    doc["version"] = serde_json::json!(2);

    // Walk the feature slots and patch only Extrude slots.
    if let Some(slots) = doc["tree"]["features"].as_array_mut() {
        for slot in slots.iter_mut() {
            if !slot_has_extrude(slot) {
                continue;
            }
            // Safe mutable navigation using get_mut (does NOT insert nulls).
            let extrude_obj = slot
                .as_object_mut()
                .and_then(|s| s.get_mut("value"))
                .and_then(Value::as_object_mut)
                .and_then(|v| v.get_mut("feature"))
                .and_then(Value::as_object_mut)
                .and_then(|f| f.get_mut("Extrude"))
                .and_then(Value::as_object_mut)
                .expect("slot_has_extrude confirmed this path exists");

            // Extract the ParameterId from end_condition.Blind.distance.
            // Use index (read-only) since we only need a clone.
            let distance_param = extrude_obj
                .get("end_condition")
                .and_then(|ec| ec.get("Blind"))
                .and_then(|b| b.get("distance"))
                .cloned()
                .expect("v3 Extrude must have end_condition.Blind.distance");

            // Rewrite to v2 flat layout.
            extrude_obj.remove("end_condition");
            extrude_obj.remove("mode");
            extrude_obj.remove("draft_angle");
            extrude_obj.remove("second_direction");
            extrude_obj.insert("distance".to_string(), distance_param);
        }
    }

    serde_json::to_vec(&doc).expect("re-serialized doc should be valid JSON")
}

/// Build a synthetic v2 document with a Revolve feature rewritten to the old
/// v2 flat layout (`axis_origin`, `axis_direction`, `angle` instead of `axis`
/// and `end_condition`).
fn make_v2_revolve_bytes() -> Vec<u8> {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: crate::feature::SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let angle_param = add_param(&mut tree, "rev_angle", 3.14159, Unit::Radians);
    tree.add_feature(
        Feature::Revolve {
            face: face_name(1),
            axis: AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            end_condition: RevolveEnd::Blind { angle: angle_param },
            mode: OperationMode::NewBody,
        },
        "Revolve 1",
    );

    let bytes = tree.save().expect("save should succeed");
    let mut doc: serde_json::Value =
        serde_json::from_slice(&bytes).expect("saved doc should be valid JSON");

    doc["version"] = serde_json::json!(2);

    if let Some(slots) = doc["tree"]["features"].as_array_mut() {
        for slot in slots.iter_mut() {
            if !slot_has_revolve(slot) {
                continue;
            }
            let revolve_obj = slot
                .as_object_mut()
                .and_then(|s| s.get_mut("value"))
                .and_then(Value::as_object_mut)
                .and_then(|v| v.get_mut("feature"))
                .and_then(Value::as_object_mut)
                .and_then(|f| f.get_mut("Revolve"))
                .and_then(Value::as_object_mut)
                .expect("slot_has_revolve confirmed this path exists");

            // Clone values before mutation.
            let origin = revolve_obj
                .get("axis")
                .and_then(|a| a.get("Literal"))
                .and_then(|l| l.get("origin"))
                .cloned()
                .expect("v3 Revolve must have axis.Literal.origin");
            let direction = revolve_obj
                .get("axis")
                .and_then(|a| a.get("Literal"))
                .and_then(|l| l.get("direction"))
                .cloned()
                .expect("v3 Revolve must have axis.Literal.direction");
            let angle = revolve_obj
                .get("end_condition")
                .and_then(|ec| ec.get("Blind"))
                .and_then(|b| b.get("angle"))
                .cloned()
                .expect("v3 Revolve must have end_condition.Blind.angle");

            // Rewrite to v2 flat layout.
            revolve_obj.remove("axis");
            revolve_obj.remove("end_condition");
            revolve_obj.remove("mode");
            revolve_obj.insert("axis_origin".to_string(), origin);
            revolve_obj.insert("axis_direction".to_string(), direction);
            revolve_obj.insert("angle".to_string(), angle);
        }
    }

    serde_json::to_vec(&doc).expect("re-serialized doc should be valid JSON")
}

// ---------------------------------------------------------------------------
// Test 13: Load v2 Extrude → migrates to ExtrudeEnd::Blind
// ---------------------------------------------------------------------------

#[test]
fn v2_extrude_roundtrip() {
    let bytes = make_v2_extrude_bytes();
    let loaded = FeatureTree::load(&bytes).expect("v2 → v3 migration should succeed");

    // Two features: Sketch + Extrude.
    assert_eq!(
        loaded.timeline().len(),
        2,
        "should have 2 features after v2 migration"
    );

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
        .expect("Extrude feature must be present");

    let node = loaded.feature(extrude_id).expect("feature must exist");

    let Feature::Extrude {
        end_condition,
        direction,
        ..
    } = &node.feature
    else {
        panic!("expected Extrude, got something else");
    };

    // direction should have been preserved (still Some after migration).
    assert!(
        direction.is_some(),
        "direction should be Some after migration"
    );

    // end_condition must be Blind and the distance parameter must resolve.
    let ExtrudeEnd::Blind { distance } = end_condition else {
        panic!("expected ExtrudeEnd::Blind, got {end_condition:?}");
    };

    let param = loaded
        .parameters()
        .get(*distance)
        .expect("distance ParameterId must still be valid after migration");
    assert!(
        (param.value - 0.08).abs() < 1e-12,
        "distance parameter value must be 0.08, got {}",
        param.value
    );
}

// ---------------------------------------------------------------------------
// Test 14: Load v2 Revolve → migrates to AxisDesc::Literal + RevolveEnd::Blind
// ---------------------------------------------------------------------------

#[test]
fn v2_revolve_roundtrip() {
    let bytes = make_v2_revolve_bytes();
    let loaded = FeatureTree::load(&bytes).expect("v2 → v3 revolve migration should succeed");

    assert_eq!(
        loaded.timeline().len(),
        2,
        "should have 2 features after v2 revolve migration"
    );

    let revolve_id = *loaded
        .timeline()
        .iter()
        .find(|&&id| {
            matches!(
                loaded.feature(id).map(|n| &n.feature),
                Some(Feature::Revolve { .. })
            )
        })
        .expect("Revolve feature must be present");

    let node = loaded.feature(revolve_id).expect("feature must exist");

    let Feature::Revolve {
        axis,
        end_condition,
        ..
    } = &node.feature
    else {
        panic!("expected Revolve feature");
    };

    // axis must be AxisDesc::Literal with the original origin/direction.
    let AxisDesc::Literal { origin, direction } = axis else {
        panic!("expected AxisDesc::Literal, got {axis:?}");
    };
    assert_eq!(*origin, [0.0, 0.0, 0.0], "axis origin must be preserved");
    assert_eq!(
        *direction,
        [0.0, 0.0, 1.0],
        "axis direction must be preserved"
    );

    // end_condition must be RevolveEnd::Blind with the original angle parameter.
    let RevolveEnd::Blind { angle } = end_condition else {
        panic!("expected RevolveEnd::Blind, got {end_condition:?}");
    };
    let param = loaded
        .parameters()
        .get(*angle)
        .expect("angle ParameterId must still be valid after migration");
    assert!(
        (param.value - 3.14159).abs() < 1e-6,
        "angle parameter value must be ~π, got {}",
        param.value
    );
}
