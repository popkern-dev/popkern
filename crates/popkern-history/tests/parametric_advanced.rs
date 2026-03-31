//! Integration tests for the parametric expression pipeline in popkern-history
//! (tests 6–10): save/load, undo/redo, delete, bounds clamping, built-ins.
//!
//! See `parametric.rs` for tests 1–5.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

// ---------------------------------------------------------------------------
// Helpers (duplicated from parametric.rs — integration test files cannot
// share code via `mod` without a wrapper crate)
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

fn feature_op_id(tree: &FeatureTree, id: FeatureId) -> u64 {
    tree.feature(id).expect("feature").operation_id.0
}

// ---------------------------------------------------------------------------
// Test 6: Save/load with expressions round-trip
// ---------------------------------------------------------------------------

#[test]
fn save_load_with_expressions_round_trip() {
    let mut tree = FeatureTree::new();

    // base = 0.02 SI (20mm)
    let base = tree
        .add_global_parameter("base", 0.02, Unit::Millimeters)
        .expect("add base");

    // Sketch gets op_id = 1
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let sketch_op = feature_op_id(&tree, sketch_id);

    let dist = tree
        .add_global_parameter("_sl_dist", 0.0, Unit::Meters)
        .expect("add dist");

    // dist = base * 3  →  0.02 * 3 = 0.06m
    let expr = Expression::BinOp(
        Box::new(Expression::Param(base)),
        BinOp::Mul,
        Box::new(Expression::Literal(3.0)),
    );
    tree.set_parameter_expression(dist, expr)
        .expect("set expression");

    tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(sketch_op),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude1",
    );

    tree.replay(Tolerance::default()).expect("initial replay");

    // Verify expression before save: dist should be 0.06m
    let dist_before = tree.parameters().get(dist).expect("dist").value;
    assert!(
        (dist_before - 0.06).abs() < 1e-10,
        "dist should be 0.06m before save, got {dist_before}"
    );

    // Save and load
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");

    assert!(
        loaded.current_brep().is_none(),
        "loaded tree should have no brep before replay"
    );

    // Parameters should have survived
    let loaded_base = loaded
        .parameters()
        .resolve_name("base")
        .expect("base param should exist after load");
    assert!(
        loaded.parameters().get(loaded_base).is_some(),
        "base parameter should survive serialization"
    );

    // Replay loaded tree
    loaded.replay(Tolerance::default()).expect("loaded replay");

    let brep = loaded.current_brep().expect("loaded brep");
    assert_eq!(
        brep.faces.len(),
        6,
        "loaded + replayed tree should produce box"
    );

    // Verify expression drove the correct height
    let max_z: f64 = brep
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (max_z - 0.06).abs() < 0.01,
        "loaded tree height should be 0.06m, got {max_z}"
    );
}

// ---------------------------------------------------------------------------
// Test 7: Undo parameter operations
// ---------------------------------------------------------------------------

#[test]
fn undo_parameter_expression_and_replay() {
    let mut tree = FeatureTree::new();

    // base_undo = 0.03 SI (30mm)
    let base = tree
        .add_global_parameter("base_undo", 0.03, Unit::Millimeters)
        .expect("add base");

    // Sketch gets op_id = 1
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let sketch_op = feature_op_id(&tree, sketch_id);

    // dist initially 0.05m
    let dist = tree
        .add_global_parameter("_undo_dist", 0.05, Unit::Meters)
        .expect("add dist");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(sketch_op),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude1",
    );

    tree.replay(Tolerance::default()).expect("initial replay");
    assert!(
        matches!(
            tree.feature(extrude_id).unwrap().last_result.as_ref(),
            Some(FeatureResult::Success { .. })
        ),
        "initial extrude should succeed"
    );

    // Baseline: no expression yet
    assert!(
        tree.parameters()
            .get(dist)
            .expect("dist")
            .expression
            .is_none(),
        "dist should have no expression initially"
    );

    // Set expression: dist = base_undo (0.03m)
    let expr = Expression::Param(base);
    tree.set_parameter_expression(dist, expr)
        .expect("set expression");

    assert!(
        tree.parameters()
            .get(dist)
            .expect("dist")
            .expression
            .is_some(),
        "expression should be set"
    );

    // Replay with expression active — extrude distance becomes 0.03m
    tree.replay(Tolerance::default())
        .expect("replay after set expr");
    assert!(
        matches!(
            tree.feature(extrude_id).unwrap().last_result.as_ref(),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed with expression"
    );

    let after_expr_z: f64 = tree
        .current_brep()
        .expect("brep")
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (after_expr_z - 0.03).abs() < 0.01,
        "height after expression should be 0.03m, got {after_expr_z}"
    );

    // Undo the set_expression — expression is cleared.
    // Note: the stored value (0.03m) is NOT automatically reverted by undo;
    // only the expression field is cleared. This is by design.
    tree.undo();
    assert!(
        tree.parameters()
            .get(dist)
            .expect("dist")
            .expression
            .is_none(),
        "expression should be cleared after undo"
    );

    // Replay without expression — uses the stored value (still 0.03m)
    tree.replay(Tolerance::default())
        .expect("replay after undo");
    assert!(
        matches!(
            tree.feature(extrude_id).unwrap().last_result.as_ref(),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed after undo"
    );

    // Redo — expression comes back and drives the parameter again
    tree.redo();
    assert!(
        tree.parameters()
            .get(dist)
            .expect("dist")
            .expression
            .is_some(),
        "expression should be restored after redo"
    );

    tree.replay(Tolerance::default())
        .expect("replay after redo");
    assert!(
        matches!(
            tree.feature(extrude_id).unwrap().last_result.as_ref(),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed after redo"
    );

    // After redo + replay: expression is active, value driven by base_undo = 0.03m
    let after_redo_z: f64 = tree
        .current_brep()
        .expect("brep")
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (after_redo_z - 0.03).abs() < 0.01,
        "height should be 0.03m after redo, got {after_redo_z}"
    );
}

// ---------------------------------------------------------------------------
// Test 8: Delete cleans up parameters
// ---------------------------------------------------------------------------

#[test]
fn delete_cleans_up_feature_parameters() {
    let mut tree = FeatureTree::new();

    // Sketch gets op_id = 1
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let sketch_op = feature_op_id(&tree, sketch_id);

    let dist = tree
        .add_global_parameter("_del_dist", 0.1, Unit::Meters)
        .expect("add dist");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(sketch_op),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude1",
    );

    tree.replay(Tolerance::default()).expect("replay");
    assert_eq!(tree.current_brep().expect("brep").faces.len(), 6);

    // Delete extrude cascade
    let del = tree
        .delete_feature(extrude_id, DeleteBehavior::Cascade)
        .expect("delete");
    assert!(
        del.deleted.contains(&extrude_id),
        "extrude should be in deleted list"
    );

    // Feature-scoped params for extrude should be removed
    let feature_params = tree.parameters().feature_params(extrude_id);
    assert!(
        feature_params.is_empty(),
        "feature-scoped params should be removed after cascade delete"
    );

    // Replay should succeed with just the sketch
    tree.replay(Tolerance::default())
        .expect("replay after delete");
    let face_count = tree.current_brep().expect("brep").faces.len();
    assert!(
        face_count < 6,
        "after deleting extrude, should have fewer than 6 faces"
    );
}

// ---------------------------------------------------------------------------
// Test 9: Parameter bounds clamping
// ---------------------------------------------------------------------------

#[test]
fn parameter_bounds_clamping() {
    let mut tree = FeatureTree::new();

    // Create param with bounds [0.001, 0.1] in SI
    let radius = tree
        .add_global_parameter("_clamp_radius", 0.05, Unit::Meters)
        .expect("add radius");
    tree.set_parameter_bounds(radius, 0.001, 0.1);

    // Expression evaluates to 0.2 (above max = 0.1)
    tree.set_parameter_expression(radius, Expression::Literal(0.2))
        .expect("set expression");

    // Need at least one feature so replay runs
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    tree.replay(Tolerance::default()).expect("replay");

    let val = tree.parameters().get(radius).expect("radius").value;
    assert!(
        (val - 0.1).abs() < 1e-10,
        "value should be clamped to max 0.1, got {val}"
    );
}

// ---------------------------------------------------------------------------
// Test 10: Built-in constants accessible via parser
// ---------------------------------------------------------------------------

#[test]
fn builtin_constants_accessible() {
    let tree = FeatureTree::new();

    // Parse "pi / 4"
    let expr =
        parse_expression("pi / 4", tree.parameters()).expect("parse 'pi / 4' should succeed");

    let result = {
        use popkern_history::MeasureContext;
        use popkern_history::expression::evaluate;
        evaluate(&expr, tree.parameters(), &MeasureContext::empty()).expect("evaluate")
    };

    let expected = std::f64::consts::PI / 4.0;
    assert!(
        (result - expected).abs() < 1e-10,
        "pi / 4 should be ≈ {expected:.7}, got {result}"
    );

    // Verify approximate value
    assert!(
        (result - 0.785_398).abs() < 1e-4,
        "pi/4 ≈ 0.785398, got {result}"
    );
}
