// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//
// Tests for two-phase parameter evaluation, suppress_if, and transaction
// rollback — the replay integration features added in Task 5.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

use crate::expression::Expression;
use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};
use crate::unit::Unit;
use crate::{FeatureResult, FeatureTree};

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

/// Create an entity name referencing a given OperationId.
fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Test 1: suppress_if = Literal(1.0) → feature skipped
// ---------------------------------------------------------------------------

#[test]
fn suppress_if_true_skips_feature() {
    let mut tree = FeatureTree::new();

    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch 1",
    );

    let dist = tree
        .parameters_mut()
        .add_global("dist", 1.0, Unit::Meters)
        .expect("add param");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude 1",
    );

    // Conditional suppression: Literal(1.0) is always truthy.
    tree.features[extrude_id].suppress_if = Some(Expression::Literal(1.0));

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Suppressed)
        ),
        "extrude should be auto-suppressed by suppress_if, got {:?}",
        result.feature_results.get(&extrude_id)
    );
    assert!(
        matches!(
            result.feature_results.get(&sketch_id),
            Some(FeatureResult::Success { .. })
        ),
        "sketch should still succeed"
    );

    // Flat sketch result — no extrusion.
    let brep = tree.current_brep().expect("brep");
    assert!(
        brep.faces.len() < 6,
        "auto-suppressed extrude means fewer than 6 faces, got {}",
        brep.faces.len()
    );
}

// ---------------------------------------------------------------------------
// Test 2: suppress_if = Literal(0.0) → feature executes normally
// ---------------------------------------------------------------------------

#[test]
fn suppress_if_false_executes_feature() {
    let mut tree = FeatureTree::new();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch 1",
    );

    let dist = tree
        .parameters_mut()
        .add_global("dist2", 1.0, Unit::Meters)
        .expect("add param");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude 1",
    );

    // Conditional suppression: Literal(0.0) is always falsy → feature runs.
    tree.features[extrude_id].suppress_if = Some(Expression::Literal(0.0));

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed when suppress_if is false, got {:?}",
        result.feature_results.get(&extrude_id)
    );

    let brep = tree.current_brep().expect("brep");
    assert_eq!(brep.faces.len(), 6, "false suppress_if → 6-face box");
}

// ---------------------------------------------------------------------------
// Test 3: Manual suppress=true overrides suppress_if=false
// ---------------------------------------------------------------------------

#[test]
fn manual_suppress_independent_of_suppress_if() {
    let mut tree = FeatureTree::new();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch 1",
    );

    let dist = tree
        .parameters_mut()
        .add_global("dist3", 1.0, Unit::Meters)
        .expect("add param");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Extrude 1",
    );

    // Manual suppress = true, suppress_if = false — manual flag still wins.
    tree.features[extrude_id].suppressed = true;
    tree.features[extrude_id].suppress_if = Some(Expression::Literal(0.0));

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Suppressed)
        ),
        "manual suppressed=true should suppress even when suppress_if is false, got {:?}",
        result.feature_results.get(&extrude_id)
    );
}

// ---------------------------------------------------------------------------
// Test 4: Expression evaluation propagates parameter before replay
// ---------------------------------------------------------------------------

#[test]
fn expression_evaluation_propagates_before_replay() {
    let mut tree = FeatureTree::new();

    // Global "height" param starts at 0.1 m.
    let height_id = tree
        .parameters_mut()
        .add_global("height", 0.1, Unit::Meters)
        .expect("add height");

    // "derived_dist" is driven by expression Param(height_id).
    let derived_id = tree
        .parameters_mut()
        .add_global("derived_dist", 0.0, Unit::Meters)
        .expect("add derived_dist");
    tree.parameters_mut()
        .get_mut(derived_id)
        .expect("derived_dist")
        .expression = Some(Expression::Param(height_id));

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch 1",
    );

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: derived_id,
        },
        "Extrude 1",
    );

    // Change height to 0.2; derived_dist should follow via expression.
    tree.parameters_mut()
        .get_mut(height_id)
        .expect("height")
        .value = 0.2;

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Success { .. })
        ),
        "extrude should succeed with propagated height, got {:?}",
        result.feature_results.get(&extrude_id)
    );

    // The derived parameter's value should have been updated to 0.2.
    let derived_val = tree
        .parameters()
        .get(derived_id)
        .expect("derived_dist")
        .value;
    assert!(
        (derived_val - 0.2).abs() < 1e-10,
        "derived_dist should be 0.2 after expression propagation, got {derived_val}"
    );

    // Extrusion with height=0.2 — BRep should exist and be a box.
    let brep = tree.current_brep().expect("brep");
    assert_eq!(brep.faces.len(), 6, "extruded box should have 6 faces");
}

// ---------------------------------------------------------------------------
// Test 5: Transaction rollback — BRep restored on feature failure
// ---------------------------------------------------------------------------

#[test]
fn transaction_rollback_restores_brep_on_failure() {
    let mut tree = FeatureTree::new();

    // Sketch produces faces; BRep has geometry after step 1.
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch 1",
    );

    let dist = tree
        .parameters_mut()
        .add_global("bad_dist", 1.0, Unit::Meters)
        .expect("add param");

    // Extrude references a bogus OperationId (99) that won't exist.
    let bad_extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(99, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist,
        },
        "Bad Extrude",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    // The sketch should succeed.
    assert!(
        matches!(
            result.feature_results.get(&sketch_id),
            Some(FeatureResult::Success { .. })
        ),
        "sketch should succeed"
    );

    // The bad extrude should fail.
    assert!(
        matches!(
            result.feature_results.get(&bad_extrude_id),
            Some(FeatureResult::Failed { .. })
        ),
        "bad extrude should fail"
    );

    // After rollback, current_brep should reflect the sketch state (pre-extrude),
    // not an empty BRep and not a partially-corrupted one.
    let brep = tree.current_brep().expect("brep after failed extrude");
    // The sketch produced a flat rectangle — at least 1 face, fewer than 6.
    assert!(
        !brep.faces.is_empty(),
        "brep should retain sketch geometry after rollback"
    );
    assert!(
        brep.faces.len() < 6,
        "brep should not have box faces after failed extrude, got {}",
        brep.faces.len()
    );
}
