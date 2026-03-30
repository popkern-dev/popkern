// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

use crate::FeatureTree;
use crate::error::FeatureResult;
use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};
use crate::node::FeatureId;

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

    // Fix all four corners so the sketch is fully constrained.
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));

    sk.solve().expect("rectangle should solve");

    serde_json::to_string(&sk).expect("rectangle should serialize")
}

/// Helper: create an EntityName referencing a given OperationId.
fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Test 1: Sketch replay produces a valid BRep
// ---------------------------------------------------------------------------

#[test]
fn sketch_replay_produces_faces() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let result = tree
        .replay(Tolerance::default())
        .expect("replay should succeed");

    // The sketch should succeed.
    for (_, outcome) in &result.feature_results {
        assert!(
            matches!(outcome, FeatureResult::Success { .. }),
            "sketch should succeed, got: {outcome:?}"
        );
    }

    // The resulting BRep should have geometry (a rectangle has 2 faces).
    let brep = tree.current_brep().expect("brep should exist");
    assert!(brep.faces.len() >= 1, "expected at least 1 face");
    assert!(brep.edges.len() >= 4, "expected at least 4 edges");
    assert!(brep.vertices.len() >= 4, "expected at least 4 vertices");
}

// ---------------------------------------------------------------------------
// Test 2: Sketch + Extrude replay produces a box
// ---------------------------------------------------------------------------

#[test]
fn sketch_and_extrude_replay_produces_box() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    // Add sketch — gets OperationId(1).
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Replay sketch first to find out what face names it produces.
    let _r = tree.replay(Tolerance::default()).expect("sketch replay");
    let _brep = tree.current_brep().expect("brep after sketch");

    // The sketch assigns Generated face names with OperationId(1).
    // Find a face to extrude.
    let face_name = entity_name(1, RoleCategory::Generated, 0);

    // Add extrude — gets OperationId(2).
    tree.add_feature(
        Feature::Extrude {
            face: face_name,
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: 1.0,
        },
        "Extrude 1",
    );

    // Replay again with both features.
    let result = tree.replay(Tolerance::default()).expect("full replay");

    // Both features should succeed.
    let successes = result
        .feature_results
        .values()
        .filter(|r| matches!(r, FeatureResult::Success { .. }))
        .count();
    assert_eq!(successes, 2, "both features should succeed");

    // An extruded rectangle should have 6 faces (box).
    let brep = tree.current_brep().expect("brep after extrude");
    assert_eq!(brep.faces.len(), 6, "extruded rectangle = 6-face box");
}

// ---------------------------------------------------------------------------
// Test 3: Suppressed feature is skipped
// ---------------------------------------------------------------------------

#[test]
fn suppressed_feature_is_skipped() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: 1.0,
        },
        "Extrude 1",
    );

    // Suppress the extrude.
    tree.features[extrude_id].suppressed = true;

    let result = tree.replay(Tolerance::default()).expect("replay");

    assert!(
        matches!(
            result.feature_results.get(&extrude_id),
            Some(FeatureResult::Suppressed)
        ),
        "extrude should be suppressed"
    );

    assert!(
        matches!(
            result.feature_results.get(&sketch_id),
            Some(FeatureResult::Success { .. })
        ),
        "sketch should succeed"
    );

    // Only the sketch BRep — no extrude means a flat face, not a box.
    let brep = tree.current_brep().expect("brep");
    assert!(
        brep.faces.len() < 6,
        "suppressed extrude means fewer than 6 faces"
    );
}

// ---------------------------------------------------------------------------
// Test 4: Missing reference produces Failed result
// ---------------------------------------------------------------------------

#[test]
fn missing_reference_produces_failed() {
    let mut tree = FeatureTree::new();

    // Add an extrude that references a bogus entity name (OperationId 99
    // does not exist in any BRep produced so far).
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(99, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: 1.0,
        },
        "Bad Extrude",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    let outcomes: Vec<_> = result.feature_results.values().collect();
    assert_eq!(outcomes.len(), 1);
    assert!(
        matches!(outcomes[0], FeatureResult::Failed { .. }),
        "bogus reference should produce Failed result, got: {outcomes:?}"
    );
}

// ---------------------------------------------------------------------------
// Test 5: Replay returns positive elapsed time
// ---------------------------------------------------------------------------

#[test]
fn replay_returns_positive_elapsed() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(
        result.elapsed > std::time::Duration::ZERO,
        "elapsed should be > 0"
    );
}

// ---------------------------------------------------------------------------
// Test 6: replay_from validates feature id
// ---------------------------------------------------------------------------

#[test]
fn replay_from_invalid_id_fails() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Use a default FeatureId that doesn't exist in the tree.
    let bogus_id = FeatureId::default();
    let result = tree.replay_from(bogus_id, Tolerance::default());
    assert!(result.is_err(), "replay_from with bogus id should fail");
}

// ---------------------------------------------------------------------------
// Test 7: Boolean returns not-yet-wired error
// ---------------------------------------------------------------------------

#[test]
fn boolean_feature_returns_not_wired() {
    let mut tree = FeatureTree::new();

    tree.add_feature(
        Feature::Boolean {
            operation: crate::feature::BoolOpKind::Union,
            tool_body: crate::node::BodyId(1),
        },
        "Union",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");
    let outcomes: Vec<_> = result.feature_results.values().collect();
    assert!(
        matches!(outcomes[0], FeatureResult::Failed { error } if error.contains("multi-body")),
        "boolean should fail with multi-body not wired"
    );
}

// ---------------------------------------------------------------------------
// Test 8: current_brep is None before replay
// ---------------------------------------------------------------------------

#[test]
fn current_brep_none_before_replay() {
    let tree = FeatureTree::new();
    assert!(tree.current_brep().is_none());
}
