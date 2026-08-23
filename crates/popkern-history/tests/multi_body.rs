//! Multi-body replay and Boolean feature wiring integration tests.
//!
//! Tests that the feature tree correctly maintains per-body BReps during
//! replay and that Boolean features consume the tool body as expected.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

// -- Helpers ------------------------------------------------------------------

/// Solved 1x1 rectangle sketch serialized to JSON.
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

/// Generated face at `op_id` index 0 (sketch face).
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

fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

/// Sketch + extrude tree. Returns (tree, sketch_id, extrude_id).
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

// -- 1. Two independent sketch+extrudes produce two bodies -------------------

#[test]
fn two_bodies_with_separate_breps() {
    let mut tree = FeatureTree::new();

    // Body 0: sketch + extrude on XY plane.
    let sk1 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    tree.replay(Tolerance::default()).expect("replay sk1");
    let sk1_face = sketch_face_name(1);

    let dist_ex1 = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    let ex1 = tree.add_feature(
        Feature::Extrude {
            face: sk1_face,
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_ex1 },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );

    // Body 1: sketch + extrude on XZ plane.
    let sk2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch2",
    );
    tree.replay(Tolerance::default()).expect("replay sk2");
    let sk2_face = sketch_face_name(3);

    let dist_ex2 = add_param(&mut tree, "_extrude2_distance", 1.0, Unit::Meters);

    let ex2 = tree.add_feature(
        Feature::Extrude {
            face: sk2_face,
            direction: Some(Vec3Desc::Literal([0.0, 1.0, 0.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_ex2 },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude2",
    );

    // Verify different bodies.
    let body0 = tree.feature(sk1).unwrap().body;
    let body1 = tree.feature(sk2).unwrap().body;
    assert_ne!(body0, body1, "two sketches should create separate bodies");
    assert_eq!(tree.feature(ex1).unwrap().body, body0);
    assert_eq!(tree.feature(ex2).unwrap().body, body1);

    // Replay and verify both body BReps are available.
    let result = tree.replay(Tolerance::default()).expect("full replay");
    assert!(
        result
            .feature_results
            .values()
            .all(|r| matches!(r, FeatureResult::Success { .. })),
        "all features should succeed"
    );

    let brep0 = tree.body_brep(body0);
    let brep1 = tree.body_brep(body1);
    assert!(brep0.is_some(), "body 0 should have a BRep");
    assert!(brep1.is_some(), "body 1 should have a BRep");
    assert_eq!(brep0.unwrap().faces.len(), 6, "body 0 should be a box");
    assert_eq!(brep1.unwrap().faces.len(), 6, "body 1 should be a box");
}

// -- 2. Boolean union of two bodies -----------------------------------------

#[test]
fn boolean_union_of_two_bodies() {
    let mut tree = FeatureTree::new();

    // Body 0: sketch + extrude on XY.
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    tree.replay(Tolerance::default()).expect("replay sk1");
    let sk1_face = sketch_face_name(1);

    let dist_ex1 = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);

    tree.add_feature(
        Feature::Extrude {
            face: sk1_face,
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_ex1 },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );

    // Body 1: sketch + extrude on XZ.
    let sk2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch2",
    );
    tree.replay(Tolerance::default()).expect("replay sk2");
    let sk2_face = sketch_face_name(3);

    let dist_ex2 = add_param(&mut tree, "_extrude2_distance", 1.0, Unit::Meters);

    tree.add_feature(
        Feature::Extrude {
            face: sk2_face,
            direction: Some(Vec3Desc::Literal([0.0, 1.0, 0.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_ex2 },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude2",
    );

    let body1 = tree.feature(sk2).unwrap().body;

    // Boolean union: merge body 1 into body 0.
    let bool_id = tree.add_feature(
        Feature::Boolean {
            operation: BoolOpKind::Union,
            tool_body: body1,
        },
        "BoolUnion",
    );

    let result = tree.replay(Tolerance::default()).expect("boolean replay");

    // The boolean feature itself should succeed.
    match result.feature_results.get(&bool_id) {
        Some(FeatureResult::Success { .. }) => {}
        other => panic!("expected Boolean Success, got: {other:?}"),
    }

    // Both extrudes build the SAME [0,1]^3 cube (the shared unit
    // square on XY extruded +z, and on XZ extruded +y), so the
    // correct union is that cube: exactly 6 faces and a valid body.
    // The historical `> 6 faces` assertion pinned the union of
    // GEOMETRICALLY BROKEN inputs — before PK-316/PK-317 the
    // replay-built boxes carried the base plane's surface on every
    // face and no side curves, and the boolean's output face count
    // reflected that garbage rather than the true shape.
    let brep = tree.current_brep().expect("should have BRep after boolean");
    brep.validate()
        .expect("union of coincident cubes validates");
    assert_eq!(
        brep.faces.len(),
        6,
        "union of two coincident unit cubes is the cube itself"
    );
}

// -- 3. Boolean with missing tool body fails gracefully ---------------------

#[test]
fn boolean_missing_tool_body_produces_failed() {
    let (mut tree, _sk, _ex) = build_sketch_extrude(1.0);

    let bool_id = tree.add_feature(
        Feature::Boolean {
            operation: BoolOpKind::Subtraction,
            tool_body: BodyId(99),
        },
        "BadBool",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    match result.feature_results.get(&bool_id) {
        Some(FeatureResult::Failed { error }) => {
            assert!(
                error.contains("not found"),
                "should mention tool body not found, got: {error}"
            );
        }
        other => panic!("expected Failed, got: {other:?}"),
    }

    // The box should still be intact (last good state).
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);
}
