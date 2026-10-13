//! Cached replays republish onto the feature nodes exactly what an
//! execution would (PK-805): the republish skips its copy only when
//! the node already holds that step's result.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

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

// ---------------------------------------------------------------------------
// Test: a cached replay restores the node's last_result (PK-805)
// ---------------------------------------------------------------------------

#[test]
fn returning_to_a_cached_distance_restores_the_extrudes_last_result() {
    // Distance 1 -> 0 -> 1: the zero extrude declines; the third replay
    // resumes the distance-1 step from the step cache and must
    // republish its Success onto the node. The republish skips the copy
    // only when the node already holds that step's result; after the
    // declined execution it does not.
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let distance = add_param(&mut tree, "_extrude1_distance", 1.0, Unit::Meters);
    let extrude = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );
    tree.replay(Tolerance::default())
        .expect("replay distance 1");
    let last = |tree: &FeatureTree| {
        format!(
            "{:?}",
            tree.feature(extrude).map(|n| n.last_result.as_ref())
        )
    };
    let at_one = last(&tree);
    tree.set_parameter_value(distance, 0.0, Unit::Meters)
        .expect("distance 0");
    let _ = tree.replay(Tolerance::default());
    let at_zero = last(&tree);
    assert_ne!(at_zero, at_one, "a zero extrude changes the last result");
    tree.set_parameter_value(distance, 1.0, Unit::Meters)
        .expect("distance 1 again");
    tree.replay(Tolerance::default())
        .expect("replay distance 1");
    assert_eq!(
        last(&tree),
        at_one,
        "the cached distance-1 result is republished"
    );
}
