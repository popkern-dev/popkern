// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//! Tests for the pattern convenience API.
#[cfg(test)]
use crate::FeatureTree;
use crate::feature::{AxisDesc, Vec3Desc};
use crate::pattern::MirrorPlaneDesc;
use crate::unit::Unit;

/// Helper: add a sketch + extrude pair and return `(sketch_id, extrude_id)`.
///
/// Pattern seed lists must start with a Sketch and contain a Boss/Cut feature,
/// so callers should pass `vec![sketch_id, extrude_id]` as the seed list.
fn add_seed(tree: &mut FeatureTree, tag: &str) -> (crate::node::FeatureId, crate::node::FeatureId) {
    use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc};

    let dist_name = format!("_dist_{tag}");
    let dist = tree
        .parameter_table
        .add_global(&dist_name, 0.005, Unit::Meters)
        .expect("add dist param");

    let sketch = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(String::new()),
        },
        &format!("Sketch_{tag}"),
    );
    let extrude = tree.add_feature(
        Feature::Extrude {
            face: {
                use popkern_core::naming::{EntityName, EntityRole, NameId, RoleCategory};
                let sketch_node = tree.feature(sketch).expect("sketch node");
                EntityName {
                    operation: sketch_node.operation_id,
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![NameId(0)],
                }
            },
            direction: None,
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::Boss,
            draft_angle: None,
            second_direction: None,
        },
        &format!("Extrude_{tag}"),
    );
    (sketch, extrude)
}

#[test]
fn add_linear_pattern_with_pitch_creates_feature() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "a");

    let pattern_id = tree
        .add_linear_pattern_with_pitch(
            vec![sketch, extrude],
            Vec3Desc::Literal([1.0, 0.0, 0.0]),
            4.0,
            10.0, // 10 mm pitch
            "lin1",
        )
        .expect("add_linear_pattern_with_pitch");

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(node.depends_on.contains(&extrude), "should depend on seed");

    // Parameters should exist in the table.
    assert!(
        tree.parameters()
            .resolve_name("lin1.primary_count")
            .is_some(),
        "primary_count param should be registered"
    );
    assert!(
        tree.parameters()
            .resolve_name("lin1.primary_total_distance")
            .is_some(),
        "primary_total_distance param should be registered"
    );

    // total_distance = pitch × (count - 1) = 10 mm × 3 = 30 mm = 0.030 m
    let total_id = tree
        .parameters()
        .resolve_name("lin1.primary_total_distance")
        .expect("total_distance param");
    let total_val = tree.parameters().get(total_id).expect("param").value;
    assert!(
        (total_val - 0.030).abs() < 1e-9,
        "total_distance should be 0.030 m, got {total_val}"
    );
}

#[test]
fn add_circular_pattern_full_circle_creates_feature() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "b");

    let pattern_id = tree
        .add_circular_pattern_full_circle(
            vec![sketch, extrude],
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            6.0,
            "circ_full",
        )
        .expect("add_circular_pattern_full_circle");

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(node.depends_on.contains(&extrude));

    let angle_id = tree
        .parameters()
        .resolve_name("circ_full.total_angle")
        .expect("total_angle param");
    let angle_val = tree.parameters().get(angle_id).expect("param").value;
    assert!(
        (angle_val - std::f64::consts::TAU).abs() < 1e-12,
        "total_angle should be TAU, got {angle_val}"
    );
}

#[test]
fn add_circular_pattern_with_pitch_creates_feature() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "c");
    let pitch = std::f64::consts::PI / 4.0; // 45 degrees

    let pattern_id = tree
        .add_circular_pattern_with_pitch(
            vec![sketch, extrude],
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            4.0,
            pitch,
            "circ_pitch",
        )
        .expect("add_circular_pattern_with_pitch");

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(node.depends_on.contains(&extrude));

    // total_angle = pitch × (count - 1) = PI/4 × 3 = 3PI/4
    let total_id = tree
        .parameters()
        .resolve_name("circ_pitch.total_angle")
        .expect("total_angle param");
    let total_val = tree.parameters().get(total_id).expect("param").value;
    let expected = pitch * 3.0;
    assert!(
        (total_val - expected).abs() < 1e-12,
        "total_angle should be {expected}, got {total_val}"
    );
}

#[test]
fn add_mirror_creates_feature() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "d");

    let mirror_id = tree
        .add_mirror(vec![sketch, extrude], MirrorPlaneDesc::XZ, true, "mir1")
        .expect("add_mirror");

    let node = tree.feature(mirror_id).expect("mirror node");
    assert!(node.depends_on.contains(&extrude));
    match &node.feature {
        crate::feature::Feature::Mirror { merge_result, .. } => {
            assert!(*merge_result);
        }
        _ => panic!("expected Mirror feature"),
    }
}

#[test]
fn skip_pattern_instance_sets_suppressed() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "e");

    let pattern_id = tree
        .add_circular_pattern_full_circle(
            vec![sketch, extrude],
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            5.0,
            "skip_test",
        )
        .expect("create pattern");

    tree.skip_pattern_instance(pattern_id, 2);

    let node = tree.feature(pattern_id).expect("pattern node");
    match &node.feature {
        crate::feature::Feature::CircularPattern {
            instance_overrides, ..
        } => {
            let ovr = instance_overrides.get(&2).expect("override at index 2");
            assert!(ovr.suppressed, "instance 2 should be suppressed");
        }
        _ => panic!("expected CircularPattern"),
    }
}

#[test]
fn override_pattern_instance_param_stores_value() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "f");

    let pattern_id = tree
        .add_linear_pattern_with_pitch(
            vec![sketch, extrude],
            Vec3Desc::Literal([0.0, 1.0, 0.0]),
            3.0,
            5.0,
            "ovr_test",
        )
        .expect("create pattern");

    // Use the count param id as a stand-in for the override target.
    let count_id = tree
        .parameters()
        .resolve_name("ovr_test.primary_count")
        .expect("count param");

    tree.override_pattern_instance_param(pattern_id, 1, count_id, 99.0);

    let node = tree.feature(pattern_id).expect("pattern node");
    match &node.feature {
        crate::feature::Feature::LinearPattern {
            instance_overrides, ..
        } => {
            let ovr = instance_overrides.get(&1).expect("override at index 1");
            assert_eq!(
                ovr.param_overrides.get(&count_id).copied(),
                Some(99.0),
                "param override should be 99.0"
            );
        }
        _ => panic!("expected LinearPattern"),
    }
}

#[test]
fn linear_pitch_total_distance_expression_re_evaluates_on_count_change() {
    use crate::expression::evaluate_all;
    use crate::measure_query::MeasureContext;

    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "expr_check");

    tree.add_linear_pattern_with_pitch(
        vec![sketch, extrude],
        Vec3Desc::Literal([1.0, 0.0, 0.0]),
        4.0,  // count = 4
        10.0, // pitch = 10 mm
        "expr_lin",
    )
    .expect("create pitch pattern");

    // Initial: total = 10mm * (4-1) = 30mm = 0.030 m
    let total_id = tree
        .parameters()
        .resolve_name("expr_lin.primary_total_distance")
        .expect("total_distance param");
    let total_val = tree.parameters().get(total_id).expect("param").value;
    assert!(
        (total_val - 0.030).abs() < 1e-9,
        "initial total should be 0.030 m, got {total_val}"
    );

    // Verify the expression is set.
    assert!(
        tree.parameters()
            .get(total_id)
            .expect("param")
            .expression
            .is_some(),
        "total_distance should have an expression"
    );

    // Change count to 6 and re-evaluate expressions.
    let count_id = tree
        .parameters()
        .resolve_name("expr_lin.primary_count")
        .expect("count param");
    tree.parameter_table
        .set_value(count_id, 6.0, Unit::Ratio)
        .expect("set count to 6");

    evaluate_all(&mut tree.parameter_table, &MeasureContext::empty());

    // total = 10mm * (6-1) = 50mm = 0.050 m
    let total_val_new = tree.parameters().get(total_id).expect("param").value;
    assert!(
        (total_val_new - 0.050).abs() < 1e-9,
        "after count=6: total should be 0.050 m, got {total_val_new}"
    );
}

#[test]
fn circular_pitch_total_angle_expression_re_evaluates_on_count_change() {
    use crate::expression::evaluate_all;
    use crate::measure_query::MeasureContext;

    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "circ_expr");

    let pitch = std::f64::consts::PI / 6.0; // 30 degrees
    tree.add_circular_pattern_with_pitch(
        vec![sketch, extrude],
        AxisDesc::Literal {
            origin: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
        4.0,   // count = 4
        pitch, // pitch = 30°
        "expr_circ",
    )
    .expect("create pitch circular pattern");

    // Initial: total = PI/6 * (4-1) = PI/2
    let total_id = tree
        .parameters()
        .resolve_name("expr_circ.total_angle")
        .expect("total_angle param");
    let total_val = tree.parameters().get(total_id).expect("param").value;
    assert!(
        (total_val - std::f64::consts::FRAC_PI_2).abs() < 1e-12,
        "initial total_angle should be PI/2, got {total_val}"
    );

    // Verify the expression is set.
    assert!(
        tree.parameters()
            .get(total_id)
            .expect("param")
            .expression
            .is_some(),
        "total_angle should have an expression"
    );

    // Change count to 7 and re-evaluate.
    let count_id = tree
        .parameters()
        .resolve_name("expr_circ.count")
        .expect("count param");
    tree.parameter_table
        .set_value(count_id, 7.0, Unit::Ratio)
        .expect("set count to 7");

    evaluate_all(&mut tree.parameter_table, &MeasureContext::empty());

    // total = PI/6 * (7-1) = PI
    let total_val_new = tree.parameters().get(total_id).expect("param").value;
    assert!(
        (total_val_new - std::f64::consts::PI).abs() < 1e-12,
        "after count=7: total_angle should be PI, got {total_val_new}"
    );
}

#[test]
fn duplicate_label_returns_error() {
    let mut tree = FeatureTree::new();
    let (sketch_g, extrude_g) = add_seed(&mut tree, "g");
    let (sketch_h, extrude_h) = add_seed(&mut tree, "h");

    tree.add_linear_pattern_with_pitch(
        vec![sketch_g, extrude_g],
        Vec3Desc::Literal([1.0, 0.0, 0.0]),
        2.0,
        5.0,
        "dup",
    )
    .expect("first insert");

    let result = tree.add_linear_pattern_with_pitch(
        vec![sketch_h, extrude_h],
        Vec3Desc::Literal([1.0, 0.0, 0.0]),
        2.0,
        5.0,
        "dup", // same label → same param names → duplicate
    );
    assert!(
        result.is_err(),
        "second insert with same label should fail with DuplicateParameterName"
    );
}
