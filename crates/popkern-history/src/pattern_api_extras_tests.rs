// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//! Tests for the additional pattern convenience API methods.
#[cfg(test)]
use crate::FeatureTree;
use crate::feature::{AxisDesc, Vec3Desc};
use crate::node::FeatureId;
use crate::pattern::PatternDirection;
use crate::unit::Unit;

/// Helper: add a sketch + extrude pair and return `(sketch_id, extrude_id)`.
fn add_seed(tree: &mut FeatureTree, tag: &str) -> (FeatureId, FeatureId) {
    use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc};

    let dist_name = format!("_dist_extras_{tag}");
    let dist = tree
        .parameter_table
        .add_global(&dist_name, 0.005, Unit::Meters)
        .expect("add dist param");

    let sketch = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
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
fn add_linear_pattern_with_total_creates_feature() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "tot1");

    let pattern_id = tree
        .add_linear_pattern_with_total(
            vec![sketch, extrude],
            Vec3Desc::Literal([1.0, 0.0, 0.0]),
            5.0,
            40.0, // 40 mm total distance
            "lin_total",
        )
        .expect("add_linear_pattern_with_total");

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(node.depends_on.contains(&extrude));

    // primary_count should be 5.
    let count_id = tree
        .parameters()
        .resolve_name("lin_total.primary_count")
        .expect("primary_count param");
    let count_val = tree.parameters().get(count_id).expect("param").value;
    assert!((count_val - 5.0).abs() < 1e-9, "count should be 5");

    // primary_total_distance = 40 mm = 0.040 m.
    let total_id = tree
        .parameters()
        .resolve_name("lin_total.primary_total_distance")
        .expect("total_distance param");
    let total_val = tree.parameters().get(total_id).expect("param").value;
    assert!(
        (total_val - 0.040).abs() < 1e-9,
        "total_distance should be 0.040 m, got {total_val}"
    );

    // No expression should be set (plain literal).
    let has_expr = tree
        .parameters()
        .get(total_id)
        .expect("param")
        .expression
        .is_some();
    assert!(!has_expr, "total_distance should have no expression");
}

#[test]
fn add_circular_pattern_with_span_creates_feature() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "span1");

    let half_pi = std::f64::consts::FRAC_PI_2;
    let pattern_id = tree
        .add_circular_pattern_with_span(
            vec![sketch, extrude],
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            3.0,
            half_pi, // 90 degrees total span
            "circ_span",
        )
        .expect("add_circular_pattern_with_span");

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(node.depends_on.contains(&extrude));

    let angle_id = tree
        .parameters()
        .resolve_name("circ_span.total_angle")
        .expect("total_angle param");
    let angle_val = tree.parameters().get(angle_id).expect("param").value;
    assert!(
        (angle_val - half_pi).abs() < 1e-12,
        "total_angle should be PI/2, got {angle_val}"
    );

    // No expression on total_angle (plain literal).
    let has_expr = tree
        .parameters()
        .get(angle_id)
        .expect("param")
        .expression
        .is_some();
    assert!(!has_expr, "total_angle should have no expression");
}

#[test]
fn add_curve_pattern_creates_feature() {
    use popkern_core::naming::{EntityName, EntityRole, NameId, RoleCategory};

    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "curve1");

    // Create a dummy EntityName to use as the curve reference.
    let (dummy_sketch, _) = add_seed(&mut tree, "dummy");
    let dummy_node = tree.feature(dummy_sketch).expect("dummy sketch node");
    let curve = EntityName {
        operation: dummy_node.operation_id,
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: vec![NameId(0)],
    };

    let pattern_id = tree
        .add_curve_pattern(vec![sketch, extrude], curve, 6.0, true, false, "curv_pat")
        .expect("add_curve_pattern");

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(node.depends_on.contains(&extrude));

    let count_id = tree
        .parameters()
        .resolve_name("curv_pat.count")
        .expect("count param");
    let count_val = tree.parameters().get(count_id).expect("param").value;
    assert!((count_val - 6.0).abs() < 1e-9, "count should be 6");
}

#[test]
fn skip_pattern_instance_2d_computes_flat_index() {
    use crate::feature::Feature;
    use crate::pattern::PatternDirection;

    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "grid1");

    // Manually build a 2D linear pattern (3×4 grid).
    let primary_count = tree
        .parameter_table
        .add_global("grid1_pat.primary_count", 3.0, Unit::Ratio)
        .expect("primary_count");
    let primary_total = tree
        .parameter_table
        .add_global("grid1_pat.primary_total_distance", 0.020, Unit::Meters)
        .expect("primary_total");
    let secondary_count = tree
        .parameter_table
        .add_global("grid1_pat.secondary_count", 4.0, Unit::Ratio)
        .expect("secondary_count");
    let secondary_total = tree
        .parameter_table
        .add_global("grid1_pat.secondary_total_distance", 0.030, Unit::Meters)
        .expect("secondary_total");

    let pat_feature = Feature::LinearPattern {
        seed_features: vec![sketch, extrude],
        primary: PatternDirection {
            direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
            count: primary_count,
            total_distance: primary_total,
        },
        secondary: Some(PatternDirection {
            direction: Vec3Desc::Literal([0.0, 1.0, 0.0]),
            count: secondary_count,
            total_distance: secondary_total,
        }),
        geometry_pattern: false,
        instance_overrides: std::collections::HashMap::new(),
    };

    let pattern_id = tree.add_feature(pat_feature, "grid1_pat");

    // Suppress (row=1, col=2) → flat = 1*4 + 2 = 6.
    tree.skip_pattern_instance_2d(pattern_id, 1, 2);

    let node = tree.feature(pattern_id).expect("pattern node");
    match &node.feature {
        Feature::LinearPattern {
            instance_overrides, ..
        } => {
            assert!(
                instance_overrides.contains_key(&6),
                "flat index 6 should be suppressed"
            );
            assert!(
                instance_overrides[&6].suppressed,
                "instance at flat index 6 should be suppressed"
            );
        }
        _ => panic!("expected LinearPattern"),
    }
}

#[test]
fn skip_pattern_instance_2d_no_secondary_is_noop() {
    let mut tree = FeatureTree::new();
    let (sketch, extrude) = add_seed(&mut tree, "noop1");

    let pattern_id = tree
        .add_linear_pattern_with_total(
            vec![sketch, extrude],
            Vec3Desc::Literal([1.0, 0.0, 0.0]),
            3.0,
            20.0,
            "noop_pat",
        )
        .expect("create pattern");

    // Should be a no-op — no secondary direction.
    tree.skip_pattern_instance_2d(pattern_id, 0, 1);

    let node = tree.feature(pattern_id).expect("pattern node");
    match &node.feature {
        crate::feature::Feature::LinearPattern {
            instance_overrides, ..
        } => {
            assert!(
                instance_overrides.is_empty(),
                "no overrides expected when there is no secondary direction"
            );
        }
        _ => panic!("expected LinearPattern"),
    }
}
