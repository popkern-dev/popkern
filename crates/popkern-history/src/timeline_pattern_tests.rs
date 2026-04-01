// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
//! Tests for pattern seed-feature dependency inference and body assignment.
#[cfg(test)]
use super::*;
use crate::FeatureTree;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, Vec3Desc};
use crate::node::BodyId;
use crate::parameter::ParameterId;
use crate::pattern::{MirrorPlaneDesc, PatternDirection};
use crate::unit::Unit;
use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};

/// Helper: create an EntityName referencing a given OperationId.
fn entity_name(op_id: u64) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: vec![NameId(0)],
    }
}

/// Helper: add a global parameter and return its id.
fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

/// Helper: add a sketch + extrude pair and return both ids.
fn add_sketch_and_extrude(
    tree: &mut FeatureTree,
    op_id: u64,
) -> (crate::node::FeatureId, crate::node::FeatureId) {
    let dist = add_param(tree, &format!("_d{op_id}"), 5.0, Unit::Meters);
    let sketch = tree.add_feature(
        Feature::Sketch {
            plane: crate::feature::SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        &format!("Sketch{op_id}"),
    );
    let extrude = tree.add_feature(
        Feature::Extrude {
            face: entity_name(op_id),
            direction: None,
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::Boss,
            draft_angle: None,
            second_direction: None,
        },
        &format!("Extrude{op_id}"),
    );
    (sketch, extrude)
}

#[test]
fn linear_pattern_depends_on_seed_features() {
    let mut tree = FeatureTree::new();
    let (_sketch, extrude_id) = add_sketch_and_extrude(&mut tree, 1);
    let count = add_param(&mut tree, "_n", 3.0, Unit::Ratio);
    let total = add_param(&mut tree, "_t", 20.0, Unit::Meters);

    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![extrude_id],
            primary: PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count,
                total_distance: total,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPattern",
    );

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(
        node.depends_on.contains(&extrude_id),
        "linear pattern must depend on its seed feature"
    );
    assert_eq!(node.body, BodyId(0), "pattern body must be BodyId(0)");
}

#[test]
fn circular_pattern_depends_on_seed_features() {
    let mut tree = FeatureTree::new();
    let (_sketch, extrude_id) = add_sketch_and_extrude(&mut tree, 1);
    let count = add_param(&mut tree, "_n", 4.0, Unit::Ratio);
    let angle = add_param(&mut tree, "_a", std::f64::consts::TAU, Unit::Radians);

    let pattern_id = tree.add_feature(
        Feature::CircularPattern {
            seed_features: vec![extrude_id],
            axis: crate::feature::AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            count,
            total_angle: angle,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "CirPattern",
    );

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(
        node.depends_on.contains(&extrude_id),
        "circular pattern must depend on its seed feature"
    );
    assert_eq!(node.body, BodyId(0), "pattern body must be BodyId(0)");
}

#[test]
fn curve_pattern_depends_on_seed_and_curve() {
    let mut tree = FeatureTree::new();
    let (_sketch, extrude_id) = add_sketch_and_extrude(&mut tree, 1);
    let count = add_param(&mut tree, "_n", 5.0, Unit::Ratio);
    // op_id 3 for the curve edge (after sketch=op1, extrude=op2)
    let curve_name = entity_name(3);

    let pattern_id = tree.add_feature(
        Feature::CurvePattern {
            seed_features: vec![extrude_id],
            curve: curve_name,
            count,
            equal_spacing: true,
            align_to_curve: false,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "CrvPattern",
    );

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(
        node.depends_on.contains(&extrude_id),
        "curve pattern must depend on its seed feature"
    );
    assert_eq!(node.body, BodyId(0), "pattern body must be BodyId(0)");
}

#[test]
fn mirror_depends_on_seed_features() {
    let mut tree = FeatureTree::new();
    let (_sketch, extrude_id) = add_sketch_and_extrude(&mut tree, 1);

    let pattern_id = tree.add_feature(
        Feature::Mirror {
            seed_features: vec![extrude_id],
            plane: MirrorPlaneDesc::XY,
            merge_result: true,
        },
        "Mirror",
    );

    let node = tree.feature(pattern_id).expect("pattern node");
    assert!(
        node.depends_on.contains(&extrude_id),
        "mirror must depend on its seed feature"
    );
    assert_eq!(node.body, BodyId(0), "mirror body must be BodyId(0)");
}

#[test]
fn mirror_on_face_extracts_plane_entity_name() {
    use popkern_core::naming::OperationId;
    let face_name = entity_name(99);
    let f = Feature::Mirror {
        seed_features: vec![],
        plane: MirrorPlaneDesc::OnFace(face_name),
        merge_result: true,
    };
    let names = feature_entity_names(&f);
    assert_eq!(names.len(), 1);
    assert_eq!(names[0].operation, OperationId(99));
}

#[test]
fn pattern_no_duplicate_deps_when_seed_also_referenced_by_name() {
    // If a seed feature is also an entity-name dep (unlikely but possible),
    // the deduplication logic must not produce duplicate FeatureIds.
    let mut tree = FeatureTree::new();
    let (_sketch, extrude_id) = add_sketch_and_extrude(&mut tree, 1);
    let extrude_op = tree.feature(extrude_id).expect("extrude").operation_id;
    let count = add_param(&mut tree, "_n", 2.0, Unit::Ratio);
    let total = add_param(&mut tree, "_t", 10.0, Unit::Meters);

    // Direction references a face produced by the same extrude op.
    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![extrude_id],
            primary: PatternDirection {
                direction: Vec3Desc::FaceNormal(EntityName {
                    operation: extrude_op,
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: 0,
                    },
                    ancestors: vec![NameId(0)],
                }),
                count,
                total_distance: total,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPatternDedup",
    );

    let node = tree.feature(pattern_id).expect("pattern node");
    let dep_count = node.depends_on.iter().filter(|&&d| d == extrude_id).count();
    assert_eq!(
        dep_count, 1,
        "extrude_id must appear exactly once in depends_on"
    );
}
