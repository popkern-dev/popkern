// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use crate::FeatureTree;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::parameter::ParameterId;
use crate::unit::Unit;
use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};

/// Build an [`EntityName`] referencing a given raw operation id.
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

fn sketch_xy() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XY,
        sketch_data: String::new(),
    }
}

fn sketch_xz() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XZ,
        sketch_data: String::new(),
    }
}

fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

// ---------------------------------------------------------------------------
// 1. dependents(A) returns B and C when A → B → C
// ---------------------------------------------------------------------------

#[test]
fn dependents_returns_transitive_downstream() {
    let mut tree = FeatureTree::new();

    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    let b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );
    let b_op = tree.feature(b).unwrap().operation_id;

    let radius_c = add_param(&mut tree, "_c_radius", 1.0, Unit::Meters);

    let c = tree.add_feature(
        Feature::Fillet {
            edges: vec![entity_name(b_op.0)],
            radius: radius_c,
        },
        "C",
    );

    let deps = tree.dependents(a);
    assert!(deps.contains(&b), "B should be a dependent of A");
    assert!(deps.contains(&c), "C should be a transitive dependent of A");
    assert!(
        !deps.contains(&a),
        "A should not appear in its own dependents"
    );
}

// ---------------------------------------------------------------------------
// 2. dependents returns empty for a leaf (no downstream consumers)
// ---------------------------------------------------------------------------

#[test]
fn dependents_empty_for_leaf_feature() {
    let mut tree = FeatureTree::new();
    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    let b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );

    // B has no dependents.
    assert!(tree.dependents(b).is_empty());
}

// ---------------------------------------------------------------------------
// 3. dependents returns in timeline order
// ---------------------------------------------------------------------------

#[test]
fn dependents_in_timeline_order() {
    let mut tree = FeatureTree::new();

    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    let b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );
    let b_op = tree.feature(b).unwrap().operation_id;

    let radius_c = add_param(&mut tree, "_c_radius", 1.0, Unit::Meters);

    let c = tree.add_feature(
        Feature::Fillet {
            edges: vec![entity_name(b_op.0)],
            radius: radius_c,
        },
        "C",
    );

    let deps = tree.dependents(a);
    // B appears before C in the timeline so it should be first.
    let b_pos = deps.iter().position(|&id| id == b).unwrap();
    let c_pos = deps.iter().position(|&id| id == c).unwrap();
    assert!(b_pos < c_pos, "B should come before C in dependency list");
}

// ---------------------------------------------------------------------------
// 4. dependencies(B) returns direct parents only
// ---------------------------------------------------------------------------

#[test]
fn dependencies_returns_direct_parents() {
    let mut tree = FeatureTree::new();

    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    let b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );
    let b_op = tree.feature(b).unwrap().operation_id;

    let radius_c = add_param(&mut tree, "_c_radius", 1.0, Unit::Meters);

    let c = tree.add_feature(
        Feature::Fillet {
            edges: vec![entity_name(b_op.0)],
            radius: radius_c,
        },
        "C",
    );

    // C directly depends on B; A should NOT appear in C's direct deps.
    let c_deps = tree.dependencies(c);
    assert!(c_deps.contains(&b), "C directly depends on B");
    assert!(!c_deps.contains(&a), "A is not a direct dependency of C");

    // B directly depends on A.
    let b_deps = tree.dependencies(b);
    assert!(b_deps.contains(&a));
}

// ---------------------------------------------------------------------------
// 5. dependencies returns empty for a root feature
// ---------------------------------------------------------------------------

#[test]
fn dependencies_empty_for_root() {
    let mut tree = FeatureTree::new();
    let a = tree.add_feature(sketch_xy(), "A");
    assert!(tree.dependencies(a).is_empty());
}

// ---------------------------------------------------------------------------
// 6. dependencies returns empty for unknown feature
// ---------------------------------------------------------------------------

#[test]
fn dependencies_empty_for_unknown_id() {
    let tree1 = FeatureTree::new();
    let mut tree2 = FeatureTree::new();
    let id = tree2.add_feature(sketch_xy(), "A");
    // `id` is from tree2; tree1 does not know about it.
    assert!(tree1.dependencies(id).is_empty());
}

// ---------------------------------------------------------------------------
// 7. impact_analysis will_rebuild = transitive dependents
// ---------------------------------------------------------------------------

#[test]
fn impact_analysis_will_rebuild_matches_dependents() {
    let mut tree = FeatureTree::new();

    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    let b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );
    let b_op = tree.feature(b).unwrap().operation_id;

    let radius_c = add_param(&mut tree, "_c_radius", 1.0, Unit::Meters);

    let c = tree.add_feature(
        Feature::Fillet {
            edges: vec![entity_name(b_op.0)],
            radius: radius_c,
        },
        "C",
    );

    let report = tree.impact_analysis(a);
    let rebuilds = &report.will_rebuild;
    assert!(rebuilds.contains(&b));
    assert!(rebuilds.contains(&c));
    assert!(!rebuilds.contains(&a));
}

// ---------------------------------------------------------------------------
// 8. impact_analysis on a feature with no dependents → empty will_rebuild
// ---------------------------------------------------------------------------

#[test]
fn impact_analysis_no_dependents_empty_rebuild() {
    let mut tree = FeatureTree::new();
    let a = tree.add_feature(sketch_xy(), "A");
    let report = tree.impact_analysis(a);
    assert!(report.will_rebuild.is_empty());
    assert!(report.at_risk.is_empty());
}

// ---------------------------------------------------------------------------
// 9. impact_analysis at_risk is empty when no last_result is cached
// ---------------------------------------------------------------------------

#[test]
fn impact_analysis_at_risk_empty_without_cached_result() {
    let mut tree = FeatureTree::new();

    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    // B depends on A but A has no cached result (never replayed).
    let _b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );

    let report = tree.impact_analysis(a);
    assert!(
        report.at_risk.is_empty(),
        "at_risk should be empty when feature has no cached result"
    );
}

// ---------------------------------------------------------------------------
// 10. Multiple independent roots — dependents does not cross-contaminate
// ---------------------------------------------------------------------------

#[test]
fn dependents_does_not_cross_roots() {
    let mut tree = FeatureTree::new();

    // Root 1
    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    let _b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B",
    );

    // Root 2 — independent
    let c = tree.add_feature(sketch_xz(), "C");

    // C has no dependents; A's dependents should not include C.
    let a_deps = tree.dependents(a);
    assert!(
        !a_deps.contains(&c),
        "C is independent, should not be a dependent of A"
    );

    let c_deps = tree.dependents(c);
    assert!(c_deps.is_empty(), "C has no dependents");
}
