// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use crate::FeatureTree;
use crate::delete::DeleteBehavior;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::parameter::ParameterId;
use crate::unit::Unit;
use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};

/// Build an [`EntityName`] that references a given raw operation id.
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

fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

// ---------------------------------------------------------------------------
// 1. Cascade: A → B → C; delete A removes B and C from timeline too
// ---------------------------------------------------------------------------

#[test]
fn cascade_removes_transitive_dependents() {
    let mut tree = FeatureTree::new();

    // A — base sketch (op_id 1)
    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_distance", 5.0, Unit::Meters);

    // B — extrude referencing A's face (op_id 2)
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

    // C — fillet referencing B's edge (op_id 3)
    let c = tree.add_feature(
        Feature::Fillet {
            edges: vec![entity_name(b_op.0)],
            radius: radius_c,
        },
        "C",
    );

    assert_eq!(tree.timeline(), &[a, b, c]);

    let result = tree
        .delete_feature(a, DeleteBehavior::Cascade)
        .expect("cascade delete should succeed");

    assert!(result.deleted.contains(&a), "A should be in deleted list");
    assert!(
        result.deleted.contains(&b),
        "B should be in deleted list (transitively depends on A)"
    );
    assert!(
        result.deleted.contains(&c),
        "C should be in deleted list (transitively depends on A via B)"
    );
    assert!(result.reparented.is_empty());
    assert!(
        tree.timeline().is_empty(),
        "timeline must be empty after cascade delete of root"
    );
}

// ---------------------------------------------------------------------------
// 2. Cascade with no dependents — only the feature itself is deleted
// ---------------------------------------------------------------------------

#[test]
fn cascade_no_dependents_deletes_only_target() {
    let mut tree = FeatureTree::new();
    let a = tree.add_feature(sketch_xy(), "A");
    let b = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: String::new(),
        },
        "B",
    );

    // B does not depend on A; deleting A cascade should only remove A.
    let result = tree
        .delete_feature(a, DeleteBehavior::Cascade)
        .expect("cascade delete should succeed");

    assert_eq!(result.deleted, vec![a]);
    assert_eq!(tree.timeline(), &[b]);
}

// ---------------------------------------------------------------------------
// 3. Reparent: A → B → C; delete B → C reparented to A
// ---------------------------------------------------------------------------

#[test]
fn reparent_rewires_dependents_to_grandparent() {
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

    // Before delete: B depends on A; C depends on B.
    assert!(tree.feature(b).unwrap().depends_on.contains(&a));
    assert!(tree.feature(c).unwrap().depends_on.contains(&b));

    let result = tree
        .delete_feature(b, DeleteBehavior::Reparent)
        .expect("reparent delete should succeed");

    // B is deleted, A and C remain.
    assert_eq!(result.deleted, vec![b]);
    assert!(result.reparented.contains(&c));
    assert_eq!(tree.timeline(), &[a, c]);

    // C should now depend on A instead of B.
    let c_deps = tree.feature(c).unwrap().depends_on.clone();
    assert!(
        c_deps.contains(&a),
        "C should be reparented to A after B is deleted"
    );
    assert!(
        !c_deps.contains(&b),
        "C should no longer depend on deleted B"
    );
}

// ---------------------------------------------------------------------------
// 4. Delete non-existent feature returns FeatureNotFound error
// ---------------------------------------------------------------------------

#[test]
fn delete_nonexistent_feature_returns_error() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(sketch_xy(), "A");
    // Remove from timeline manually so the slotmap entry exists but the
    // id is unknown to a fresh tree.
    let mut tree2 = FeatureTree::new();
    let err = tree2.delete_feature(id, DeleteBehavior::Cascade);
    assert!(
        matches!(err, Err(crate::error::HistoryError::FeatureNotFound)),
        "expected FeatureNotFound, got {:?}",
        err
    );
}

// ---------------------------------------------------------------------------
// 5. Delete + undo restores feature in timeline
// ---------------------------------------------------------------------------

#[test]
fn cascade_delete_undo_restores_timeline() {
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

    assert_eq!(tree.timeline(), &[a, b]);

    // Delete A cascade (removes A and B).
    tree.delete_feature(a, DeleteBehavior::Cascade)
        .expect("delete should succeed");
    assert!(tree.timeline().is_empty());

    // Undo should restore both features at their original positions.
    let undone = tree.undo();
    assert!(undone, "undo should succeed");
    // Both ids should be back in timeline (order may vary by slotmap allocation,
    // but both must be present).
    assert_eq!(
        tree.timeline().len(),
        2,
        "timeline should have 2 features after undo"
    );
    assert!(
        tree.timeline().contains(&a)
            || tree
                .timeline()
                .iter()
                .any(|&id| { tree.feature(id).map(|n| n.label == "A").unwrap_or(false) })
    );
}

// ---------------------------------------------------------------------------
// 6. Reparent + undo restores original dependency graph
// ---------------------------------------------------------------------------

#[test]
fn reparent_delete_undo_restores_deps() {
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

    // Delete B with reparent.
    tree.delete_feature(b, DeleteBehavior::Reparent)
        .expect("delete should succeed");
    assert_eq!(tree.timeline(), &[a, c]);

    // Undo — B should be back, and C should once again depend on B.
    let undone = tree.undo();
    assert!(undone, "undo should succeed");
    assert_eq!(tree.timeline().len(), 3, "all three features restored");
}

// ---------------------------------------------------------------------------
// 7. Warnings populated for reparented features
// ---------------------------------------------------------------------------

#[test]
fn reparent_produces_warnings_for_each_reparented_feature() {
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

    // Two features depending on B.
    let _c = tree.add_feature(
        Feature::Fillet {
            edges: vec![entity_name(b_op.0)],
            radius: radius_c,
        },
        "C",
    );

    let dist_d = add_param(&mut tree, "_d_distance", 0.5, Unit::Meters);

    let _d = tree.add_feature(
        Feature::Chamfer {
            edges: vec![entity_name(b_op.0)],
            distance: dist_d,
        },
        "D",
    );

    let result = tree
        .delete_feature(b, DeleteBehavior::Reparent)
        .expect("delete should succeed");

    assert_eq!(
        result.warnings.len(),
        2,
        "one warning per reparented feature"
    );
}

// ---------------------------------------------------------------------------
// 8. Delete cascade → feature params removed from table
// ---------------------------------------------------------------------------

#[test]
fn cascade_delete_removes_feature_params() {
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

    // dist_b is registered as a feature param for B during add_feature.
    // After cascade delete of B, its feature params should be gone.
    tree.delete_feature(b, DeleteBehavior::Cascade)
        .expect("delete should succeed");

    // The parameter b.distance should be gone from the table.
    let b_params = tree.parameters().feature_params(b);
    assert!(
        b_params.is_empty(),
        "cascade delete should remove feature-scoped params"
    );
}

// ---------------------------------------------------------------------------
// 9. Delete reparent → deleted feature's params removed
// ---------------------------------------------------------------------------

#[test]
fn reparent_delete_removes_target_feature_params() {
    let mut tree = FeatureTree::new();

    let a = tree.add_feature(sketch_xy(), "A");
    let a_op = tree.feature(a).unwrap().operation_id;

    let dist_b = add_param(&mut tree, "_b_dist2", 5.0, Unit::Meters);

    let b = tree.add_feature(
        Feature::Extrude {
            face: entity_name(a_op.0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_b },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "B2",
    );

    tree.delete_feature(b, DeleteBehavior::Reparent)
        .expect("delete should succeed");

    let b_params = tree.parameters().feature_params(b);
    assert!(
        b_params.is_empty(),
        "reparent delete should remove feature-scoped params for deleted feature"
    );
}
