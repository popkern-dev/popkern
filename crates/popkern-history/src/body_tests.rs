// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use crate::node::FeatureId;
use slotmap::SlotMap;

// ---------------------------------------------------------------------------
// BodyManager unit tests
// ---------------------------------------------------------------------------

#[test]
fn create_body_returns_sequential_ids() {
    let mut mgr = BodyManager::new();
    let b0 = mgr.create_body();
    let b1 = mgr.create_body();
    let b2 = mgr.create_body();
    assert_eq!(b0, BodyId(0));
    assert_eq!(b1, BodyId(1));
    assert_eq!(b2, BodyId(2));
}

#[test]
fn assign_feature_and_body_features_round_trip() {
    let mut mgr = BodyManager::new();
    let body = mgr.create_body();

    // Use a SlotMap to get real FeatureIds.
    let mut sm: SlotMap<FeatureId, u32> = SlotMap::with_key();
    let fid_a = sm.insert(1);
    let fid_b = sm.insert(2);

    mgr.assign_feature(body, fid_a);
    mgr.assign_feature(body, fid_b);

    let features = mgr.body_features(body);
    assert_eq!(features.len(), 2);
    assert!(features.contains(&fid_a));
    assert!(features.contains(&fid_b));
}

#[test]
fn body_features_unknown_body_returns_empty() {
    let mgr = BodyManager::new();
    let features = mgr.body_features(BodyId(99));
    assert!(features.is_empty());
}

#[test]
fn consume_body_removes_it() {
    let mut mgr = BodyManager::new();
    let body = mgr.create_body();

    let mut sm: SlotMap<FeatureId, u32> = SlotMap::with_key();
    let fid = sm.insert(1);
    mgr.assign_feature(body, fid);

    mgr.consume_body(body);

    // The body should no longer appear in the list.
    assert!(!mgr.bodies().contains(&body));
    // body_features on a consumed body returns empty (not found).
    assert!(mgr.body_features(body).is_empty());
}

#[test]
fn clear_breps_removes_snapshots() {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let mut mgr = BodyManager::new();
    let body = mgr.create_body();

    let tol = Tolerance::default();
    mgr.set_brep(body, BRep::new(tol));
    assert!(mgr.get_brep(body).is_some());

    mgr.clear_breps();
    assert!(mgr.get_brep(body).is_none());
}

#[test]
fn clear_features_empties_lists_but_keeps_bodies() {
    let mut mgr = BodyManager::new();
    let body = mgr.create_body();

    let mut sm: SlotMap<FeatureId, u32> = SlotMap::with_key();
    let fid = sm.insert(1);
    mgr.assign_feature(body, fid);

    mgr.clear_features();

    // Body still exists.
    assert!(mgr.bodies().contains(&body));
    // Feature list is empty.
    assert!(mgr.body_features(body).is_empty());
}

#[test]
fn bodies_lists_all_created_bodies() {
    let mut mgr = BodyManager::new();
    let b0 = mgr.create_body();
    let b1 = mgr.create_body();

    let all = mgr.bodies();
    assert_eq!(all.len(), 2);
    assert!(all.contains(&b0));
    assert!(all.contains(&b1));
}

// ---------------------------------------------------------------------------
// FeatureTree integration tests
// ---------------------------------------------------------------------------

#[test]
fn two_independent_sketches_create_two_bodies() {
    use crate::FeatureTree;
    use crate::feature::{Feature, SketchPlaneDesc};

    let mut tree = FeatureTree::new();

    let _s1 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        "Sketch 1",
    );
    let _s2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: String::new(),
        },
        "Sketch 2",
    );

    let bodies = tree.bodies();
    assert_eq!(
        bodies.len(),
        2,
        "two independent sketches should create two bodies"
    );
    assert!(bodies.contains(&BodyId(0)));
    assert!(bodies.contains(&BodyId(1)));
}

#[test]
fn extrude_inherits_sketch_body() {
    use crate::FeatureTree;
    use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};
    use crate::unit::Unit;
    use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};

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

    let mut tree = FeatureTree::new();

    // Sketch gets OperationId(1).
    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        "Sketch 1",
    );

    let dist_param = tree
        .parameters_mut()
        .add_global("_extrude1_distance", 10.0, Unit::Meters)
        .expect("add param");

    // Extrude references a face from OperationId(1) → depends on sketch.
    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist_param,
        },
        "Extrude 1",
    );

    let sketch_body = tree.feature(sketch_id).expect("sketch exists").body;
    let extrude_body = tree.feature(extrude_id).expect("extrude exists").body;
    assert_eq!(
        sketch_body, extrude_body,
        "extrude should inherit the same body as its sketch dependency"
    );
}

#[test]
fn body_features_returns_correct_list() {
    use crate::FeatureTree;
    use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};
    use crate::unit::Unit;
    use popkern_core::naming::{EntityName, EntityRole, NameId, OperationId, RoleCategory};

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

    let mut tree = FeatureTree::new();

    let sketch_id = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        "Sketch 1",
    );

    let dist_param = tree
        .parameters_mut()
        .add_global("_extrude1_distance", 5.0, Unit::Meters)
        .expect("add param");

    let extrude_id = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: dist_param,
        },
        "Extrude 1",
    );

    let sketch_body = tree.feature(sketch_id).expect("sketch exists").body;
    let features = tree.body_features(sketch_body);

    assert!(
        features.contains(&sketch_id),
        "body_features should contain the sketch"
    );
    assert!(
        features.contains(&extrude_id),
        "body_features should contain the extrude"
    );
    assert_eq!(features.len(), 2);
}

#[test]
fn body_brep_returns_none_before_replay() {
    use crate::FeatureTree;
    use crate::feature::{Feature, SketchPlaneDesc};

    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        "Sketch 1",
    );

    // No replay has happened — body_brep should be None.
    assert!(tree.body_brep(BodyId(0)).is_none());
}
