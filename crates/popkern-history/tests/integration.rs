//! Integration tests for the `popkern-history` crate.
//!
//! Exercises the full pipeline through the public API: sketch, extrude,
//! fillet, timeline scrubbing, save/load, undo/redo, delete cascade,
//! suppress/unsuppress, and broken-reference handling.

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

fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

/// Generated face at `op_id` index 0 (sketch face).
fn sketch_face_name(op_id: u64) -> EntityName {
    entity_name(op_id, RoleCategory::Generated, 0)
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
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
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

/// First named edge in the BRep (for fillet/chamfer targets).
fn find_named_edge(tree: &FeatureTree) -> EntityName {
    let brep = tree.current_brep().expect("brep should exist");
    for (_h, edge) in &brep.edges {
        if let Some(nid) = edge.name
            && let Some(name) = brep.name_pool.get(nid)
        {
            return name.clone();
        }
    }
    panic!("no named edge found in BRep");
}

// -- 1. Sketch -> extrude -> fillet pipeline ---------------------------------

#[test]
fn sketch_extrude_fillet_pipeline() {
    let (mut tree, _sk, _ex) = build_sketch_extrude(1.0);
    tree.replay(Tolerance::default()).expect("initial replay");
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);

    let edge_name = find_named_edge(&tree);
    let radius_param = add_param(&mut tree, "_fillet1_radius", 0.1, Unit::Meters);
    let fillet_id = tree.add_feature(
        Feature::Fillet {
            edges: vec![edge_name],
            radius: radius_param,
            recipe: Default::default(),
        },
        "Fillet1",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");
    let sk = tree.timeline()[0];
    let ex = tree.timeline()[1];
    assert!(matches!(
        result.feature_results.get(&sk),
        Some(FeatureResult::Success { .. })
    ));
    assert!(matches!(
        result.feature_results.get(&ex),
        Some(FeatureResult::Success { .. })
    ));
    assert!(
        result.feature_results.contains_key(&fillet_id),
        "fillet evaluated"
    );
    assert!(tree.current_brep().unwrap().faces.len() >= 6);
}

// -- 2. Modify extrude distance -> rebuild -----------------------------------

#[test]
fn modify_extrude_distance_and_rebuild() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(1.0);
    tree.replay(Tolerance::default()).expect("replay");
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);

    let new_dist_param = add_param(&mut tree, "_extrude1_distance_v2", 2.0, Unit::Meters);

    tree.modify_feature(
        extrude_id,
        Feature::Extrude {
            face: sketch_face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind {
                distance: new_dist_param,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
    );
    tree.replay_from(extrude_id, Tolerance::default())
        .expect("replay_from");

    let brep = tree.current_brep().unwrap();
    assert_eq!(brep.faces.len(), 6, "still 6 faces");
    let max_z = brep
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (max_z - 2.0).abs() < 0.01,
        "height should be 2.0, got {max_z}"
    );
}

// -- 3. Timeline scrubbing with brep_at --------------------------------------

#[test]
fn timeline_scrubbing() {
    let (mut tree, _sk, _ex) = build_sketch_extrude(1.0);
    let tol = Tolerance::default();

    // After sketch: flat face.
    let brep1 = tree.brep_at(1, tol).expect("brep_at(1)");
    assert!(!brep1.faces.is_empty() && brep1.faces.len() < 6);

    // After extrude: box.
    let brep2 = tree.brep_at(2, tol).expect("brep_at(2)");
    assert_eq!(brep2.faces.len(), 6);
}

// -- 4. Rollback + insert ----------------------------------------------------

#[test]
fn rollback_and_insert_feature() {
    let (mut tree, sk, ex) = build_sketch_extrude(1.0);
    let sk2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch2",
    );
    assert_eq!(tree.timeline().len(), 3);

    tree.rollback_to(2);
    assert_eq!(tree.rollback_position(), Some(2));

    let ins = tree
        .insert_feature_at_rollback(
            Feature::Sketch {
                plane: SketchPlaneDesc::YZ,
                sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
            },
            "InsertedSketch",
        )
        .expect("insert");

    tree.rollback_to_end();
    let tl = tree.timeline();
    assert_eq!(tl.len(), 4);
    assert_eq!(tl[0], sk);
    assert_eq!(tl[1], ex);
    assert_eq!(tl[2], ins);
    assert_eq!(tl[3], sk2);

    tree.replay(Tolerance::default()).expect("replay");
    assert!(tree.feature(ins).unwrap().last_result.is_some());
}

// -- 5. Delete cascade -------------------------------------------------------

#[test]
fn delete_extrude_cascades() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(1.0);
    tree.replay(Tolerance::default()).expect("replay");

    let edge_name = find_named_edge(&tree);
    let radius_param = add_param(&mut tree, "_fillet1_radius", 0.1, Unit::Meters);
    let _fillet_id = tree.add_feature(
        Feature::Fillet {
            edges: vec![edge_name],
            radius: radius_param,
            recipe: Default::default(),
        },
        "Fillet1",
    );
    assert_eq!(tree.timeline().len(), 3);

    let del = tree
        .delete_feature(extrude_id, DeleteBehavior::Cascade)
        .expect("delete");
    assert!(del.deleted.contains(&extrude_id));

    tree.replay(Tolerance::default())
        .expect("replay after delete");
    assert!(tree.current_brep().unwrap().faces.len() < 6);
}

// -- 6. Suppress / unsuppress ------------------------------------------------

#[test]
fn suppress_and_unsuppress_extrude() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(1.0);
    tree.replay(Tolerance::default()).expect("replay");
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);

    tree.suppress(extrude_id);
    tree.replay(Tolerance::default())
        .expect("replay suppressed");
    assert!(tree.current_brep().unwrap().faces.len() < 6);

    tree.unsuppress(extrude_id);
    tree.replay(Tolerance::default())
        .expect("replay unsuppressed");
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);
}

// -- 7. Save / load round-trip -----------------------------------------------

#[test]
fn save_load_round_trip() {
    let (mut tree, _sk, _ex) = build_sketch_extrude(1.0);
    tree.replay(Tolerance::default()).expect("replay");
    let orig_faces = tree.current_brep().unwrap().faces.len();
    let orig_len = tree.timeline().len();

    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    assert!(loaded.current_brep().is_none(), "no brep before replay");

    loaded.replay(Tolerance::default()).expect("loaded replay");
    assert_eq!(loaded.current_brep().unwrap().faces.len(), orig_faces);
    assert_eq!(loaded.timeline().len(), orig_len);
}

// -- 8. Broken reference -> Failed -------------------------------------------

#[test]
fn broken_reference_produces_failed() {
    let (mut tree, _sk, extrude_id) = build_sketch_extrude(1.0);
    let bogus = entity_name(99, RoleCategory::Lateral, 0);
    let radius_param = add_param(&mut tree, "_broken_fillet_radius", 0.1, Unit::Meters);
    let fillet_id = tree.add_feature(
        Feature::Fillet {
            edges: vec![bogus],
            radius: radius_param,
            recipe: Default::default(),
        },
        "BrokenFillet",
    );

    let result = tree.replay(Tolerance::default()).expect("replay");

    let sk = tree.timeline()[0];
    assert!(matches!(
        result.feature_results.get(&sk),
        Some(FeatureResult::Success { .. })
    ));
    assert!(matches!(
        result.feature_results.get(&extrude_id),
        Some(FeatureResult::Success { .. })
    ));
    match result.feature_results.get(&fillet_id) {
        Some(FeatureResult::Failed { error }) => assert!(!error.is_empty()),
        other => panic!("expected Failed, got: {other:?}"),
    }
    assert_eq!(
        tree.current_brep().unwrap().faces.len(),
        6,
        "last good state"
    );
}

// -- 9. Tree undo / redo -----------------------------------------------------

#[test]
fn tree_undo_redo() {
    let (mut tree, _sk, _ex) = build_sketch_extrude(1.0);
    assert_eq!(tree.timeline().len(), 2);

    tree.undo();
    assert_eq!(tree.timeline().len(), 1);
    assert!(tree.can_redo());

    tree.replay(Tolerance::default())
        .expect("replay after undo");
    assert!(tree.current_brep().unwrap().faces.len() < 6);

    tree.redo();
    assert_eq!(tree.timeline().len(), 2);
    tree.replay(Tolerance::default())
        .expect("replay after redo");
    assert_eq!(tree.current_brep().unwrap().faces.len(), 6);
}

// -- 10. Dependency and impact queries ---------------------------------------

#[test]
fn dependency_and_impact_queries() {
    let (mut tree, sketch_id, extrude_id) = build_sketch_extrude(1.0);
    tree.replay(Tolerance::default()).expect("replay");

    let edge_name = find_named_edge(&tree);
    let radius_param = add_param(&mut tree, "_fillet1_radius", 0.1, Unit::Meters);
    let fillet_id = tree.add_feature(
        Feature::Fillet {
            edges: vec![edge_name],
            radius: radius_param,
            recipe: Default::default(),
        },
        "Fillet1",
    );
    tree.replay(Tolerance::default()).expect("replay");

    assert!(tree.dependencies(extrude_id).contains(&sketch_id));
    assert!(tree.dependents(sketch_id).contains(&extrude_id));
    assert!(!tree.dependencies(fillet_id).is_empty());
    assert!(
        tree.impact_analysis(sketch_id)
            .will_rebuild
            .contains(&extrude_id)
    );
}

// -- 11. Multiple sketches create separate bodies ----------------------------

#[test]
fn multiple_sketches_create_separate_bodies() {
    let mut tree = FeatureTree::new();
    let sk1 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let sk2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch2",
    );
    assert_ne!(
        tree.feature(sk1).unwrap().body,
        tree.feature(sk2).unwrap().body
    );

    tree.replay(Tolerance::default()).expect("replay");
    assert!(!tree.current_brep().unwrap().faces.is_empty());
}

// -- 12. Save/load preserves rollback index ----------------------------------

#[test]
fn save_load_preserves_rollback_index() {
    let (mut tree, _sk, _ex) = build_sketch_extrude(1.0);
    tree.rollback_to(1);
    assert_eq!(tree.rollback_position(), Some(1));

    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");
    assert_eq!(loaded.rollback_position(), Some(1));
}

// Multi-body and Boolean integration tests are in multi_body.rs.
