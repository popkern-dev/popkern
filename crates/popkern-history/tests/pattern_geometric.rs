//! Geometric integration tests for the feature-pattern pipeline.
//!
//! Unlike `pattern_integration.rs` and `pattern_integration_2.rs` (which test
//! serialization round-trips), these tests exercise the full BRep pipeline:
//! Sketch -> Extrude -> Pattern -> replay -> face count checks.

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

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

/// Solved 0.3x0.3 small square sketch (for boss features on top faces).
fn small_square_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.2, 0.2);
    let p1 = sk.add_point(0.5, 0.2);
    let p2 = sk.add_point(0.5, 0.5);
    let p3 = sk.add_point(0.2, 0.5);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    sk.solve().expect("small square should solve");
    serde_json::to_string(&sk).expect("small square should serialize")
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

/// Find the first face in the BRep whose name has the given operation and
/// role category. Returns the full EntityName (including ancestors).
fn find_face_name(
    tree: &FeatureTree,
    body: BodyId,
    op_id: u64,
    category: RoleCategory,
) -> EntityName {
    let brep = tree.body_brep(body).expect("body must have BRep");
    for (_h, face) in &brep.faces {
        if let Some(nid) = face.name
            && let Some(name) = brep.name_pool.get(nid)
            && name.operation == OperationId(op_id)
            && name.role.category == category
        {
            return name.clone();
        }
    }
    panic!("no face found with operation={op_id}, category={category:?} in body {body:?}");
}

/// Build a base box: Sketch(XY) + Extrude(NewBody), replay, and return the
/// tree along with identifiers. Returns (tree, sketch_id, extrude_id).
fn build_base_box(distance: f64) -> (FeatureTree, FeatureId, FeatureId) {
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
    tree.replay(Tolerance::default()).expect("base box replay");
    (tree, sketch_id, extrude_id)
}

/// Build a base box + Boss on top face. Returns (tree, base_body, sk2, ex2).
///
/// After this call, the tree has been replayed and the Boss has been
/// auto-merged into the base body.
fn build_box_with_boss() -> (FeatureTree, BodyId, FeatureId, FeatureId) {
    let (mut tree, sk1, _ex1) = build_base_box(1.0);
    let base_body = tree.feature(sk1).expect("sk1").body;

    // Find the cap face (top face) of the extrude dynamically.
    // Extrude1 has OperationId(2); its cap face includes ancestors.
    let cap_face = find_face_name(&tree, base_body, 2, RoleCategory::Cap);

    // Boss feature: Sketch on top face + Extrude(Boss).
    let sk2 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::OnFace(cap_face),
            sketch_data: SketchPayload::Constraint(small_square_sketch_json()),
        },
        "Sketch2",
    );
    tree.replay(Tolerance::default()).expect("replay sk2");

    // Sketch2 is OperationId(3). Its sketch face name is Generated(3, 0).
    let dist_ex2 = add_param(&mut tree, "_extrude2_dist", 0.3, Unit::Meters);
    let ex2 = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(3),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist_ex2 },
            mode: OperationMode::Boss,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude2",
    );
    tree.replay(Tolerance::default()).expect("replay boss");

    (tree, base_body, sk2, ex2)
}

// ---------------------------------------------------------------------------
// Test 1: Baseline -- sketch + extrude produces a 6-face box
// ---------------------------------------------------------------------------

#[test]
fn sketch_extrude_produces_valid_box() {
    let (tree, _sk, _ex) = build_base_box(1.0);

    let brep = tree.current_brep().expect("should have BRep");
    assert_eq!(brep.faces.len(), 6, "box should have 6 faces");
}

// ---------------------------------------------------------------------------
// Test 2: Linear pattern with Boss seed produces multiple instances
// ---------------------------------------------------------------------------

#[test]
fn linear_pattern_boss_seed_produces_multiple_instances() {
    let (mut tree, base_body, sk2, ex2) = build_box_with_boss();

    // Note: OnFace sketch + Boss extrude share BodyId(0) with the base box.
    // The boss extrude modifies the same BRep rather than creating a separate
    // tool body, so we measure faces before adding the pattern.
    let faces_before_pattern = tree
        .body_brep(base_body)
        .expect("base body before pattern")
        .faces
        .len();

    // LinearPattern on the boss seed features: count=3 along X.
    let count_param = add_param(&mut tree, "_pat_count", 3.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_pat_total", 2.0, Unit::Meters);

    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![sk2, ex2],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPattern1",
    );

    let result = tree.replay(Tolerance::default()).expect("replay pattern");

    // The pattern feature should succeed.
    match result.feature_results.get(&pattern_id) {
        Some(FeatureResult::Success { .. }) => {}
        other => panic!("expected LinearPattern Success, got: {other:?}"),
    }

    // After pattern with 3 instances, face count should increase (instances
    // 1..N replay the seed chain and boolean the result into the body).
    let final_brep = tree.body_brep(base_body).expect("base body after pattern");
    assert!(
        final_brep.faces.len() > faces_before_pattern,
        "3-instance pattern should add more faces (got {} vs {faces_before_pattern} before)",
        final_brep.faces.len()
    );
}

// ---------------------------------------------------------------------------
// Test 3: Pattern count=1 is just the original (no extra instances)
// ---------------------------------------------------------------------------

#[test]
fn pattern_count_1_is_just_original() {
    let (mut tree, base_body, sk2, ex2) = build_box_with_boss();

    let faces_after_boss = tree.body_brep(base_body).expect("base body").faces.len();

    // Pattern with count=1: only instance 0 (the original), no transforms.
    let count_param = add_param(&mut tree, "_pc1_count", 1.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_pc1_total", 0.0, Unit::Meters);

    let pattern_id = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![sk2, ex2],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPattern1",
    );

    let result = tree.replay(Tolerance::default()).expect("replay pattern");
    match result.feature_results.get(&pattern_id) {
        Some(FeatureResult::Success { .. }) => {}
        other => panic!("expected Success for count=1 pattern, got: {other:?}"),
    }

    // With count=1, the pattern should not add any faces beyond the boss.
    let final_faces = tree
        .body_brep(base_body)
        .expect("base body after pattern")
        .faces
        .len();
    assert_eq!(
        final_faces, faces_after_boss,
        "count=1 pattern should not change face count (got {final_faces}, expected {faces_after_boss})"
    );
}

// ---------------------------------------------------------------------------
// Test 4: Pattern with skip instance reduces the effective count
// ---------------------------------------------------------------------------

#[test]
fn pattern_skip_instance_reduces_count() {
    let (mut tree, base_body, sk2, ex2) = build_box_with_boss();

    // Pattern A: count=4, no skips.
    let count_a = add_param(&mut tree, "_patA_count", 4.0, Unit::Ratio);
    let total_a = add_param(&mut tree, "_patA_total", 3.0, Unit::Meters);

    let pat_a = tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![sk2, ex2],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_a,
                total_distance: total_a,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPatternA",
    );

    let result_a = tree.replay(Tolerance::default()).expect("replay patA");
    match result_a.feature_results.get(&pat_a) {
        Some(FeatureResult::Success { .. }) => {}
        other => panic!("expected PatternA Success, got: {other:?}"),
    }
    let faces_4_instances = tree
        .body_brep(base_body)
        .expect("base body after 4-instance pattern")
        .faces
        .len();

    // Modify to skip instance 2.
    let mut overrides: std::collections::HashMap<usize, pattern::InstanceOverride> =
        Default::default();
    overrides.insert(
        2,
        pattern::InstanceOverride {
            suppressed: true,
            param_overrides: Default::default(),
        },
    );

    tree.modify_feature(
        pat_a,
        Feature::LinearPattern {
            seed_features: vec![sk2, ex2],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_a,
                total_distance: total_a,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: overrides,
        },
    );

    let result_b = tree.replay(Tolerance::default()).expect("replay with skip");
    match result_b.feature_results.get(&pat_a) {
        Some(FeatureResult::Success { .. }) => {}
        other => panic!("expected modified pattern Success, got: {other:?}"),
    }
    let faces_3_instances = tree
        .body_brep(base_body)
        .expect("base body after skip")
        .faces
        .len();

    // Skipping one instance should reduce the face count.
    assert!(
        faces_3_instances < faces_4_instances,
        "skipping instance 2 should reduce face count ({faces_3_instances} < {faces_4_instances})"
    );
}

// ---------------------------------------------------------------------------
// Test 5: Mirror pattern produces symmetric geometry
// ---------------------------------------------------------------------------

#[test]
fn mirror_pattern_produces_symmetric_geometry() {
    let (mut tree, base_body, sk2, ex2) = build_box_with_boss();

    let faces_before_mirror = tree.body_brep(base_body).expect("base body").faces.len();

    // Mirror about YZ plane (X = 0, normal = +X).
    let mirror_id = tree.add_feature(
        Feature::Mirror {
            seed_features: vec![sk2, ex2],
            plane: pattern::MirrorPlaneDesc::YZ,
            merge_result: true,
        },
        "Mirror1",
    );

    let result = tree.replay(Tolerance::default()).expect("replay mirror");
    match result.feature_results.get(&mirror_id) {
        Some(FeatureResult::Success { .. }) => {}
        other => panic!("expected Mirror Success, got: {other:?}"),
    }

    let mirror_brep = tree.body_brep(base_body).expect("base body after mirror");

    // Mirror should add faces from the reflected boss.
    assert!(
        mirror_brep.faces.len() > faces_before_mirror,
        "mirror should increase face count ({} should be > {faces_before_mirror})",
        mirror_brep.faces.len()
    );

    // Check vertex symmetry about the YZ plane (X coordinates should have
    // pairs: for every vertex at +x there should be one near -x).
    let xs: Vec<f64> = mirror_brep
        .vertices
        .iter()
        .map(|(_, v)| v.position.x)
        .collect();
    // Only the SEEDS mirror: every boss vertex (x in [0.2, 0.5]) has a
    // counterpart at -x; the box itself stays where it is. (Before
    // PK-387 the seeds lived in the host body, so the pattern's tool
    // snapshot was the whole box + boss and the "mirror" cloned the
    // host — every vertex had a counterpart by accident.)
    let boss_xs: Vec<f64> = xs
        .iter()
        .copied()
        .filter(|x| (0.2..=0.5).contains(x))
        .collect();
    assert!(!boss_xs.is_empty(), "boss vertices present");
    for &x in &boss_xs {
        let has_mirror = xs.iter().any(|&other| (other + x).abs() < 0.1);
        assert!(
            has_mirror,
            "boss vertex at x={x} should have a mirrored counterpart near x={}",
            -x
        );
    }
    assert!(
        xs.iter().all(|&x| (-0.5 - 1e-9..=1.0 + 1e-9).contains(&x)),
        "the box is not mirrored: x range {:?}",
        (
            xs.iter().cloned().fold(f64::MAX, f64::min),
            xs.iter().cloned().fold(f64::MIN, f64::max)
        )
    );
}

// ---------------------------------------------------------------------------
// Test 6: Pattern save/load/replay produces matching face counts
// ---------------------------------------------------------------------------

#[test]
fn pattern_save_load_replay_matches() {
    let (mut tree, base_body, sk2, ex2) = build_box_with_boss();

    // LinearPattern with count=2.
    let count_param = add_param(&mut tree, "_pat_count", 2.0, Unit::Ratio);
    let total_param = add_param(&mut tree, "_pat_total", 1.5, Unit::Meters);

    tree.add_feature(
        Feature::LinearPattern {
            seed_features: vec![sk2, ex2],
            primary: pattern::PatternDirection {
                direction: Vec3Desc::Literal([1.0, 0.0, 0.0]),
                count: count_param,
                total_distance: total_param,
            },
            secondary: None,
            geometry_pattern: false,
            instance_overrides: Default::default(),
        },
        "LinPattern1",
    );

    tree.replay(Tolerance::default()).expect("replay original");
    let orig_faces = tree
        .body_brep(base_body)
        .expect("original BRep")
        .faces
        .len();
    let orig_timeline_len = tree.timeline().len();

    // Save and load round-trip.
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    assert!(
        loaded.current_brep().is_none(),
        "loaded tree should have no BRep before replay"
    );

    loaded.replay(Tolerance::default()).expect("replay loaded");

    // The loaded tree should have the same timeline length.
    assert_eq!(
        loaded.timeline().len(),
        orig_timeline_len,
        "timeline length should match after round-trip"
    );

    // The face count after replaying the loaded tree should match.
    let loaded_brep = loaded.current_brep().expect("loaded BRep");
    assert_eq!(
        loaded_brep.faces.len(),
        orig_faces,
        "face count should match after save/load/replay ({} vs {orig_faces})",
        loaded_brep.faces.len()
    );
}
