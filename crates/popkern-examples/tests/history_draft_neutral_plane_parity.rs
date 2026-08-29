//! PK-363 acceptance: a `Feature::Draft` whose neutral plane is a
//! literal `{origin, normal}` — the recipe `Part::draft` derives
//! implicitly (perpendicular to the pull direction, through the
//! part's lowest point along it) — replays to the same body as the
//! `popkern-part` call on a shared fixture, survives `.pkn`
//! save/load, and takes the in-place angle-tweak route. The face arm
//! and the literal arm agree when they spell the same plane.
//!
//! Parity is asserted through numeric gates — face/edge/vertex
//! census, sorted vertex positions, body volume — never visual
//! inspection.

use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::*;
use popkern_part::prelude::*;
use popkern_sketch::{Constraint, Sketch};

const DRAFT_ANGLE: f64 = 0.1;

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
    for p in [p0, p1, p2, p3] {
        sk.add_constraint(Constraint::Fixed(p));
    }
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

fn name(op: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

fn add_param(tree: &mut FeatureTree, label: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(label, value, unit)
        .expect("add param")
}

/// Unit cube via sketch + extrude (features 1 and 2), replayed.
/// Returns the tree, the cube, and the four lateral face names.
fn cube_tree() -> (FeatureTree, BRep, Vec<EntityName>) {
    let mut tree = FeatureTree::new();
    tree.set_cache_strategy(CacheStrategy::EveryStep);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let dist = add_param(&mut tree, "_extrude_distance", 1.0, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );
    tree.replay(Tolerance::default()).expect("cube replays");
    let brep = tree.current_brep().expect("cube body").clone();
    assert_eq!(brep.faces.len(), 6);
    // Lateral faces by stable name, read off the body (their names
    // carry ancestry, so they are looked up rather than spelled).
    let mut lateral = Vec::new();
    for (h, f) in &brep.faces {
        let pts = brep.face_outer_loop_positions(h);
        let planar_z = pts.iter().all(|p| (p.z - pts[0].z).abs() < 1e-12);
        if !planar_z {
            let nid = f.name.expect("cube faces are named");
            lateral.push(brep.name_pool.get(nid).expect("name pool").clone());
        }
    }
    lateral.sort_by_key(|n| n.role.index);
    assert_eq!(lateral.len(), 4, "a cube has four side walls");
    (tree, brep, lateral)
}

fn assert_replay_clean(result: &ReplayResult) {
    let failures: Vec<_> = result
        .feature_results
        .iter()
        .filter(|(_, r)| !matches!(r, FeatureResult::Success { .. }))
        .collect();
    assert!(failures.is_empty(), "replay failures: {failures:?}");
}

fn sorted_vertices(brep: &BRep) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    pts
}

/// Assert two bodies are the same geometry to within `tol`.
fn assert_same_body(label: &str, a: &BRep, b: &BRep, tol: f64) {
    assert_eq!(a.faces.len(), b.faces.len(), "{label}: face count");
    assert_eq!(a.edges.len(), b.edges.len(), "{label}: edge count");
    assert_eq!(a.vertices.len(), b.vertices.len(), "{label}: vertex count");
    for (p, q) in sorted_vertices(a).iter().zip(sorted_vertices(b).iter()) {
        for k in 0..3 {
            assert!((p[k] - q[k]).abs() < tol, "{label}: vertex {p:?} vs {q:?}");
        }
    }
    let t = Tolerance::default();
    let va = popkern_measure::body_volume(a, &t).expect("volume a");
    let vb = popkern_measure::body_volume(b, &t).expect("volume b");
    assert!(
        (va - vb).abs() < tol * va.abs().max(1.0),
        "{label}: volume {va} vs {vb}"
    );
}

/// Add a draft feature on the lateral faces, save/load, replay;
/// returns the loaded tree, the feature id, and the replayed body.
fn replay_draft(
    mut tree: FeatureTree,
    faces: Vec<EntityName>,
    angle: ParameterId,
    neutral_plane: NeutralPlaneDesc,
) -> (FeatureTree, FeatureId, BRep) {
    tree.add_feature(
        Feature::Draft {
            faces,
            angle,
            neutral_plane,
        },
        "Draft1",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    loaded.set_cache_strategy(CacheStrategy::EveryStep);
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let draft_id = *loaded.timeline().last().expect("draft feature");
    let body = loaded.current_brep().expect("body").clone();
    body.validate().expect("replayed drafted body validates");
    (loaded, draft_id, body)
}

/// `Part::draft` on the same cube: pull +Z, the four side walls.
fn part_draft(cube: &BRep, angle: f64) -> Part {
    Part::from_brep(cube.clone(), Tolerance::default(), "cube")
        .expect("part")
        .draft(angle, Dir::Z, |f: FaceQuery| {
            f.filter(|info| info.normal.is_some_and(|n| n.z.abs() < 0.5))
        })
        .expect("part draft")
}

#[test]
fn literal_plane_draft_replays_to_part_parity() {
    let (mut tree, cube, lateral) = cube_tree();
    let angle = add_param(&mut tree, "draft", DRAFT_ANGLE, Unit::Radians);
    // Part::draft's implicit plane: perpendicular to +Z through the
    // lowest vertex along Z (z = 0), spelled literally.
    let (_, _, replayed) = replay_draft(
        tree,
        lateral,
        angle,
        NeutralPlaneDesc::Literal {
            origin: [0.3, -0.2, 0.0],
            normal: [0.0, 0.0, 3.0],
        },
    );
    let part = part_draft(&cube, DRAFT_ANGLE);
    assert_eq!(replayed.faces.len(), 6, "drafted cube keeps six faces");
    let base = popkern_measure::body_volume(&cube, &Tolerance::default()).expect("v");
    let drafted = popkern_measure::body_volume(&replayed, &Tolerance::default()).expect("v");
    assert!(
        (drafted - base).abs() > 1e-3,
        "the draft must change the body: {base} vs {drafted}"
    );
    assert_same_body(
        "literal-plane draft vs Part::draft",
        &replayed,
        part.brep(),
        1e-9,
    );
}

/// The face arm and the literal arm are the same plane spelled two
/// ways: the top cap (`z = 1`, outward +Z) versus `{[0,0,1], +Z}`.
#[test]
fn face_arm_and_literal_arm_agree_on_the_same_plane() {
    let (mut tree, _cube, lateral) = cube_tree();
    let angle = add_param(&mut tree, "draft", DRAFT_ANGLE, Unit::Radians);
    let cube = tree.current_brep().expect("cube");
    let top = cube
        .faces
        .iter()
        .find_map(|(h, f)| {
            let pts = cube.face_outer_loop_positions(h);
            let is_top = !pts.is_empty() && pts.iter().all(|p| (p.z - 1.0).abs() < 1e-12);
            is_top.then_some(cube.name_pool.get(f.name?)?.clone())
        })
        .expect("top cap is named");
    let (_, _, via_face) = replay_draft(
        tree.clone(),
        lateral.clone(),
        angle,
        NeutralPlaneDesc::Face(top),
    );
    let (_, _, via_literal) = replay_draft(
        tree,
        lateral,
        angle,
        NeutralPlaneDesc::Literal {
            origin: [0.0, 0.0, 1.0],
            normal: [0.0, 0.0, 1.0],
        },
    );
    assert_same_body("face arm vs literal arm", &via_face, &via_literal, 1e-9);
}

/// An angle tweak on a literal-plane draft takes the in-place route
/// and lands on the same body a full replay at the new angle gives —
/// which is also `Part::draft` at the new angle.
#[test]
fn literal_plane_angle_tweak_matches_full_replay_and_part() {
    let (mut tree, cube, lateral) = cube_tree();
    let angle = add_param(&mut tree, "draft", DRAFT_ANGLE, Unit::Radians);
    let plane = NeutralPlaneDesc::Literal {
        origin: [0.0, 0.0, 0.0],
        normal: [0.0, 0.0, 1.0],
    };
    let (mut loaded, draft_id, _) = replay_draft(tree, lateral.clone(), angle, plane.clone());
    let new_angle = 0.17;
    // The in-place tweak recipe: a fresh parameter carrying the new
    // value, swapped in through `modify_feature` (same faces, same
    // plane), which refreshes the cached state via
    // `modify_draft_angle` on the literal plane.
    let angle_new = add_param(&mut loaded, "draft2", new_angle, Unit::Radians);
    loaded.modify_feature(
        draft_id,
        Feature::Draft {
            faces: lateral.clone(),
            angle: angle_new,
            neutral_plane: plane.clone(),
        },
    );
    let result = loaded
        .replay(Tolerance::default())
        .expect("replay after tweak");
    assert_replay_clean(&result);
    let tweaked = loaded.current_brep().expect("tweaked body").clone();

    let (mut fresh, _, _) = cube_tree();
    let fresh_angle = add_param(&mut fresh, "draft", new_angle, Unit::Radians);
    let (_, _, full) = replay_draft(fresh, lateral, fresh_angle, plane);
    assert_same_body("in-place tweak vs full replay", &tweaked, &full, 1e-9);
    let part = part_draft(&cube, new_angle);
    assert_same_body("in-place tweak vs Part::draft", &tweaked, part.brep(), 1e-9);
}

/// A degenerate literal normal declines transactionally: the feature
/// reports failure and the body is left as the undrafted cube.
#[test]
fn degenerate_literal_normal_declines_transactionally() {
    let (mut tree, cube, lateral) = cube_tree();
    let angle = add_param(&mut tree, "draft", DRAFT_ANGLE, Unit::Radians);
    tree.add_feature(
        Feature::Draft {
            faces: lateral,
            angle,
            neutral_plane: NeutralPlaneDesc::Literal {
                origin: [0.0; 3],
                normal: [0.0; 3],
            },
        },
        "Draft1",
    );
    let result = tree.replay(Tolerance::default()).expect("replay returns");
    let failed = result
        .feature_results
        .values()
        .any(|r| matches!(r, FeatureResult::Failed { error, .. } if error.contains("degenerate")));
    assert!(
        failed,
        "degenerate normal must be reported: {:?}",
        result.feature_results
    );
    assert_same_body(
        "body untouched after decline",
        tree.current_brep().expect("body"),
        &cube,
        1e-12,
    );
}
