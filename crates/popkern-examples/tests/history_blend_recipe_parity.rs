//! PK-361 acceptance: fillet/chamfer recipes stored in a feature tree
//! replay to the same body as the equivalent `popkern-part` call on a
//! shared fixture, and survive a `.pkn` save/load with the options
//! honoured.
//!
//! Parity is asserted through numeric gates — face/edge/vertex census,
//! sorted vertex positions, body volume, and the measured continuity
//! evidence over the blend faces — never through visual inspection.

use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_history::*;
use popkern_part::dressup::Continuity;
use popkern_part::prelude::*;
use popkern_sketch::{Constraint, Sketch};

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

fn add_param(tree: &mut FeatureTree, name: &str, value: f64) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, Unit::Meters)
        .expect("add param")
}

/// Unit cube via sketch + extrude, replayed. Returns the tree and the
/// pre-fillet body.
fn cube_tree() -> (FeatureTree, BRep) {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let dist = add_param(&mut tree, "_extrude_distance", 1.0);
    tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(1),
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
    (tree, brep)
}

/// Stable names of the cube's edges whose tangent is along ±Z.
fn vertical_edge_names(brep: &BRep) -> Vec<EntityName> {
    let mut names = Vec::new();
    for (handle, edge) in &brep.edges {
        let tangent = popkern_measure::edge_tangent(brep, handle).expect("edge tangent");
        if tangent.z.abs() > 0.99 {
            let nid = edge.name.expect("cube edges are named");
            names.push(brep.name_pool.get(nid).expect("name pool").clone());
        }
    }
    assert_eq!(names.len(), 4, "a cube has four vertical edges");
    names
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

fn new_faces(pre: &BRep, post: &BRep) -> Vec<FaceHandle> {
    let before: std::collections::BTreeSet<_> = pre.faces.keys().collect();
    post.faces.keys().filter(|f| !before.contains(f)).collect()
}

fn evidence(
    pre: &BRep,
    post: &BRep,
    requested: Continuity,
) -> popkern_blend::BlendContinuityEvidence {
    popkern_blend::measure_blend_continuity_evidence(
        post,
        &new_faces(pre, post),
        requested,
        &popkern_blend::ContinuityGates::commercial(),
        33,
    )
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

/// Replay a fillet recipe on the cube's vertical edges, through a
/// `.pkn` save/load so the options are proven to persist; returns the
/// pre-fillet body and the replayed result.
fn replay_fillet(recipe: FilletRecipeDesc, radius: f64) -> (BRep, BRep) {
    let (mut tree, cube) = cube_tree();
    let edges = vertical_edge_names(&cube);
    let r = add_param(&mut tree, "_r", radius);
    tree.add_feature(
        Feature::Fillet {
            edges,
            radius: r,
            recipe,
        },
        "Fillet1",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let body = loaded.current_brep().expect("body").clone();
    body.validate().expect("replayed body validates");
    (cube, body)
}

fn part_fillet(
    cube: &BRep,
    recipe: impl FnOnce(popkern_part::dressup::FilletRecipe) -> popkern_part::dressup::FilletRecipe,
) -> Part {
    let part = Part::from_brep(cube.clone(), Tolerance::default(), "cube").expect("part");
    part.fillet_with(recipe, |e| e.straight().direction(Dir::Z))
        .expect("part fillet serves")
}

/// G2 was unreachable from a saved model before PK-361 (replay
/// hardcoded `FilletOptions::default()`). Now the recipe drives the
/// curvature-continuous construction and the measured evidence says
/// so on both sides.
#[test]
fn g2_recipe_replays_to_part_parity_with_measured_g2() {
    let recipe = FilletRecipeDesc {
        continuity: ContinuityKind::G2,
        ..Default::default()
    };
    let (cube, replayed) = replay_fillet(recipe, 0.1);
    let part = part_fillet(&cube, |f| f.radius(0.1).continuity(Continuity::G2));
    assert_same_body("G2 fillet", &replayed, part.brep(), 1.0e-9);

    let replay_evidence = evidence(&cube, &replayed, Continuity::G2);
    let part_evidence = part
        .last_blend_report()
        .and_then(|r| r.continuity.as_ref())
        .expect("part records evidence");
    assert_eq!(
        replay_evidence.achieved,
        popkern_blend::ContinuityLevel::G2,
        "replayed G2 recipe must measure G2: {replay_evidence:?}"
    );
    assert_eq!(replay_evidence.achieved, part_evidence.achieved);
    assert_eq!(replay_evidence.blend_faces, part_evidence.blend_faces);

    // The default recipe on the same fixture is a different (G1) body.
    let (_, g1) = replay_fillet(FilletRecipeDesc::default(), 0.1);
    let g1_evidence = evidence(&cube, &g1, Continuity::G2);
    assert_ne!(
        g1_evidence.achieved,
        popkern_blend::ContinuityLevel::G2,
        "the default recipe must not silently deliver G2"
    );
}

#[test]
fn variable_radius_recipe_replays_to_part_parity() {
    let (mut tree, cube) = cube_tree();
    let edges = vertical_edge_names(&cube);
    let r0 = add_param(&mut tree, "_r_start", 0.08);
    let r1 = add_param(&mut tree, "_r_end", 0.16);
    tree.add_feature(
        Feature::Fillet {
            edges,
            radius: r0,
            recipe: FilletRecipeDesc {
                radius_law: FilletRadiusLaw::Linear { end: r1 },
                ..Default::default()
            },
        },
        "VariableFillet",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let replayed = loaded.current_brep().expect("body").clone();
    replayed.validate().expect("validates");

    let part = part_fillet(&cube, |f| f.variable_radius(0.08, 0.16));
    assert_same_body("variable fillet", &replayed, part.brep(), 1.0e-9);

    // The law is honoured, not collapsed to constant: a constant
    // 0.08 fillet is a different body.
    let (_, constant) = replay_fillet(FilletRecipeDesc::default(), 0.08);
    let t = Tolerance::default();
    let v_var = popkern_measure::body_volume(&replayed, &t).expect("v");
    let v_const = popkern_measure::body_volume(&constant, &t).expect("v");
    assert!(
        (v_var - v_const).abs() > 1.0e-6,
        "variable law must change the body: {v_var} vs {v_const}"
    );
}

#[test]
fn conic_profile_recipe_replays_to_part_parity() {
    let recipe = FilletRecipeDesc {
        profile: BlendProfileKind::Conic { rho: 0.4 },
        ..Default::default()
    };
    let (cube, replayed) = replay_fillet(recipe, 0.1);
    let part = part_fillet(&cube, |f| f.radius(0.1).conic(0.4));
    assert_same_body("conic fillet", &replayed, part.brep(), 1.0e-9);

    let (_, circular) = replay_fillet(FilletRecipeDesc::default(), 0.1);
    let t = Tolerance::default();
    let v_conic = popkern_measure::body_volume(&replayed, &t).expect("v");
    let v_circ = popkern_measure::body_volume(&circular, &t).expect("v");
    assert!(
        (v_conic - v_circ).abs() > 1.0e-6,
        "rho=0.4 conic must differ from the circular section: {v_conic} vs {v_circ}"
    );
}

/// Chamfer kinds: the equal bevel matches `Part::chamfer`; the
/// two-distance and distance-angle kinds resolve their reference face
/// by stable name and match the direct kernel call with the same
/// specs (popkern-part exposes only the equal bevel today).
#[test]
fn chamfer_kinds_replay_to_kernel_parity() {
    let (mut tree, cube) = cube_tree();
    let edges = vertical_edge_names(&cube);
    let d = add_param(&mut tree, "_d", 0.1);
    tree.add_feature(
        Feature::Chamfer {
            edges: edges.clone(),
            distance: d,
            kind: ChamferKindDesc::Equal,
        },
        "Chamfer1",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let replayed = loaded.current_brep().expect("body").clone();
    let part = Part::from_brep(cube.clone(), Tolerance::default(), "cube")
        .expect("part")
        .chamfer(0.1, |e| e.straight().direction(Dir::Z))
        .expect("part chamfer");
    assert_same_body("equal chamfer", &replayed, part.brep(), 1.0e-9);

    // Two-distance on one vertical edge, referenced to the face whose
    // normal is +X (the x=1 side).
    let edge_name = edges
        .iter()
        .find(|n| {
            let h = cube
                .edges
                .iter()
                .find(|(_, e)| e.name.map(|id| cube.name_pool.get(id) == Some(*n)) == Some(true))
                .map(|(h, _)| h)
                .expect("edge");
            edge_endpoint_max_x(&cube, h) > 0.5 && edge_endpoint_max_y(&cube, h) > 0.5
        })
        .expect("x=1,y=1 vertical edge")
        .clone();
    let (face_name, face_handle) = {
        let (h, f) = cube
            .faces
            .iter()
            .find(|(h, _)| {
                let n = popkern_measure::face_normal(&cube, *h).expect("normal");
                n.x > 0.99
            })
            .expect("+X face");
        (
            cube.name_pool
                .get(f.name.expect("named face"))
                .expect("name")
                .clone(),
            h,
        )
    };
    let (mut tree, _) = cube_tree();
    let d1 = add_param(&mut tree, "_d1", 0.1);
    let d2 = add_param(&mut tree, "_d2", 0.2);
    tree.add_feature(
        Feature::Chamfer {
            edges: vec![edge_name.clone()],
            distance: d1,
            kind: ChamferKindDesc::TwoDistance {
                face: face_name.clone(),
                distance_2: d2,
            },
        },
        "TwoDistance",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let replayed = loaded.current_brep().expect("body").clone();
    let edge_handle = cube
        .edges
        .iter()
        .find(|(_, e)| e.name.map(|id| cube.name_pool.get(id) == Some(&edge_name)) == Some(true))
        .map(|(h, _)| h)
        .expect("edge handle");
    let direct = popkern_blend::chamfer_edges(
        &cube,
        &[popkern_blend::ChamferSpec {
            edge: edge_handle,
            kind: popkern_blend::ChamferKind::TwoDistance {
                face_1: face_handle,
                dist_1: 0.1,
                dist_2: 0.2,
            },
        }],
        &popkern_blend::ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("direct chamfer")
    .brep;
    assert_same_body("two-distance chamfer", &replayed, &direct, 1.0e-9);
    // Asymmetric: swapping which face carries 0.2 is a different body.
    let equal = popkern_blend::chamfer_edges(
        &cube,
        &[popkern_blend::ChamferSpec {
            edge: edge_handle,
            kind: popkern_blend::ChamferKind::Equal(0.1),
        }],
        &popkern_blend::ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("equal chamfer")
    .brep;
    let t = Tolerance::default();
    let v_two = popkern_measure::body_volume(&replayed, &t).expect("v");
    let v_eq = popkern_measure::body_volume(&equal, &t).expect("v");
    assert!(
        (v_two - v_eq).abs() > 1.0e-6,
        "two-distance must differ from equal"
    );

    // Distance-angle on the same edge and face.
    let (mut tree, _) = cube_tree();
    let d1 = add_param(&mut tree, "_d1", 0.1);
    let ang = add_param(&mut tree, "_ang", std::f64::consts::FRAC_PI_6);
    tree.add_feature(
        Feature::Chamfer {
            edges: vec![edge_name.clone()],
            distance: d1,
            kind: ChamferKindDesc::DistanceAngle {
                face: face_name,
                angle: ang,
            },
        },
        "DistanceAngle",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let replayed = loaded.current_brep().expect("body").clone();
    let direct = popkern_blend::chamfer_edges(
        &cube,
        &[popkern_blend::ChamferSpec {
            edge: edge_handle,
            kind: popkern_blend::ChamferKind::DistanceAngle {
                face: face_handle,
                distance: 0.1,
                angle: std::f64::consts::FRAC_PI_6,
            },
        }],
        &popkern_blend::ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("direct chamfer")
    .brep;
    assert_same_body("distance-angle chamfer", &replayed, &direct, 1.0e-9);
}

fn edge_endpoint_max_x(brep: &BRep, edge: popkern_core::topology::EdgeHandle) -> f64 {
    edge_endpoints(brep, edge)
        .iter()
        .map(|p| p[0])
        .fold(f64::NEG_INFINITY, f64::max)
}

fn edge_endpoint_max_y(brep: &BRep, edge: popkern_core::topology::EdgeHandle) -> f64 {
    edge_endpoints(brep, edge)
        .iter()
        .map(|p| p[1])
        .fold(f64::NEG_INFINITY, f64::max)
}

fn edge_endpoints(brep: &BRep, edge: popkern_core::topology::EdgeHandle) -> Vec<[f64; 3]> {
    let he = brep.edges[edge].half_edge;
    let twin = brep.half_edges[he].twin;
    [brep.half_edges[he].origin, brep.half_edges[twin].origin]
        .iter()
        .map(|&v| {
            let p = brep.vertices[v].position;
            [p.x, p.y, p.z]
        })
        .collect()
}
