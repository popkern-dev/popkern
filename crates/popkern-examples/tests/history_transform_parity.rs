//! PK-362 acceptance: a `Feature::Transform` stored in a feature tree
//! replays to the same body as the equivalent `popkern-part` call
//! (`translate`, `rotate_z`, `transformed`) on a shared fixture, and
//! survives a `.pkn` save/load with the motion honoured.
//!
//! Parity is asserted through numeric gates — face/edge/vertex census,
//! sorted vertex positions, body volume — never visual inspection.
//! Stable-name survival is asserted by filleting edges named BEFORE
//! the transform, on a body carrying rotational charts (cylinder
//! fillets), which exercises `BRep::transform`'s chart re-anchoring
//! under a strict `validate_all`.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::*;
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

fn add_param(tree: &mut FeatureTree, name: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, unit)
        .expect("add param")
}

/// Unit cube via sketch + extrude, replayed. Returns the tree and the
/// pre-transform body.
fn cube_tree() -> (FeatureTree, BRep) {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let dist = add_param(&mut tree, "_extrude_distance", 1.0, Unit::Meters);
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
    // PK-376: the sketch→extrude body is flag-honest under the strict
    // validator (its base cap used to carry a dishonest flag).
    let errs = brep.validate_all();
    assert!(
        errs.is_empty(),
        "sketch→extrude cube validates strictly: {errs:?}"
    );
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

/// Add a transform feature, save/load, replay; returns the replayed
/// body (validated with the strict path).
fn replay_transform(mut tree: FeatureTree, motion: RigidMotionDesc) -> BRep {
    tree.add_feature(Feature::Transform { motion }, "Transform1");
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let body = loaded.current_brep().expect("body").clone();
    let errs = body.validate_all();
    assert!(errs.is_empty(), "replayed body validates: {errs:?}");
    body
}

fn cube_part(cube: &BRep) -> Part {
    Part::from_brep(cube.clone(), Tolerance::default(), "cube").expect("part")
}

#[test]
fn translate_replays_to_part_parity() {
    let (mut tree, cube) = cube_tree();
    let dx = add_param(&mut tree, "dx", 2.0, Unit::Meters);
    let dy = add_param(&mut tree, "dy", -1.5, Unit::Meters);
    let dz = add_param(&mut tree, "dz", 0.25, Unit::Meters);
    let replayed = replay_transform(tree, RigidMotionDesc::Translate { dx, dy, dz });
    let part = cube_part(&cube)
        .translate(2.0, -1.5, 0.25)
        .expect("part translate");
    assert_same_body("translate", &replayed, part.brep(), 1.0e-12);
    // The motion actually moved the body.
    assert!((sorted_vertices(&replayed)[0][0] - 2.0).abs() < 1e-12);
}

#[test]
fn rotate_about_literal_axis_replays_to_part_parity() {
    let (mut tree, cube) = cube_tree();
    let ang = add_param(&mut tree, "ang", 0.6, Unit::Radians);
    let replayed = replay_transform(
        tree,
        RigidMotionDesc::Rotate {
            axis: AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            angle: ang,
        },
    );
    let part = cube_part(&cube).rotate_z(0.6).expect("part rotate_z");
    assert_same_body("rotate_z", &replayed, part.brep(), 1.0e-12);
}

/// The rotation axis may follow a model edge: rotating the cube about
/// one of its own vertical edges is the same body `Part::transformed`
/// produces with the explicit `rotation_around_axis`.
#[test]
fn rotate_about_edge_axis_replays_to_part_parity() {
    let (mut tree, cube) = cube_tree();
    let edge = vertical_edge_names(&cube)
        .into_iter()
        .next()
        .expect("vertical edge");
    // Locate the edge's line for the Part-side transform.
    let edge_handle = cube
        .edges
        .iter()
        .find(|(_, e)| {
            e.name
                .map(|nid| cube.name_pool.get(nid) == Some(&edge))
                .unwrap_or(false)
        })
        .map(|(h, _)| h)
        .expect("edge handle");
    let he = &cube.half_edges[cube.edges[edge_handle].half_edge];
    let start = cube.vertices[he.origin].position;
    let end = cube.vertices[cube.half_edges[he.twin].origin].position;
    let ang = add_param(&mut tree, "ang", 1.1, Unit::Radians);
    let replayed = replay_transform(
        tree,
        RigidMotionDesc::Rotate {
            axis: AxisDesc::Edge(edge),
            angle: ang,
        },
    );
    let t = Transform3::rotation_around_axis(start, end - start, 1.1).expect("axis transform");
    let part = cube_part(&cube)
        .transformed(&t, "rotate about edge".into())
        .expect("part transformed");
    assert_same_body("rotate about edge", &replayed, part.brep(), 1.0e-12);
    // The axis line is fixed under the motion.
    let fixed = t.apply_point(&start);
    assert!((fixed - start).length() < 1e-12);
}

#[test]
fn literal_matrix_replays_to_part_parity() {
    let (tree, cube) = cube_tree();
    let t = Transform3::rotation_around_axis(
        Point3::new(0.3, -0.2, 0.7),
        Vec3::new(1.0, 1.0, 0.5),
        0.9,
    )
    .expect("axis transform")
    .compose(&Transform3::translation(1.0, 2.0, 3.0));
    let replayed = replay_transform(tree, RigidMotionDesc::Matrix(t.matrix));
    let part = cube_part(&cube)
        .transformed(&t, "general".into())
        .expect("part transformed");
    assert_same_body("matrix", &replayed, part.brep(), 1.0e-12);
}

/// A non-similarity matrix declines transactionally at replay: the
/// feature reports `Failed`, the tree does not panic, and the body
/// left in the tree is untouched (the pre-transform cube).
#[test]
fn shear_matrix_declines_transactionally() {
    let (mut tree, cube) = cube_tree();
    tree.add_feature(
        Feature::Transform {
            motion: RigidMotionDesc::Matrix([
                [1.0, 0.4, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
        },
        "Shear",
    );
    let result = tree.replay(Tolerance::default()).expect("replay runs");
    let failed = result
        .feature_results
        .values()
        .filter_map(|r| match r {
            FeatureResult::Failed { error } => Some(error.clone()),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(failed.len(), 1, "exactly the shear fails: {failed:?}");
    assert!(failed[0].contains("similarity"), "{}", failed[0]);
    let body = tree.current_brep().expect("body");
    assert_same_body("untouched cube", body, &cube, 1.0e-12);
}

/// Stable names survive the transform: edges named on the cube are
/// filleted AFTER a rotation, on both sides, and the rotated cylinder
/// charts pass the strict validator. Parity vs `Part::rotate_z` then
/// `Part::fillet`.
#[test]
fn names_survive_transform_and_rotational_charts_validate() {
    let (mut tree, cube) = cube_tree();
    let edges = vertical_edge_names(&cube);
    let ang = add_param(&mut tree, "ang", 0.8, Unit::Radians);
    let r = add_param(&mut tree, "r", 0.15, Unit::Meters);
    tree.add_feature(
        Feature::Transform {
            motion: RigidMotionDesc::Rotate {
                axis: AxisDesc::Literal {
                    origin: [0.0, 0.0, 0.0],
                    direction: [0.0, 0.0, 1.0],
                },
                angle: ang,
            },
        },
        "Rotate",
    );
    tree.add_feature(
        Feature::Fillet {
            edges,
            radius: r,
            recipe: Default::default(),
        },
        "Fillet1",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let replayed = loaded.current_brep().expect("body").clone();
    let errs = replayed.validate_all();
    assert!(
        errs.is_empty(),
        "rotated + filleted body validates strictly: {errs:?}"
    );
    assert_eq!(replayed.faces.len(), 10, "four cylinder fillet faces");

    let part = cube_part(&cube)
        .rotate_z(0.8)
        .expect("rotate")
        .fillet(0.15, |e| e.straight().direction(Dir::Z))
        .expect("fillet after rotate");
    assert_same_body("rotate then fillet", &replayed, part.brep(), 1.0e-9);
}
