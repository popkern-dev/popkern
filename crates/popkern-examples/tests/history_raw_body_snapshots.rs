//! `Feature::RawBody` snapshot gate (design row 1g, PK-367): bodies the
//! part API builds without a parametric recipe — imports, polyhedra,
//! smooth-mesh solids, engraved results — record as snapshots, survive
//! save/load, replay to the identical body, and expose their entities
//! to downstream features by stable name.

use popkern_core::naming::{EntityName, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::{Feature, FeatureTree, FilletRecipeDesc, ReplayResult, SnapshotOrigin, Unit};
use popkern_part::{
    EngravePlacement, HoledPolygon, Part, SmoothMeshOptions, engrave_polygons, face_frame,
    smooth_solid_from_obj,
};

fn sorted_vertices(brep: &BRep) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    pts
}

fn assert_replay_clean(result: &ReplayResult) {
    assert!(
        result.report.is_clean(),
        "replay report: {:?}",
        result.report
    );
}

/// Record `part` as a snapshot, save, load, replay; return the served
/// body and the tree (for downstream tests).
fn record_snapshot(part: &Part, origin: SnapshotOrigin) -> (FeatureTree, BRep) {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::RawBody {
            snapshot: Box::new(part.brep().clone()),
            origin,
        },
        "Snapshot",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let body = loaded.current_brep().expect("body").clone();
    body.validate().expect("replayed snapshot validates");
    (loaded, body)
}

fn assert_same_body(label: &str, a: &BRep, b: &BRep) {
    assert_eq!(a.faces.len(), b.faces.len(), "{label}: faces");
    assert_eq!(a.edges.len(), b.edges.len(), "{label}: edges");
    assert_eq!(a.vertices.len(), b.vertices.len(), "{label}: vertices");
    assert_eq!(
        sorted_vertices(a),
        sorted_vertices(b),
        "{label}: vertex positions"
    );
}

fn cube_points(size: f64) -> Vec<popkern_core::geometry::Point3> {
    let mut pts = Vec::new();
    for z in [0.0, size] {
        for y in [0.0, size] {
            for x in [0.0, size] {
                pts.push(popkern_core::geometry::Point3::new(x, y, z));
            }
        }
    }
    pts
}

#[test]
fn polyhedron_records_as_snapshot_and_serves_a_downstream_fillet() {
    let part = Part::polyhedron(&cube_points(1.0), Tolerance::default()).expect("polyhedron");
    let (mut tree, replayed) = record_snapshot(&part, SnapshotOrigin::Polyhedron);
    assert_same_body("polyhedron", part.brep(), &replayed);

    // Fillet the four vertical edges by their snapshot names.
    let mut edges: Vec<EntityName> = Vec::new();
    for (handle, edge) in &replayed.edges {
        let tangent = popkern_measure::edge_tangent(&replayed, handle).expect("tangent");
        if tangent.z.abs() > 0.99 {
            let name = replayed
                .name_pool
                .get(edge.name.expect("snapshot edges are named"))
                .expect("pool")
                .clone();
            assert_eq!(name.role.category, RoleCategory::Lateral);
            edges.push(name);
        }
    }
    assert_eq!(edges.len(), 4);
    let radius = tree
        .add_global_parameter("r", 0.1, Unit::Meters)
        .expect("radius");
    tree.add_feature(
        Feature::Fillet {
            edges,
            radius,
            recipe: FilletRecipeDesc::default(),
        },
        "Rounds",
    );
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let filleted = tree.current_brep().expect("body");
    filleted.validate().expect("filleted snapshot validates");
    assert_eq!(filleted.faces.len(), 10, "6 faces + 4 fillet strips");
    let evidence = result
        .report
        .feature(tree.feature_by_label("Rounds").expect("fillet"))
        .and_then(|r| r.blend.as_ref())
        .expect("fillet evidence");
    assert!(!evidence.downgraded(), "{evidence:?}");
}

#[test]
fn imported_body_records_as_snapshot() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let part = Part::from_brep(cube, Tolerance::default(), "unit cube").expect("import");
    let (_, replayed) = record_snapshot(
        &part,
        SnapshotOrigin::Import {
            source: "unit cube".into(),
        },
    );
    assert_same_body("import", part.brep(), &replayed);
}

#[test]
fn smooth_mesh_solid_records_as_snapshot() {
    let cage = "\
v -1 -1 -1
v 1 -1 -1
v 1 1 -1
v -1 1 -1
v -1 -1 1
v 1 -1 1
v 1 1 1
v -1 1 1
f 1 4 3 2
f 5 6 7 8
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";
    let (part, _report) =
        smooth_solid_from_obj(cage, Tolerance::default(), SmoothMeshOptions::default())
            .expect("smooth cube");
    let (_, replayed) = record_snapshot(&part, SnapshotOrigin::SmoothMesh);
    assert_same_body("smooth mesh", part.brep(), &replayed);
}

#[test]
fn engraved_result_records_as_snapshot() {
    let block = Part::polyhedron(&cube_points(10.0), Tolerance::default()).expect("block");
    let square = HoledPolygon {
        outer: vec![(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)],
        holes: Vec::new(),
    };
    let engraved = block
        .edit_polyhedron("engrave a square pocket", |points, faces| {
            let top = (0..faces.len())
                .find(|&i| face_frame(points, &faces[i]).normal.z > 0.9)
                .expect("top face");
            let placement = EngravePlacement::centered(
                &face_frame(points, &faces[top]),
                std::slice::from_ref(&square),
                0.4,
            );
            engrave_polygons(
                points,
                faces,
                top,
                std::slice::from_ref(&square),
                &placement,
                0.5,
            )
        })
        .expect("engrave");
    assert!(engraved.brep().faces.len() > 6, "pocket adds faces");
    let (_, replayed) = record_snapshot(&engraved, SnapshotOrigin::Engrave);
    assert_same_body("engrave", engraved.brep(), &replayed);
}
