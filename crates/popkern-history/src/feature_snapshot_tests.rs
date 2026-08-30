//! Tests for the `RawBody` snapshot feature: materialisation, name
//! ownership, determinism across save/load, downstream addressing, and
//! the decline path.

use super::*;
use crate::feature::Feature;
use crate::feature_transform::RigidMotionDesc;
use crate::unit::Unit;
use crate::{FeatureResult, FeatureTree, ReplayOutcome};
use popkern_core::naming::OperationId;

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("unit cube")
}

fn raw_body(snapshot: BRep) -> Feature {
    Feature::RawBody {
        snapshot: Box::new(snapshot),
        origin: SnapshotOrigin::Import {
            source: "unit cube".into(),
        },
    }
}

fn names_of(brep: &BRep) -> Vec<EntityName> {
    let mut names: Vec<EntityName> = brep
        .faces
        .values()
        .map(|f| f.name)
        .chain(brep.edges.values().map(|e| e.name))
        .chain(brep.vertices.values().map(|v| v.name))
        .map(|id| {
            brep.name_pool
                .get(id.expect("every entity is named"))
                .expect("pool")
                .clone()
        })
        .collect();
    names.sort_by_key(|n| (n.role.category as u8, n.role.index));
    names
}

#[test]
fn snapshot_materialises_as_a_new_body_owning_every_name() {
    let mut tree = FeatureTree::new();
    let mut snapshot = cube();
    // Names the snapshot arrives with belong to nobody in this document.
    snapshot.name_pool.intern(
        EntityName {
            operation: OperationId(77),
            role: EntityRole {
                category: RoleCategory::Generated,
                index: 0,
            },
            ancestors: Vec::new(),
        },
        None,
    );
    let id = tree.add_feature(raw_body(snapshot), "");
    assert_eq!(tree.feature(id).expect("node").label, "raw_body-1");
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(result.report.is_clean(), "{:?}", result.report);
    let body = tree.current_brep().expect("body");
    assert_eq!(
        (body.faces.len(), body.edges.len(), body.vertices.len()),
        (6, 12, 8)
    );
    body.validate().expect("served body validates");
    let op = tree.feature(id).expect("node").operation_id;
    let names = names_of(body);
    assert_eq!(names.len(), 26);
    assert!(
        names
            .iter()
            .all(|n| n.operation == op && n.ancestors.is_empty())
    );
    assert!(
        body.name_pool.names.iter().all(|n| n.operation == op),
        "inherited names are discarded"
    );
    let faces: Vec<u32> = names
        .iter()
        .filter(|n| n.role.category == RoleCategory::Generated)
        .map(|n| n.role.index)
        .collect();
    assert_eq!(faces, (0..6).collect::<Vec<_>>());
    match tree.feature(id).and_then(|n| n.last_result.as_ref()) {
        Some(FeatureResult::Success { history }) => assert_eq!(history.generated.len(), 26),
        other => panic!("expected success, got {other:?}"),
    }
}

#[test]
fn snapshot_names_are_identical_across_save_load_and_replay() {
    let mut tree = FeatureTree::new();
    tree.add_feature(raw_body(cube()), "Imported");
    tree.replay(Tolerance::default()).expect("replay");
    let direct = tree.current_brep().expect("body").clone();

    let bytes = tree.save().expect("save");
    assert!(String::from_utf8_lossy(&bytes).contains("\"RawBody\""));
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    loaded.replay(Tolerance::default()).expect("replay");
    let replayed = loaded.current_brep().expect("body");

    assert_eq!(names_of(replayed), names_of(&direct));
    let pos = |b: &BRep| -> Vec<[f64; 3]> {
        let mut p: Vec<[f64; 3]> = b
            .vertices
            .values()
            .map(|v| [v.position.x, v.position.y, v.position.z])
            .collect();
        p.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
        p
    };
    assert_eq!(pos(replayed), pos(&direct), "bit-identical geometry");
}

#[test]
fn downstream_features_address_snapshot_entities_by_name() {
    let mut tree = FeatureTree::new();
    let raw = tree.add_feature(raw_body(cube()), "Imported");
    let op = tree.feature(raw).expect("node").operation_id;
    let dx = tree
        .add_global_parameter("dx", 2.0, Unit::Meters)
        .expect("dx");
    let zero = tree
        .add_global_parameter("zero", 0.0, Unit::Meters)
        .expect("zero");
    let moved = tree.add_feature(
        Feature::Transform {
            motion: RigidMotionDesc::Translate {
                dx,
                dy: zero,
                dz: zero,
            },
        },
        "Shift",
    );
    // A reference measurement on a snapshot face resolves through its name.
    let top = EntityName {
        operation: op,
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    };
    let reference = tree.add_feature(
        Feature::Reference {
            name: "area".into(),
            query: crate::MeasureQuery::FaceArea(top),
        },
        "Area",
    );
    assert_eq!(
        tree.feature(moved).expect("node").body,
        tree.feature(raw).expect("node").body
    );
    assert_eq!(tree.feature(reference).expect("node").depends_on, vec![raw]);
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(result.report.is_clean(), "{:?}", result.report);
    let body = tree.current_brep().expect("body");
    let min_x = body
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::INFINITY, f64::min);
    assert_eq!(min_x, 2.0, "translated snapshot");
    let area = tree
        .parameters()
        .resolve_name("area.value")
        .and_then(|p| tree.parameters().get(p))
        .map(|p| p.value)
        .expect("reference value");
    assert!(
        (area - 1.0).abs() < 1e-12,
        "unit cube face area, got {area}"
    );
}

#[test]
fn empty_snapshot_declines_and_reports_failed() {
    let mut tree = FeatureTree::new();
    let id = tree.add_feature(raw_body(BRep::new(Tolerance::default())), "Empty");
    let result = tree
        .replay(Tolerance::default())
        .expect("replay is not fatal");
    match &result.report.feature(id).expect("report").outcome {
        ReplayOutcome::Failed { error } => assert!(error.contains("no faces"), "{error}"),
        other => panic!("expected Failed, got {other:?}"),
    }
}
