//! Tests for the label table, feature-label uniqueness, the build
//! tolerance attribute, and the naming-scheme envelope field.

use super::*;
use crate::feature::{ExtrudeEnd, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::feature_sketch::SketchPayload;
use crate::unit::Unit;
use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

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
    sk.solve().expect("solve");
    serde_json::to_string(&sk).expect("serialize")
}

fn sketch() -> Feature {
    Feature::Sketch {
        plane: SketchPlaneDesc::XY,
        sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
    }
}

fn name(op: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

/// Sketch (op 1) + extrude (op 2) block. Returns the tree and extrude id.
fn block(sketch_label: &str, extrude_label: &str) -> (FeatureTree, FeatureId) {
    let mut tree = FeatureTree::new();
    tree.add_feature(sketch(), sketch_label);
    let dist = tree
        .parameters_mut()
        .add_global("dist", 1.0, Unit::Meters)
        .expect("dist");
    let extrude = tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        extrude_label,
    );
    (tree, extrude)
}

fn sorted_vertices(tree: &FeatureTree) -> Vec<[f64; 3]> {
    let brep = tree.current_brep().expect("body");
    let mut pts: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    pts
}

#[test]
fn duplicate_or_empty_labels_receive_unique_auto_labels() {
    let mut tree = FeatureTree::new();
    let a = tree.add_feature(sketch(), "Sketch");
    let b = tree.add_feature(sketch(), "Sketch");
    let c = tree.add_feature(sketch(), "");
    let labels: Vec<String> = tree.feature_labels().into_iter().map(|(l, _)| l).collect();
    assert_eq!(labels, ["Sketch", "sketch-1", "sketch-2"]);
    assert_eq!(tree.feature_by_label("Sketch"), Some(a));
    assert_eq!(tree.feature_by_label("sketch-1"), Some(b));
    assert_eq!(tree.feature_by_label("sketch-2"), Some(c));
    assert_eq!(tree.labels().auto_counters.get("sketch"), Some(&3));
}

#[test]
fn auto_label_counters_persist_and_never_reuse_numbers() {
    let mut tree = FeatureTree::new();
    tree.add_feature(sketch(), "");
    let second = tree.add_feature(sketch(), "");
    tree.delete_feature(second, crate::DeleteBehavior::Cascade)
        .expect("delete");
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let third = loaded.add_feature(sketch(), "");
    assert_eq!(
        loaded.feature(third).map(|n| n.label.as_str()),
        Some("sketch-3"),
        "a deleted feature's number is never minted again"
    );
}

#[test]
fn rename_edits_the_label_only_and_rekeys_scoped_parameters() {
    let (mut tree, extrude) = block("Base Sketch", "Main Block");
    let before = tree.feature(extrude).expect("node").clone();
    assert!(
        tree.parameters()
            .resolve_name("main_block.distance")
            .is_some()
    );
    tree.replay(Tolerance::default()).expect("replay");
    let body_before = sorted_vertices(&tree);

    tree.rename_feature(extrude, "Tower").expect("rename");

    let after = tree.feature(extrude).expect("node");
    assert_eq!(after.label, "Tower");
    assert_eq!(after.operation_id, before.operation_id);
    assert_eq!(after.depends_on, before.depends_on);
    assert_eq!(after.body, before.body);
    assert_eq!(tree.feature_by_label("Tower"), Some(extrude));
    assert_eq!(tree.feature_by_label("Main Block"), None);
    assert!(
        tree.parameters().resolve_name("tower.distance").is_some(),
        "scoped parameters re-key under the new label"
    );
    assert!(
        tree.parameters()
            .resolve_name("main_block.distance")
            .is_none()
    );

    tree.replay(Tolerance::default())
        .expect("replay after rename");
    assert_eq!(sorted_vertices(&tree), body_before, "geometry untouched");

    assert!(matches!(
        tree.rename_feature(extrude, "Base Sketch"),
        Err(HistoryError::DuplicateLabel(_))
    ));
    assert!(matches!(
        tree.rename_feature(extrude, ""),
        Err(HistoryError::EmptyLabel)
    ));

    assert!(tree.undo());
    assert_eq!(tree.feature(extrude).expect("node").label, "Main Block");
    assert!(
        tree.parameters()
            .resolve_name("main_block.distance")
            .is_some()
    );
    assert!(tree.redo());
    assert_eq!(tree.feature(extrude).expect("node").label, "Tower");
    assert!(tree.parameters().resolve_name("tower.distance").is_some());
}

#[test]
fn entity_labels_are_an_indirection_over_canonical_names() {
    let mut tree = FeatureTree::new();
    let top = name(2, RoleCategory::Generated, 0);
    let side = name(2, RoleCategory::Lateral, 3);
    tree.label_entity("top", top.clone()).expect("label");
    tree.label_entity("top", top.clone())
        .expect("rebinding the same entity is a no-op");
    assert!(matches!(
        tree.label_entity("top", side.clone()),
        Err(HistoryError::DuplicateLabel(_))
    ));
    assert!(matches!(
        tree.label_entity("", side.clone()),
        Err(HistoryError::EmptyLabel)
    ));
    assert_eq!(tree.entity_by_label("top"), Some(&top));
    assert_eq!(tree.entity_label(&top), Some("top"));

    tree.rename_entity_label("top", "lid").expect("rename");
    assert_eq!(tree.entity_by_label("lid"), Some(&top), "name unchanged");
    assert_eq!(tree.entity_by_label("top"), None);
    assert!(matches!(
        tree.rename_entity_label("missing", "x"),
        Err(HistoryError::LabelNotFound(_))
    ));

    assert!(tree.undo(), "rename is one undo step");
    assert_eq!(tree.entity_by_label("top"), Some(&top));
    assert_eq!(tree.entity_by_label("lid"), None);
    assert!(tree.redo());
    assert_eq!(tree.entity_by_label("lid"), Some(&top));

    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    assert_eq!(loaded.labels(), tree.labels());
    assert_eq!(loaded.unlabel_entity("lid").expect("unlabel"), top);
    assert!(matches!(
        loaded.unlabel_entity("lid"),
        Err(HistoryError::LabelNotFound(_))
    ));
    assert!(loaded.undo());
    assert_eq!(loaded.entity_by_label("lid"), Some(&top));
}

#[test]
fn build_tolerance_round_trips_and_reaches_replay() {
    let (mut tree, _) = block("S", "E");
    assert_eq!(tree.build_tolerance(), None);
    let result = tree.replay_document().expect("replay");
    assert_eq!(result.tolerance, Tolerance::default());

    let loose = Tolerance::default()
        .with_linear(3.5e-7)
        .with_angular(2.0e-4);
    tree.set_build_tolerance(Some(loose));
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    assert_eq!(loaded.build_tolerance(), Some(loose), "bit-exact");
    let result = loaded.replay_document().expect("replay");
    assert_eq!(result.tolerance, loose, "document tolerance reached replay");
    let explicit = loaded.replay(Tolerance::default()).expect("replay");
    assert_eq!(
        explicit.tolerance,
        Tolerance::default(),
        "explicit override wins"
    );

    assert!(tree.undo());
    assert_eq!(tree.build_tolerance(), None);
    assert!(tree.redo());
    assert_eq!(tree.build_tolerance(), Some(loose));
}

#[test]
fn naming_scheme_rides_the_envelope() {
    // A new tree mints under the current scheme (3, PK-731) and says so
    // in the envelope and in the tree.
    let tree = FeatureTree::new();
    assert_eq!(tree.naming_scheme(), 3);
    let json = String::from_utf8(tree.save().expect("save")).expect("utf8");
    assert!(json.contains("\"naming_scheme\":3"), "{json}");

    let newer = json.replace("\"naming_scheme\":3", "\"naming_scheme\":99");
    assert!(matches!(
        FeatureTree::load(newer.as_bytes()),
        Err(HistoryError::UnsupportedNamingScheme {
            file_scheme: 99,
            supported: 3
        })
    ));

    // A document written before either field existed is scheme 1, keeps
    // scheme 1's rules, and saves back as scheme 1.
    let unrecorded = json
        .replace("\"naming_scheme\":3,", "")
        .replace(",\"naming_scheme\":3", "");
    assert!(!unrecorded.contains("naming_scheme"), "{unrecorded}");
    let old =
        FeatureTree::load(unrecorded.as_bytes()).expect("documents without the field are scheme 1");
    assert_eq!(old.naming_scheme(), 1);
    let resaved = String::from_utf8(old.save().expect("save")).expect("utf8");
    assert!(resaved.contains("\"naming_scheme\":1"), "{resaved}");
}

#[test]
fn load_makes_version_3_duplicate_labels_unique_deterministically() {
    let mut tree = FeatureTree::new();
    tree.add_feature(sketch(), "Sketch 1");
    tree.add_feature(sketch(), "Sketch 1");
    let json = String::from_utf8(tree.save().expect("save")).expect("utf8");
    assert!(json.contains("\"label\":\"sketch-1\""));
    // A version-3 writer stored the duplicate verbatim.
    let unrecorded = json.replace("\"label\":\"sketch-1\"", "\"label\":\"Sketch 1\"");
    let loaded = FeatureTree::load(unrecorded.as_bytes()).expect("load");
    let labels: Vec<String> = loaded
        .feature_labels()
        .into_iter()
        .map(|(l, _)| l)
        .collect();
    assert_eq!(labels[0], "Sketch 1", "the first keeps its label");
    assert_ne!(labels[1], "Sketch 1");
    assert!(labels[1].starts_with("sketch-"), "{labels:?}");
    let again = FeatureTree::load(unrecorded.as_bytes()).expect("load");
    assert_eq!(again.feature_labels(), loaded.feature_labels());
}
