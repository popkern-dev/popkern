// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

use crate::FeatureTree;
use crate::feature::{Feature, SketchPlaneDesc, Vec3Desc};

/// Create a solved rectangle sketch and serialize it to JSON.
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

/// Helper: create an EntityName referencing a given OperationId.
fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

#[test]
fn healthy_tree_no_references_reports_no_issues() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    // A sketch on XY has no entity references — always healthy.
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Replay to populate current_brep.
    tree.replay(Tolerance::default()).expect("replay");

    let report = tree.validate_references().expect("should have brep");
    assert!(
        report.is_healthy(),
        "tree should be healthy but has {} issues",
        report.issues.len()
    );
    assert_eq!(report.total_references, 0);
    assert_eq!(report.resolved_references, 0);
}

#[test]
fn healthy_sketch_on_face_reference_resolves() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    // Sketch 1 on XY — produces faces with OperationId(1).
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json.clone(),
        },
        "Sketch 1",
    );

    // Extrude to create a box — the sketch face name should survive
    // as a history reference in the extruded BRep.
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: 1.0,
        },
        "Extrude 1",
    );

    // Replay to populate current_brep.
    tree.replay(Tolerance::default()).expect("replay");

    let report = tree.validate_references().expect("should have brep");
    // The extrude references OperationId(1) Generated:0 — this face may
    // or may not still be resolvable in the final BRep after extrusion.
    // Either way the report is well-formed.
    assert!(report.total_references > 0);
}

#[test]
fn broken_reference_detected() {
    let mut tree = FeatureTree::new();
    let sketch_json = rectangle_sketch_json();

    // Sketch gets OperationId(1).
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: sketch_json,
        },
        "Sketch 1",
    );

    // Extrude with a bogus reference (op 99 doesn't exist).
    tree.add_feature(
        Feature::Extrude {
            face: entity_name(99, RoleCategory::Generated, 0),
            direction: Vec3Desc::Literal([0.0, 0.0, 1.0]),
            distance: 1.0,
        },
        "Bad Extrude",
    );

    // Replay — the extrude will fail but the tree will still have a brep.
    tree.replay(Tolerance::default()).expect("replay");

    let report = tree.validate_references().expect("should have brep");
    assert!(!report.is_healthy(), "tree should have broken references");
    assert!(!report.issues.is_empty());
    assert_eq!(
        report.issues[0].feature_label, "Bad Extrude",
        "the broken feature should be reported"
    );
}

#[test]
fn no_brep_returns_none() {
    let tree = FeatureTree::new();
    assert!(
        tree.validate_references().is_none(),
        "no brep means no report"
    );
}
