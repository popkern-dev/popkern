// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;
use crate::FeatureTree;
use crate::feature::Feature;
use crate::unit::Unit;
use popkern_core::geometry::Point3;
use popkern_core::naming::{EntityRole, OperationId, RoleCategory};

fn name(op: u64, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole {
            category: RoleCategory::Generated,
            index,
        },
        ancestors: Vec::new(),
    }
}

fn param(tree: &mut FeatureTree, label: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(label, value, unit)
        .expect("param")
}

fn rot_z_matrix(angle: f64) -> [[f64; 4]; 4] {
    Transform3::rotation_z(angle).matrix
}

/// Every motion arm survives save/load bit-for-bit (PK-362
/// acceptance: round-trips through `.pkn`).
#[test]
fn every_motion_round_trips_through_pkn() {
    let mut tree = FeatureTree::new();
    let dx = param(&mut tree, "dx", 1.0, Unit::Meters);
    let dy = param(&mut tree, "dy", 2.0, Unit::Meters);
    let dz = param(&mut tree, "dz", 3.0, Unit::Meters);
    let ang = param(&mut tree, "ang", 0.7, Unit::Radians);
    let motions = [
        RigidMotionDesc::Translate { dx, dy, dz },
        RigidMotionDesc::Rotate {
            axis: AxisDesc::Edge(name(1, 4)),
            angle: ang,
        },
        RigidMotionDesc::Rotate {
            axis: AxisDesc::Literal {
                origin: [0.5, 0.5, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            angle: ang,
        },
        RigidMotionDesc::Matrix(rot_z_matrix(0.3)),
    ];
    for (i, m) in motions.iter().enumerate() {
        tree.add_feature(Feature::Transform { motion: m.clone() }, &format!("T{i}"));
    }
    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");
    for (i, m) in motions.iter().enumerate() {
        match &loaded.feature(loaded.timeline()[i]).expect("node").feature {
            Feature::Transform { motion } => assert_eq!(motion, m),
            other => panic!("expected transform, got {other:?}"),
        }
    }
}

/// Motion parameters and the axis edge participate in dependency
/// reporting exactly like every other feature's parameters/references.
#[test]
fn motion_parameters_and_axis_edge_are_reported() {
    let mut tree = FeatureTree::new();
    let dx = param(&mut tree, "dx", 1.0, Unit::Meters);
    let dy = param(&mut tree, "dy", 2.0, Unit::Meters);
    let dz = param(&mut tree, "dz", 3.0, Unit::Meters);
    let ang = param(&mut tree, "ang", 0.7, Unit::Radians);
    let t = RigidMotionDesc::Translate { dx, dy, dz };
    assert_eq!(t.param_ids(), vec![("dx", dx), ("dy", dy), ("dz", dz)]);
    assert!(t.entity_names().is_empty());
    let r = RigidMotionDesc::Rotate {
        axis: AxisDesc::Edge(name(3, 1)),
        angle: ang,
    };
    assert_eq!(r.param_ids(), vec![("angle", ang)]);
    assert_eq!(r.entity_names(), vec![&name(3, 1)]);
    let m = RigidMotionDesc::Matrix(rot_z_matrix(0.1));
    assert!(m.param_ids().is_empty());
    assert!(m.entity_names().is_empty());
}

/// Resolution reproduces the kernel constructors, and declines the
/// two illegal inputs transactionally: a non-similarity matrix and a
/// degenerate literal axis.
#[test]
fn resolve_matches_kernel_constructors_and_declines_illegal_motions() {
    let mut tree = FeatureTree::new();
    let dx = param(&mut tree, "dx", 1.0, Unit::Meters);
    let dy = param(&mut tree, "dy", -2.0, Unit::Meters);
    let dz = param(&mut tree, "dz", 0.5, Unit::Meters);
    let ang = param(&mut tree, "ang", 0.7, Unit::Radians);
    let brep = BRep::new(popkern_core::tolerance::Tolerance::default());
    let table = tree.parameters();

    let t = RigidMotionDesc::Translate { dx, dy, dz }
        .resolve(&brep, table)
        .expect("translate resolves");
    assert_eq!(t, Transform3::translation(1.0, -2.0, 0.5));

    let r = RigidMotionDesc::Rotate {
        axis: AxisDesc::Literal {
            origin: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 2.0],
        },
        angle: ang,
    }
    .resolve(&brep, table)
    .expect("rotate resolves");
    assert!(r.approx_eq(&Transform3::rotation_z(0.7), 1e-12));

    let about = RigidMotionDesc::Rotate {
        axis: AxisDesc::Literal {
            origin: [1.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
        angle: ang,
    }
    .resolve(&brep, table)
    .expect("rotate about offset axis resolves");
    let fixed = about.apply_point(&Point3::new(1.0, 0.0, 0.0));
    assert!((fixed - Point3::new(1.0, 0.0, 0.0)).length() < 1e-12);

    let shear = RigidMotionDesc::Matrix([
        [1.0, 0.5, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
    .resolve(&brep, table);
    let err = shear.expect_err("shear must decline").to_string();
    assert!(err.contains("similarity"), "{err}");

    let degenerate = RigidMotionDesc::Rotate {
        axis: AxisDesc::Literal {
            origin: [0.0; 3],
            direction: [0.0; 3],
        },
        angle: ang,
    }
    .resolve(&brep, table);
    assert!(degenerate.is_err());
}

/// A tree saved before the variant existed contains no `Transform`
/// features and loads unchanged (no envelope bump in phase 1b).
#[test]
fn pre_transform_document_loads_unchanged() {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: crate::feature::SketchPlaneDesc::XY,
            sketch_data: "{}".into(),
        },
        "Sketch1",
    );
    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");
    assert_eq!(loaded.timeline().len(), 1);
}
