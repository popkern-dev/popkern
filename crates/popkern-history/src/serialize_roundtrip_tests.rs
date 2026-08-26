// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};

use crate::FeatureTree;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc};
use crate::parameter::ParameterId;
use crate::unit::Unit;

/// Build an `EntityName` referencing the generated face at `op_id`.
fn face_name(op_id: u64) -> EntityName {
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

// ---------------------------------------------------------------------------
// Test 7: Round-trip for all feature types
// ---------------------------------------------------------------------------

#[test]
fn save_load_all_feature_types() {
    use crate::feature::{AxisDesc, BoolOpKind, RevolveEnd, ShellDirectionKind, Vec3Desc};
    use crate::node::BodyId;

    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: "base".into(),
        },
        "Sketch",
    );

    let dist = add_param(&mut tree, "_d", 10.0, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );

    let r = add_param(&mut tree, "_r", 2.0, Unit::Millimeters);
    tree.add_feature(
        Feature::Fillet {
            edges: vec![face_name(2)],
            radius: r,
        },
        "Fillet",
    );

    let cd = add_param(&mut tree, "_cd", 1.5, Unit::Millimeters);
    tree.add_feature(
        Feature::Chamfer {
            edges: vec![face_name(2)],
            distance: cd,
        },
        "Chamfer",
    );

    let st = add_param(&mut tree, "_st", 0.5, Unit::Millimeters);
    tree.add_feature(
        Feature::Shell {
            faces_to_remove: vec![face_name(2)],
            thickness: st,
            direction: ShellDirectionKind::Outward,
            join_style: crate::feature::JoinStyleKind::Intersection,
            trim_policy: crate::feature::TrimPolicyKind::Decline,
        },
        "Shell",
    );

    let tt = add_param(&mut tree, "_tt", 0.2, Unit::Millimeters);
    tree.add_feature(
        Feature::Thicken {
            thickness: tt,
            join_style: crate::feature::JoinStyleKind::default(),
            trim_policy: crate::feature::TrimPolicyKind::default(),
        },
        "Thicken",
    );
    let od = add_param(&mut tree, "_od", 0.1, Unit::Millimeters);
    tree.add_feature(
        Feature::OffsetSolid {
            distance: od,
            join_style: crate::feature::JoinStyleKind::default(),
            trim_policy: crate::feature::TrimPolicyKind::default(),
        },
        "OffsetSolid",
    );

    let da = add_param(&mut tree, "_da", 0.1, Unit::Radians);
    tree.add_feature(
        Feature::Draft {
            faces: vec![face_name(2)],
            angle: da,
            neutral_plane: face_name(2),
        },
        "Draft",
    );

    tree.add_feature(
        Feature::Boolean {
            operation: BoolOpKind::Subtraction,
            tool_body: BodyId(1),
        },
        "Boolean",
    );

    let ra = add_param(&mut tree, "_ra", 3.5, Unit::Radians);
    tree.add_feature(
        Feature::Revolve {
            face: face_name(2),
            axis: AxisDesc::Literal {
                origin: [0.0; 3],
                direction: [0.0, 1.0, 0.0],
            },
            end_condition: RevolveEnd::Blind { angle: ra },
            mode: OperationMode::Boss,
        },
        "Revolve",
    );

    let tw = add_param(&mut tree, "_tw", 0.0, Unit::Radians);
    tree.add_feature(
        Feature::Sweep {
            face: face_name(2),
            path_edge: face_name(2),
            twist: tw,
            mode: OperationMode::Cut,
        },
        "Sweep",
    );

    tree.add_feature(
        Feature::Loft {
            sections: vec![face_name(2), face_name(3)],
            closed: true,
            mode: OperationMode::NewBody,
        },
        "Loft",
    );

    let pr = add_param(&mut tree, "_pr", 3.0, Unit::Millimeters);
    tree.add_feature(
        Feature::Pipe {
            path_edge: face_name(2),
            radius: pr,
            mode: OperationMode::NewBody,
        },
        "Pipe",
    );

    let bytes = tree.save().expect("save all features");
    let loaded = FeatureTree::load(&bytes).expect("load all features");

    assert_eq!(loaded.timeline().len(), tree.timeline().len());
    let labels = [
        "Sketch",
        "Extrude",
        "Fillet",
        "Chamfer",
        "Shell",
        "Thicken",
        "OffsetSolid",
        "Draft",
        "Boolean",
        "Revolve",
        "Sweep",
        "Loft",
        "Pipe",
    ];
    for (i, &fid) in loaded.timeline().iter().enumerate() {
        assert_eq!(loaded.feature(fid).unwrap().label, labels[i]);
    }
    let ids = loaded.timeline().to_vec();
    let f = |i: usize| &loaded.feature(ids[i]).unwrap().feature;
    assert!(matches!(f(2), Feature::Fillet { .. }));
    assert!(matches!(f(4), Feature::Shell { .. }));
    assert!(matches!(
        f(5),
        Feature::Thicken {
            join_style: crate::feature::JoinStyleKind::Arc,
            ..
        }
    ));
    assert!(matches!(f(6), Feature::OffsetSolid { .. }));
    // The Shell variant's knobs survive the round trip.
    assert!(matches!(
        f(4),
        Feature::Shell {
            join_style: crate::feature::JoinStyleKind::Intersection,
            trim_policy: crate::feature::TrimPolicyKind::Decline,
            ..
        }
    ));
    assert!(matches!(f(8), Feature::Boolean { .. }));
    assert!(matches!(f(9), Feature::Revolve { .. }));
    assert!(matches!(f(11), Feature::Loft { .. }));
    assert!(matches!(f(12), Feature::Pipe { .. }));
}

// ---------------------------------------------------------------------------
// Test 8: Suppressed + config state round-trip
// ---------------------------------------------------------------------------

#[test]
fn save_load_preserves_suppressed_and_config() {
    let mut tree = FeatureTree::new();
    let s0 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: String::new(),
        },
        "S0",
    );
    let _s1 = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XZ,
            sketch_data: String::new(),
        },
        "S1",
    );
    tree.suppress(s0);

    let cfg_id = tree.add_config("Big", "big").expect("add config");
    let pid = tree
        .add_global_parameter("width", 5.0, Unit::Meters)
        .expect("param");
    tree.set_config_override(cfg_id, pid, crate::config::ParamOverride::Literal(20.0))
        .expect("override");
    tree.set_active_config(Some(cfg_id));

    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");

    let ids = loaded.timeline().to_vec();
    assert!(loaded.feature(ids[0]).unwrap().suppressed);
    assert!(!loaded.feature(ids[1]).unwrap().suppressed);
    assert_eq!(loaded.configs().len(), 1);
    let cfg = loaded.config(loaded.configs()[0]).unwrap();
    assert_eq!(cfg.name, "Big");
    assert_eq!(
        loaded.config_overrides(loaded.configs()[0]).unwrap().len(),
        1
    );
    assert!(loaded.active_config().is_some());
}
