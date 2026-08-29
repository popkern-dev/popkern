// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;
use crate::FeatureTree;
use crate::feature::Feature;
use crate::unit::Unit;
use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;

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

/// Both arms survive a `.pkn` save/load bit-for-bit.
#[test]
fn both_arms_round_trip_through_pkn() {
    let mut tree = FeatureTree::new();
    let angle = tree
        .parameters_mut()
        .add_global("draft", 0.1, Unit::Radians)
        .expect("param");
    let planes = [
        NeutralPlaneDesc::Face(name(1, 2)),
        NeutralPlaneDesc::Literal {
            origin: [0.0, 0.0, -0.5],
            normal: [0.0, 0.0, 2.0],
        },
    ];
    for (i, p) in planes.iter().enumerate() {
        tree.add_feature(
            Feature::Draft {
                faces: vec![name(1, 3)],
                angle,
                neutral_plane: p.clone(),
            },
            &format!("D{i}"),
        );
    }
    let loaded = FeatureTree::load(&tree.save().expect("save")).expect("load");
    for (i, p) in planes.iter().enumerate() {
        match &loaded.feature(loaded.timeline()[i]).expect("node").feature {
            Feature::Draft { neutral_plane, .. } => assert_eq!(neutral_plane, p),
            other => panic!("expected Draft, got {other:?}"),
        }
    }
}

/// The face arm serialises as the bare `EntityName` object that
/// version-3 documents carry, so a pre-PK-363 draft loads unchanged.
#[test]
fn face_arm_is_the_bare_entity_name_wire_form() {
    let n = name(7, 1);
    let bare = serde_json::to_value(&n).expect("name json");
    let desc = serde_json::to_value(NeutralPlaneDesc::Face(n.clone())).expect("desc json");
    assert_eq!(bare, desc, "face arm must not add a tag");
    let back: NeutralPlaneDesc = serde_json::from_value(bare).expect("v3 shape loads");
    assert_eq!(back, NeutralPlaneDesc::Face(n));
    let lit: NeutralPlaneDesc =
        serde_json::from_str(r#"{"origin":[1.0,2.0,3.0],"normal":[0.0,1.0,0.0]}"#)
            .expect("literal shape loads");
    assert_eq!(
        lit,
        NeutralPlaneDesc::Literal {
            origin: [1.0, 2.0, 3.0],
            normal: [0.0, 1.0, 0.0],
        }
    );
}

/// Entity reporting: only the face arm names a model entity.
#[test]
fn entity_name_reports_only_the_face_arm() {
    let n = name(1, 2);
    assert_eq!(NeutralPlaneDesc::Face(n.clone()).entity_name(), Some(&n));
    let lit = NeutralPlaneDesc::Literal {
        origin: [0.0; 3],
        normal: [0.0, 0.0, 1.0],
    };
    assert!(lit.entity_name().is_none());
    let mut face = NeutralPlaneDesc::Face(n.clone());
    assert!(face.entity_name_mut().is_some());
    let mut lit = lit;
    assert!(lit.entity_name_mut().is_none());
}

/// Literal resolve normalises the normal and declines a zero normal
/// with a diagnostic naming the input.
#[test]
fn literal_resolve_normalises_and_declines_degenerate_normal() {
    let brep = BRep::new(Tolerance::default());
    let (o, n) = NeutralPlaneDesc::Literal {
        origin: [1.0, 2.0, 3.0],
        normal: [0.0, 0.0, 4.0],
    }
    .resolve(&brep)
    .expect("resolves");
    assert_eq!(o, Point3::new(1.0, 2.0, 3.0));
    assert!((n - Vec3::new(0.0, 0.0, 1.0)).length() < 1e-15);
    let err = NeutralPlaneDesc::Literal {
        origin: [0.0; 3],
        normal: [0.0; 3],
    }
    .resolve(&brep)
    .expect_err("zero normal declines");
    assert!(err.to_string().contains("degenerate"), "{err}");
    let err = NeutralPlaneDesc::Face(name(9, 9))
        .resolve(&brep)
        .expect_err("missing face declines");
    let _ = err.to_string();
}
