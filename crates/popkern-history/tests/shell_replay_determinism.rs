//! PK-328 row 482: history features carry EVERY shell/offset knob and
//! replay is deterministic — the same tree replays to a
//! fingerprint-identical body across repeated replays and across a
//! save/load round trip, on a curved+join fixture (the outward shell
//! of a box grows cylinder strips and sphere corners, so the join
//! machinery and its options are on the replayed path).

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::*;

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

/// FNV-1a over the serialized body — the workspace-standard
/// whole-body fingerprint (identity pin idiom, PK-276 row 457).
fn body_fingerprint(brep: &BRep) -> u64 {
    let bytes = serde_json::to_vec(brep).expect("output BRep must serialize");
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Box extrude → OffsetSolid(+0.15, Arc/Trim) → Thicken(0.1) tree:
/// every new variant and knob on the replayed path.
fn build_offset_family_tree() -> FeatureTree {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: rectangle_sketch_json(),
        },
        "Sketch1",
    );
    let dist = tree
        .parameters_mut()
        .add_global("_extrude_distance", 0.5, Unit::Meters)
        .expect("param");
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
    let grow = tree
        .parameters_mut()
        .add_global("_grow", 0.15, Unit::Meters)
        .expect("param");
    tree.add_feature(
        Feature::OffsetSolid {
            distance: grow,
            join_style: JoinStyleKind::Arc,
            trim_policy: TrimPolicyKind::Trim,
        },
        "Grow1",
    );
    let wall = tree
        .parameters_mut()
        .add_global("_wall", 0.1, Unit::Meters)
        .expect("param");
    tree.add_feature(
        Feature::Thicken {
            thickness: wall,
            join_style: JoinStyleKind::Arc,
            trim_policy: TrimPolicyKind::Trim,
        },
        "Thicken1",
    );
    tree
}

#[test]
fn offset_family_replay_is_fingerprint_deterministic() {
    let mut tree = build_offset_family_tree();
    tree.replay(Tolerance::default()).expect("first replay");
    let first = body_fingerprint(tree.current_brep().expect("brep"));
    // The grown rounded box gained join geometry: strips + corners.
    let faces = tree.current_brep().expect("brep").faces.len();
    assert!(faces > 12, "join faces must be present, got {faces}");

    tree.replay(Tolerance::default()).expect("second replay");
    let second = body_fingerprint(tree.current_brep().expect("brep"));
    assert_eq!(first, second, "repeated replay must be bit-identical");

    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    loaded.replay(Tolerance::default()).expect("loaded replay");
    let third = body_fingerprint(loaded.current_brep().expect("brep"));
    assert_eq!(
        first, third,
        "save/load replay must be fingerprint-identical"
    );
}

#[test]
fn shell_knobs_change_the_replayed_body() {
    // The knobs must be LOAD-BEARING through replay: the same tree
    // with Intersection joins replays to a DIFFERENT body than with
    // Arc joins (sharp grow vs rounded grow of the box).
    let build = |join: JoinStyleKind| {
        let mut tree = FeatureTree::new();
        tree.add_feature(
            Feature::Sketch {
                plane: SketchPlaneDesc::XY,
                sketch_data: rectangle_sketch_json(),
            },
            "Sketch1",
        );
        let dist = tree
            .parameters_mut()
            .add_global("_d", 0.5, Unit::Meters)
            .expect("param");
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
        let grow = tree
            .parameters_mut()
            .add_global("_g", 0.15, Unit::Meters)
            .expect("param");
        tree.add_feature(
            Feature::OffsetSolid {
                distance: grow,
                join_style: join,
                trim_policy: TrimPolicyKind::Trim,
            },
            "Grow1",
        );
        tree.replay(Tolerance::default()).expect("replay");
        let brep = tree.current_brep().expect("brep").clone();
        (body_fingerprint(&brep), brep.faces.len())
    };
    let (arc_fp, arc_faces) = build(JoinStyleKind::Arc);
    let (miter_fp, miter_faces) = build(JoinStyleKind::Intersection);
    assert_ne!(arc_fp, miter_fp, "join style must be load-bearing");
    assert!(
        arc_faces > miter_faces,
        "rounded grow carries strip/corner faces ({arc_faces}) beyond the sharp grow ({miter_faces})"
    );
}
