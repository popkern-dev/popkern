//! PK-661: the session's history timeline — rows in timeline order with
//! script positions and changed bodies, step meshes keyed so untouched
//! bodies are never re-sent, created-face masks (pattern instances
//! included), and rows attributed to project migrations.

use std::path::Path;

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_history::*;
use popkern_live::{HostOptions, Session};

fn tracked() -> Session {
    Session::new(HostOptions {
        track_sources: true,
        ..HostOptions::default()
    })
}

fn gallery(name: &str) -> String {
    std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../web/live/examples/{name}.rhai")),
    )
    .expect("gallery script")
}

#[test]
fn rows_follow_the_timeline_with_lines_and_changed_bodies() {
    let mut session = tracked();
    let envelope = session.run(&gallery("boss_plate"), "boss_plate.rhai");
    assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
    let parts = session.timeline().expect("timeline");
    assert_eq!(parts.len(), 1);
    let steps = &parts[0].steps;
    let labels: Vec<&str> = steps.iter().map(|s| s.label.as_str()).collect();
    let tree_labels: Vec<String> = session
        .model()
        .expect("model")
        .features()
        .into_iter()
        .map(|f| f.label)
        .collect();
    assert_eq!(labels, tree_labels, "timeline order");
    let lines: Vec<u32> = steps
        .iter()
        .map(|s| s.source.as_ref().expect("tracked").line)
        .collect();
    // plate sketch (5), extrude (6), boss sketch (7), extrude (8),
    // union (9), rounds (from line 10).
    assert_eq!(&lines[..5], &[5, 6, 7, 8, 9]);
    assert!(lines[5] >= 10 && lines[5] <= 14, "rounds at {}", lines[5]);
    assert!(steps.iter().all(|s| {
        s.source
            .as_ref()
            .is_some_and(|src| src.script == "boss_plate.rhai")
    }));
    // The boss sketch builds a tool body; the union changes the plate.
    assert_eq!(steps[2].bodies_changed, [1]);
    assert_eq!(steps[4].bodies_changed, [0]);
    assert!(steps.iter().all(|s| !s.failed && !s.suppressed));
}

#[test]
fn untouched_bodies_keep_their_key_and_are_not_resent() {
    let mut session = tracked();
    session.run(&gallery("boss_plate"), "boss_plate.rhai");
    // After the plate extrude (step 2): one body, meshed.
    let plate = session.step_bodies(None, 2, &[]).expect("step 2");
    assert_eq!(plate.len(), 1);
    assert!(plate[0].buffers && plate[0].mesh.is_some());
    let known: Vec<String> = plate.iter().map(|b| b.key.clone()).collect();
    // After the boss sketch (step 3): the plate is untouched — same
    // key, no buffers sent — and the new tool body is meshed.
    let boss = session.step_bodies(None, 3, &known).expect("step 3");
    let plate_row = boss.iter().find(|b| b.body == 0).expect("plate");
    assert_eq!(plate_row.key, known[0]);
    assert!(!plate_row.buffers && plate_row.mesh.is_none());
    assert!(plate_row.created_faces.is_empty());
    let tool_row = boss.iter().find(|b| b.body == 1).expect("tool");
    assert!(tool_row.buffers);
    assert!(!tool_row.created_faces.is_empty(), "the sketch face is new");
    // The rounds step marks exactly the rounds' own faces.
    let last = session.model().expect("model").features().len();
    let rows = session.step_bodies(None, last, &[]).expect("last step");
    let plate = &rows[0];
    let mesh = plate.mesh.as_ref().expect("mesh");
    let rounds: Vec<u32> = mesh
        .face_labels
        .iter()
        .enumerate()
        .filter(|(_, l)| l.starts_with("rounds/"))
        .map(|(i, _)| i as u32)
        .collect();
    assert!(!rounds.is_empty());
    assert_eq!(
        plate.created_faces, rounds,
        "created = faces named under the fillet"
    );
}

#[test]
fn suppressed_rows_are_marked() {
    let mut session = tracked();
    let source = r#"
let m = model();
let b = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(0.5));
b.fillet(edges().straight().direction("z"), 0.1, #{ label: "rounds" });
m.feature("rounds").suppress();
m
"#;
    assert!(session.run(source, "s.rhai").outcome.is_served());
    let steps = &session.timeline().expect("timeline")[0].steps;
    let rounds = steps.iter().find(|s| s.label == "rounds").expect("rounds");
    assert!(rounds.suppressed);
    assert!(rounds.bodies_changed.is_empty());
}

#[test]
fn project_rows_name_their_migration() {
    let dir =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../popkern-project/tests/fixtures/boss_plate");
    let mut files = Vec::new();
    for sub in [
        "popkern.toml",
        "migrations/001_plate.rhai",
        "migrations/002_boss.rhai",
        "migrations/003_rounds.rhai",
    ] {
        files.push((
            sub.to_string(),
            std::fs::read(dir.join(sub)).expect("fixture"),
        ));
    }
    let mut session = tracked();
    assert!(session.run_project(files).outcome.is_served());
    let steps = &session.timeline().expect("timeline")[0].steps;
    let by: Vec<(&str, Option<&str>, Option<&str>)> = steps
        .iter()
        .map(|s| {
            (
                s.label.as_str(),
                s.migration.as_deref(),
                s.source.as_ref().map(|src| src.script.as_str()),
            )
        })
        .collect();
    assert_eq!(
        by[0],
        ("plate_sketch", Some("001_plate"), Some("001_plate.rhai"))
    );
    assert_eq!(by[2].1, Some("002_boss"));
    assert_eq!(by.last().expect("rounds").1, Some("003_rounds"));
    assert_eq!(by.last().expect("rounds").2, Some("003_rounds.rhai"));
}

fn square(x0: f64, x1: f64) -> String {
    let mut sk = Sketch::new();
    let p = [
        sk.add_point(x0, x0),
        sk.add_point(x1, x0),
        sk.add_point(x1, x1),
        sk.add_point(x0, x1),
    ];
    for i in 0..4 {
        sk.add_line(p[i], p[(i + 1) % 4]);
        sk.add_constraint(Constraint::Fixed(p[i]));
    }
    sk.solve().expect("square");
    serde_json::to_string(&sk).expect("json")
}

fn generated(op: u64) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    }
}

/// A box with an on-face boss, linearly patterned three times.
fn patterned_document() -> Vec<u8> {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(square(0.0, 3.0)),
        },
        "base_sketch",
    );
    let d1 = tree
        .parameters_mut()
        .add_global("d1", 1.0, Unit::Meters)
        .expect("p");
    let extrude = |face, distance, mode| Feature::Extrude {
        face,
        direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
        end_condition: ExtrudeEnd::Blind { distance },
        mode,
        draft_angle: None,
        second_direction: None,
    };
    tree.add_feature(extrude(generated(1), d1, OperationMode::NewBody), "base");
    tree.replay(Tolerance::default()).expect("box");
    let brep = tree.current_brep().expect("box").clone();
    let cap = brep
        .faces
        .iter()
        .filter_map(|(_, f)| brep.name_pool.get(f.name?))
        .find(|n| n.operation == OperationId(2) && n.role.category == RoleCategory::Cap)
        .cloned()
        .expect("cap");
    let sk = tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::OnFace(cap),
            sketch_data: SketchPayload::Constraint(square(0.2, 0.5)),
        },
        "boss_sketch",
    );
    tree.replay(Tolerance::default()).expect("sketch");
    let d2 = tree
        .parameters_mut()
        .add_global("d2", 0.3, Unit::Meters)
        .expect("p");
    let boss = tree.add_feature(extrude(generated(3), d2, OperationMode::Boss), "boss");
    tree.replay(Tolerance::default()).expect("boss");
    tree.add_linear_pattern_with_pitch(
        vec![sk, boss],
        Vec3Desc::Literal([1.0, 0.0, 0.0]),
        3.0,
        // Pitch is in display millimetres: 0.8 m apart. (0.8 placed the
        // instances 0.8 mm off the seed — one merged boss — which the
        // per-boolean renaming of every face hid until naming scheme 4.)
        800.0,
        "bosses",
    )
    .expect("pattern");
    tree.replay(Tolerance::default()).expect("pattern replay");
    tree.save().expect("bytes")
}

#[test]
fn a_pattern_step_marks_its_instances_faces() {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("timeline_pattern.pkn");
    std::fs::write(&path, patterned_document()).expect("write");
    let mut session = tracked();
    let script = format!(
        "let m = open({:?});\nm.replay();\nm\n",
        path.display().to_string()
    );
    let envelope = session.run(&script, "open.rhai");
    assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
    let model = session.model_mut().expect("model");
    let n = model.features().len();
    let pattern = model.features().last().expect("pattern").id;
    let node = model.tree().feature(pattern).expect("node").clone();
    let state = model.step_state(n).expect("state");
    let instance_ops: Vec<u64> = node.instance_op_ids.values().map(|o| o.0).collect();
    assert!(
        !instance_ops.is_empty(),
        "the pattern minted instance operations"
    );
    let (body, brep) = state.bodies().next().expect("body");
    let named = |ops: &[u64]| -> usize {
        brep.faces
            .iter()
            .filter(|(_, f)| {
                f.name
                    .and_then(|id| brep.name_pool.get(id))
                    .is_some_and(|name| ops.contains(&name.operation.0))
            })
            .count()
    };
    let rows = session.step_bodies(None, n, &[]).expect("pattern step");
    let row = rows
        .iter()
        .find(|r| r.body == u64::from(body.0))
        .expect("row");
    let mut ops = vec![node.operation_id.0];
    ops.extend(&instance_ops);
    assert_eq!(row.created_faces.len(), named(&ops));
    // The merge names every instance's faces under the pattern's own
    // operation: the two new bosses (five faces each) are all created
    // by this step, plus the cap they modified.
    assert!(
        row.created_faces.len() >= 10,
        "both instances' faces are marked: {}",
        row.created_faces.len()
    );
}
