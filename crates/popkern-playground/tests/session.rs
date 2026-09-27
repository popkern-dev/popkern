//! PK-410: the playground engine natively — envelope shapes, the one
//! class partition shared with `popkern run`, budgets, parameter
//! edits, and the byte-form exports.

use std::path::{Path, PathBuf};
use std::time::Duration;

use popkern_playground::mesh_buffers::{fnv1a_f32, welded_open_edges};
use popkern_playground::{HostOptions, Outcome, RunOptions, Session, WallBudget};
use popkern_script::Model;

fn boss_plate() -> String {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../popkern-rhai/examples/boss_plate.rhai");
    std::fs::read_to_string(path).expect("boss plate example")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("popkern_playground");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

#[test]
fn boss_plate_serves_with_one_watertight_body_identical_to_the_native_mesh() {
    let mut session = Session::new(HostOptions::default());
    let envelope = session.run(&boss_plate(), "boss_plate.rhai");
    assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
    assert_eq!(envelope.features, 6);
    assert_eq!(envelope.printed.len(), 1);
    assert!(envelope.clock_available);
    assert_eq!(envelope.wall_budget, WallBudget::None);
    assert_eq!(envelope.bodies.len(), 1, "one alive body");
    let body = &envelope.bodies[0];
    assert_eq!(body.open_edges, 0, "watertight");
    assert_eq!(body.positions.len(), body.vertex_count * 3);
    assert_eq!(body.indices.len(), body.triangle_count * 3);
    assert_eq!(body.face_ids.len(), body.triangle_count);
    assert_eq!(body.face_labels.len(), body.face_count);
    assert!(
        body.face_labels.iter().all(|l| l != "unnamed"),
        "every face carries a stable name: {:?}",
        body.face_labels
    );
    assert!(
        body.face_ids
            .iter()
            .all(|&f| (f as usize) < body.face_count)
    );
    assert!(body.bounds[0][2] < body.bounds[1][2]);
    // The same document tessellated through the recorder directly, at
    // the display quality the page receives.
    let model = session.model_mut().expect("model");
    let brep = model.to_brep().expect("brep");
    let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    assert_eq!(mesh.triangles.len(), body.triangle_count);
    assert_eq!(welded_open_edges(&mesh), 0);
    let native: Vec<f32> = mesh
        .vertices
        .iter()
        .flat_map(|v| {
            [
                v.position.x as f32,
                v.position.y as f32,
                v.position.z as f32,
            ]
        })
        .collect();
    assert_eq!(fnv1a_f32(&native), fnv1a_f32(&body.positions));
    // JSON form: outcome flattened, buffers omitted, report counts present.
    let json = envelope.to_json();
    assert_eq!(json["outcome"], "served");
    assert_eq!(json["report"]["failures"], 0);
    assert_eq!(json["bodies"][0]["open_edges"], 0);
    assert!(json["bodies"][0].get("positions").is_none());
    assert!(json["params"].as_array().is_some_and(|p| !p.is_empty()));
}

#[test]
fn failures_partition_exactly_as_the_cli_does() {
    let mut session = Session::new(HostOptions::default());
    let syntax = session.run("let m = model(\n", "syntax.rhai");
    assert!(
        matches!(&syntax.outcome, Outcome::ScriptError { kind, .. } if kind == "SyntaxError"),
        "{:?}",
        syntax.outcome
    );
    assert_eq!(syntax.outcome.exit_code(), 2);
    assert!(syntax.bodies.is_empty());
    let decline = session.run(
        r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges(), 0.9);
m"#,
        "decline.rhai",
    );
    match &decline.outcome {
        Outcome::KernelDecline { kind, payload, .. } => {
            assert_eq!(kind, "FeatureFailed");
            assert_eq!(payload["error"]["kind"], "BlendDeclined");
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(decline.outcome.exit_code(), 3);
    let none = session.run("let x = 1;\n", "none.rhai");
    assert!(
        matches!(&none.outcome, Outcome::GateFailure { step, .. } if step == "model"),
        "{:?}",
        none.outcome
    );
    assert_eq!(none.outcome.exit_code(), 4);
    // No model yet: exports and edits say so instead of panicking.
    assert!(matches!(
        session.export_pkn(),
        Err(Outcome::GateFailure { .. })
    ));
    assert!(matches!(
        session.to_rhai(),
        Err(Outcome::GateFailure { .. })
    ));
    let edit = session.set_param("nope", 1.0);
    assert!(matches!(edit.outcome, Outcome::GateFailure { .. }));
}

#[test]
fn a_failed_run_keeps_the_previous_model() {
    let mut session = Session::new(HostOptions::default());
    assert!(session.run(&boss_plate(), "boss.rhai").outcome.is_served());
    let broken = session.run("let m = model(\n", "broken.rhai");
    assert!(!broken.outcome.is_served());
    assert_eq!(session.model().map(|m| m.features().len()), Some(6));
    assert!(session.export_pkn().is_ok(), "exports still answer");
}

#[test]
fn budgets_stop_runaway_scripts_with_the_budget_named() {
    let ops = RunOptions {
        max_operations: Some(200_000),
        timeout_ms: None,
        ..RunOptions::default()
    }
    .host();
    let mut session = Session::new(ops);
    let envelope = session.run(
        "let m = model();\nlet i = 0;\nloop { i += 1; }\n",
        "spin.rhai",
    );
    match &envelope.outcome {
        Outcome::ScriptError { kind, payload, .. } => {
            assert_eq!(kind, "BudgetExceeded");
            assert_eq!(payload["budget"], "operations");
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(envelope.wall_budget, WallBudget::None);
    let wall = RunOptions {
        max_operations: None,
        timeout_ms: Some(300.0),
        ..RunOptions::default()
    }
    .host();
    assert_eq!(wall.timeout, Some(Duration::from_millis(300)));
    let mut session = Session::new(wall);
    let envelope = session.run(
        r#"let m = model();
let i = 0;
loop { m.profile(xy_at(i * 1.0)).rect(1.0, 1.0).sketch("s" + i).extrude(blind(0.5)); i += 1; }
"#,
        "pattern.rhai",
    );
    match &envelope.outcome {
        Outcome::ScriptError { kind, payload, .. } => {
            assert_eq!(kind, "BudgetExceeded");
            assert_eq!(payload["budget"], "wall_clock");
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(
        envelope.wall_budget,
        WallBudget::Armed,
        "native clock present"
    );
    // Without a clock the same request is reported as unenforceable,
    // never silently ignored.
    assert_eq!(WallBudget::describe(true, false), WallBudget::Unavailable);
    assert_eq!(WallBudget::describe(false, false), WallBudget::None);
    // Malformed budgets are dropped, not misread.
    let bad = RunOptions {
        max_operations: None,
        timeout_ms: Some(f64::NAN),
        ..RunOptions::default()
    }
    .host();
    assert_eq!(bad.timeout, None);
}

#[test]
fn set_param_replays_and_reports_bounds_and_unknown_names() {
    let mut session = Session::new(HostOptions::default());
    let first = session.run(&boss_plate(), "boss.rhai");
    let before = first.bodies[0].triangle_count;
    // The boss (z 0.2..1.0) owns the top of the bounds; the plate
    // growing to 0.5 shows in the mesh and the parameter row instead.
    let edited = session.set_param("extrude-1.distance", 0.5);
    assert!(edited.outcome.is_served(), "{:?}", edited.outcome);
    assert_eq!(edited.bodies[0].open_edges, 0, "still watertight");
    assert!(edited.bodies[0].triangle_count > 0 && before > 0);
    let model = session.model_mut().expect("model");
    let brep = model.to_brep().expect("brep");
    let top_of_plate = brep
        .vertices
        .values()
        .map(|v| v.position.z)
        .filter(|z| (*z - 0.5).abs() < 1e-9)
        .count();
    assert!(top_of_plate >= 4, "plate top moved to z=0.5");
    let row = edited
        .params
        .iter()
        .find(|p| p.name == "extrude-1.distance")
        .expect("param row");
    assert!((row.value - 0.5).abs() < 1e-12);
    let unknown = session.set_param("nope", 1.0);
    match &unknown.outcome {
        Outcome::ScriptError { kind, .. } => assert_eq!(kind, "ArgumentError"),
        other => panic!("{other:?}"),
    }
    assert_eq!(unknown.features, 6, "the model is still described");
}

#[test]
fn byte_exports_match_the_path_forms_and_reopen() {
    let mut session = Session::new(HostOptions::default());
    assert!(session.run(&boss_plate(), "boss.rhai").outcome.is_served());
    let step = session.export_step("boss").expect("step bytes");
    let pkn = session.export_pkn().expect("pkn bytes");
    let model = session.model_mut().expect("model");
    let step_path = tmp("boss.step");
    let pkn_path = tmp("boss.pkn");
    model.to_step(&step_path).expect("to_step");
    model.save(&pkn_path).expect("save");
    assert_eq!(std::fs::read(&step_path).expect("step"), step, "STEP bytes");
    assert_eq!(std::fs::read(&pkn_path).expect("pkn"), pkn, ".pkn bytes");
    let mut reopened = Model::from_pkn_bytes(&pkn).expect("reopen");
    assert_eq!(reopened.replay().expect("replay").failures().count(), 0);
    assert_eq!(reopened.features().len(), 6);
    let code = session.to_rhai().expect("rhai");
    assert!(code.contains("fillet"), "{code}");
    // OCC oracle on the bytes when the venv is present.
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_check.py");
    if oracle.exists() && script.exists() {
        let output = std::process::Command::new(&oracle)
            .arg(&script)
            .arg(&step_path)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(report.contains("\"passed\": true"), "OCC oracle: {report}");
    } else {
        eprintln!("[pk410] OCC oracle absent; STEP bytes checked by popkern only");
    }
}
