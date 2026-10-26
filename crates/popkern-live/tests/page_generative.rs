//! PK-956: a generative design in popkern live, end to end.
//!
//! Natively a session solves the design as it builds, and the designs
//! its run used export as a bake carries them. Under node the worker
//! (`web/live/test/generative.test.mjs`) declines the same build
//! unsolved (`NeedsSolve`), serves it once seeded with the native
//! export (the wasm cache key is the native one, bit for bit, and the
//! body the same), re-solves a variant on request with a progress
//! message per update, and serves that from its cache on the next
//! plain build. The design's history (PK-959) travels as the bake's
//! design file: once seeded, the worker meshes the design among its
//! records exactly as the native session does, and a film evolves it.
//! Needs the web bundle (`tools/build-live.sh --web`) and node; the
//! native leg runs regardless.

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_live::{HostOptions, Session};
use serde_json::Value;

/// A 12 × 4 × 4 mm cantilever, clamped at x = 0 and pulled down at
/// x = 12, carved to `fraction` of its volume (small enough to solve
/// in seconds).
fn script(fraction: f64) -> String {
    format!(
        r#"let m = model();
let b = m.profile(yz()).rect(4.0, 4.0).sketch("block").extrude(blind(12.0));
m.label_face("wall", faces().facing_neg("x"));
m.label_face("tip", faces().facing("x"));
b = b.generative(#{{
    supports: [#{{ faces: faces().labeled("wall") }}],
    loads: [#{{ faces: faces().labeled("tip"), force: 10.0, direction: [0.0, 0.0, -1.0] }}],
    keep: [faces().labeled("wall"), faces().labeled("tip")],
    volume_fraction: {fraction:?},
    filter: 1.5,
    keep_depth: 1.0,
    cell: 1.0,
    max_iterations: 150,
    label: "carve",
}});
m
"#
    )
}

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

/// The generative row of an envelope's report.
fn row(envelope: &Value) -> Value {
    envelope["report"]["features"]
        .as_array()
        .and_then(|f| f.iter().find(|r| r["operation"] == "generative"))
        .cloned()
        .unwrap_or(Value::Null)
}

#[test]
fn a_bake_carries_its_design_to_the_browser_which_otherwise_solves_only_on_request() {
    // Native: solve as it builds; the run's design exports.
    let source = script(0.5);
    let _ = popkern_generative::export_touched();
    let mut session = Session::new(HostOptions::default());
    let native = session.run(&source, "generative.rhai").to_json();
    assert_eq!(native["outcome"], "served", "{native}");
    let native_row = row(&native);
    assert_eq!(native_row["generative"]["cached"], false, "{native_row}");
    let touched = popkern_generative::touched_designs();
    let designs = popkern_generative::export_touched().expect("the run's design");
    let histories = popkern_generative::export_histories(&touched).expect("its history");
    let design = native_row["generative"]["design"]
        .as_str()
        .unwrap_or("")
        .to_string();
    assert_eq!(touched, vec![design.clone()]);
    // Natively, a film evolves the design the build solved.
    let film = session
        .film(
            "let film = storyboard();\nfilm.step(\"Carve\", 4.0, |s| { s.evolve(\"carve\"); });\nfilm",
            "films/carve.rhai",
            None,
        )
        .unwrap_or_else(|outcome| panic!("{outcome:?}"));
    assert_eq!(film.context.designs.len(), 1);
    assert_eq!(film.context.designs[0].design, design);
    let parsed: Value = serde_json::from_str(&designs).expect("JSON");
    assert_eq!(parsed["entries"].as_array().map(Vec::len), Some(1));
    // A second native build is served from the cache.
    let again = session.run(&source, "generative.rhai").to_json();
    assert_eq!(row(&again)["generative"]["cached"], true);
    assert_eq!(row(&again)["generative"], {
        let mut r = native_row["generative"].clone();
        r["cached"] = Value::Bool(true);
        r
    });

    let page = root().join("web/live");
    if !page.join("pkg/popkern_live.js").exists() {
        eprintln!("[pk956] SKIP browser leg: no web bundle (run tools/build-live.sh --web)");
        return;
    }
    if Command::new("node").arg("--version").output().is_err() {
        eprintln!("[pk956] SKIP browser leg: node is not installed");
        return;
    }
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("pk956");
    std::fs::create_dir_all(&dir).expect("tmp");
    let (script_path, designs_path, variant_path, histories_path) = (
        dir.join("generative.rhai"),
        dir.join("designs.json"),
        dir.join("variant.rhai"),
        dir.join("designs.bin"),
    );
    std::fs::write(&histories_path, &histories).expect("write");
    std::fs::write(&script_path, &source).expect("write");
    std::fs::write(&designs_path, &designs).expect("write");
    std::fs::write(&variant_path, script(0.45)).expect("write");
    let output = Command::new("node")
        .arg(page.join("test/generative.test.mjs"))
        .args([&script_path, &designs_path, &variant_path, &histories_path])
        .output()
        .expect("node");
    assert!(
        output.status.success(),
        "generative.test.mjs failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let out: Value = serde_json::from_slice(&output.stdout).expect("harness JSON");
    // Unsolved, the browser declines and asks for a re-solve.
    assert_eq!(out["deferred"]["outcome"], "kernel_decline", "{out}");
    assert_eq!(out["deferred"]["generative_kind"], "NeedsSolve", "{out}");
    // Seeded with the native design, it serves it unsolved: same key,
    // same body.
    assert_eq!(out["seed"]["added"], 1, "{out}");
    assert_eq!(out["seeded"]["outcome"], "served", "{out}");
    let seeded = &out["seeded"]["generative"];
    assert_eq!(seeded["cached"], true, "{out}");
    assert_eq!(seeded["faces"], native_row["generative"]["faces"], "{out}");
    assert_eq!(
        seeded["compliance_nmm"], native_row["generative"]["compliance_nmm"],
        "{out}"
    );
    // A re-solve of a design nobody solved: served, one progress
    // message per update, all for that request.
    assert_eq!(out["solve"]["outcome"], "served", "{out}");
    let solved = &out["solve"]["generative"];
    assert_eq!(solved["cached"], false, "{out}");
    assert_eq!(out["progress"]["count"], solved["iterations"], "{out}");
    assert_eq!(out["progress"]["ids"].as_array().map(Vec::len), Some(1));
    assert_eq!(out["progress"]["last"]["count"], solved["iterations"]);
    // And the next plain build serves it from the worker's cache.
    assert_eq!(out["cached"]["outcome"], "served", "{out}");
    assert_eq!(out["cached"]["generative"]["cached"], true, "{out}");

    // The design converging (PK-959): a seeded design's history is not
    // in the worker until the bake's design file is; a re-solved one's
    // is. The worker's meshes are the native session's, to the bit
    // (the same f32 positions sum the same in the same order).
    assert_eq!(out["history_missing"], false, "{out}");
    assert_eq!(out["history_seed"]["added"], 1, "{out}");
    assert_eq!(out["variant_history"], true, "{out}");
    let frames = out["history"]["frames"].as_u64().unwrap_or(0) as usize;
    assert!(frames >= 10, "{out}");
    assert_eq!(
        out["history"]["updates"],
        native_row["generative"]["iterations"]
    );
    let native_mesh = |at: Option<f64>| {
        let mesh = popkern_live::designs::mesh(&design, at).expect("native history");
        let sum: f64 = mesh.positions.iter().map(|&v| f64::from(v)).sum();
        (mesh.indices.len() / 3, sum)
    };
    let last = (frames - 1) as f64;
    for (name, at) in [
        ("first", Some(0.0)),
        ("middle", Some(last / 2.0)),
        ("last", Some(last)),
        ("space", None),
    ] {
        let (triangles, sum) = native_mesh(at);
        let wasm = &out["meshes"][name];
        assert_eq!(
            wasm["triangles"].as_u64(),
            Some(triangles as u64),
            "{name}: {out}"
        );
        assert_eq!(wasm["sum"].as_f64(), Some(sum), "{name}: {out}");
    }
    // The design grows from the kept faces to the final design.
    let triangles = |name: &str| out["meshes"][name]["triangles"].as_u64().unwrap_or(0);
    assert!(triangles("last") > 0 && triangles("space") > 0, "{out}");
    assert_eq!(
        out["film"]["designs"][0]["design"],
        design.as_str(),
        "{out}"
    );
    assert_eq!(out["film"]["evolve"][0]["target"], "design carve", "{out}");
}

/// The stable names of the faces lying on the plane x = `x` (within
/// 0.01 mm: a design solved on the fallback grid keeps its faces within
/// a thousandth of a cell), largest first (the display mesh's face
/// labels, as the Simulate panel's picks read them): on a design body,
/// the kept face's patches.
fn faces_at_x(envelope: &popkern_live::Envelope, x: f64) -> Vec<String> {
    let body = &envelope.bodies[0];
    let mut count = vec![0usize; body.face_count];
    let mut off = vec![false; body.face_count];
    for (t, &face) in body.face_ids.iter().enumerate() {
        count[face as usize] += 1;
        for k in 0..3 {
            let v = body.indices[3 * t + k] as usize;
            if (f64::from(body.positions[3 * v]) - x).abs() > 0.01 {
                off[face as usize] = true;
            }
        }
    }
    let mut on: Vec<usize> = (0..body.face_count)
        .filter(|&f| !off[f] && count[f] > 0)
        .collect();
    on.sort_by_key(|&f| std::cmp::Reverse(count[f]));
    assert!(!on.is_empty(), "no face on x = {x}");
    on.iter().map(|&f| body.face_labels[f].clone()).collect()
}

/// A cantilever setup holding the faces `held` and pulling the faces
/// `loaded` down with 10 N, every face by its stable name.
fn setup(held: &[String], loaded: &[String]) -> String {
    let mut faces = String::new();
    let mut list = |prefix: &str, names: &[String]| {
        let keys: Vec<String> = (0..names.len())
            .map(|i| format!("\"{prefix}_{i}\""))
            .collect();
        for (i, n) in names.iter().enumerate() {
            faces += &format!("{prefix}_{i} = {{ name = \"{n}\" }}\n");
        }
        keys.join(", ")
    };
    let (h, l) = (list("held", held), list("tip", loaded));
    format!(
        "units = \"mm\"\nmaterial = \"PLA (FDM)\"\n\n[faces]\n{faces}\n[[support]]\nfaces = [{h}]\nkind = \"fixed\"\n\n[[load]]\nfaces = [{l}]\nforce_n = 10\ndirection = [0, 0, -1]\n"
    )
}

/// Simulate on a design (PK-956 acceptance): a setup naming the kept
/// wall and tip by their stable names solves on the carved body; the
/// names survive a re-solve at another volume fraction (the kept faces'
/// names come from the design space's, not from the design); the body
/// converges under refinement (it is the design the solver found, no
/// pinched hinge); less material bends more, and the design bends more
/// than the solid it was carved from, within a few times.
#[test]
fn simulate_solves_on_the_design_by_the_kept_faces_names() {
    let mut session = Session::new(HostOptions::default());
    let carved = session.run(&script(0.5), "generative.rhai");
    assert!(carved.outcome.is_served(), "{:?}", carved.outcome);
    let (wall, tip) = (faces_at_x(&carved, 0.0), faces_at_x(&carved, 12.0));
    assert!(
        wall.iter().chain(&tip).all(|n| n != "unnamed"),
        "{wall:?} {tip:?}"
    );
    let solve = |session: &mut Session, held: &[String], loaded: &[String]| {
        let out = session.simulate(&setup(held, loaded)).expect("simulate");
        assert_eq!(out.summary["status"], "solved", "{:#}", out.summary);
        out.summary["result"]["max_displacement_mm"]
            .as_f64()
            .expect("deflection")
    };
    let design = solve(&mut session, &wall, &tip);
    // The body is what the design solved: resolving it finer converges
    // (a one-node hinge left as a pinch read 0.25, 0.36, 3.27 mm at 1,
    // 0.5 and 0.25 mm cells; closed, 0.22, 0.245, 0.256).
    let at = |session: &mut Session, cell: f64| {
        let text = format!(
            "{}\n[options]\ntarget_cell_size = {cell}\n",
            setup(&wall, &tip)
        );
        let out = session.simulate(&text).expect("simulate");
        out.summary["result"]["max_displacement_mm"]
            .as_f64()
            .expect("deflection")
    };
    let (half, quarter) = (at(&mut session, 0.5), at(&mut session, 0.25));
    assert!(
        (quarter - half).abs() < 0.1 * half,
        "the body converges under refinement: {half} mm at 0.5, {quarter} at 0.25"
    );
    let other = session.run(&script(0.45), "generative.rhai");
    assert!(other.outcome.is_served(), "{:?}", other.outcome);
    // The largest patch of each kept face keeps its name across
    // designs (the others are as many as the design makes).
    assert_eq!(
        (
            faces_at_x(&other, 0.0)[0].clone(),
            faces_at_x(&other, 12.0)[0].clone()
        ),
        (wall[0].clone(), tip[0].clone()),
        "kept faces keep their names across designs"
    );
    let (wall, tip) = (faces_at_x(&other, 0.0), faces_at_x(&other, 12.0));
    let thinner = solve(&mut session, &wall, &tip);
    assert!(
        thinner > design,
        "less material bends more: {thinner} vs {design}"
    );
    // The solid block (the design space alone).
    let solid_source = script(0.5)
        .split("b = b.generative")
        .next()
        .map(|head| format!("{head}m\n"))
        .expect("script head");
    let solid = session.run(&solid_source, "solid.rhai");
    assert!(solid.outcome.is_served(), "{:?}", solid.outcome);
    let block = {
        let (w, t) = (faces_at_x(&solid, 0.0), faces_at_x(&solid, 12.0));
        solve(&mut session, &w, &t)
    };
    eprintln!(
        "[pk956] tip deflection: solid {block:.4} mm, design (50%) {design:.4} mm, design (45%) {thinner:.4} mm"
    );
    assert!(
        design > block && design < 4.0 * block,
        "{design} vs solid {block}"
    );
}

/// The stable names of the faces whose every vertex lies within 0.8 mm
/// of the cylinder of `radius` about the y-parallel axis through
/// (`x`, `z`): the eye bore's patches (a faceted bore's corners sit
/// where the cells' crossings fall, within half a 1.25 mm cell).
fn faces_on_eye(envelope: &popkern_live::Envelope, x: f64, z: f64, radius: f64) -> Vec<String> {
    let body = &envelope.bodies[0];
    let mut off = vec![false; body.face_count];
    let mut seen = vec![false; body.face_count];
    for (t, &face) in body.face_ids.iter().enumerate() {
        seen[face as usize] = true;
        for k in 0..3 {
            let v = body.indices[3 * t + k] as usize;
            let (px, pz) = (
                f64::from(body.positions[3 * v]),
                f64::from(body.positions[3 * v + 2]),
            );
            let r = ((px - x).powi(2) + (pz - z).powi(2)).sqrt();
            if (r - radius).abs() > 0.8 {
                off[face as usize] = true;
            }
        }
    }
    (0..body.face_count)
        .filter(|&f| seen[f] && !off[f])
        .map(|f| body.face_labels[f].clone())
        .collect()
}

/// The gallery bracket's strength (PK-956 acceptance, minutes: run with
/// `--ignored`): held on its wall and loaded down through its eye, the
/// design (a quarter of the material) against the design space block.
#[test]
#[ignore = "heavy: solves the gallery bracket natively (minutes); run with --ignored for the example's Simulate figures"]
fn the_gallery_bracket_against_its_design_space() {
    let source = std::fs::read_to_string(
        root().join("crates/popkern-rhai/examples/generative_bracket.rhai"),
    )
    .expect("example");
    let block_source = source
        .split("let bracket = space.generative")
        .next()
        .map(|head| format!("{head}m\n"))
        .expect("script head");
    let mut session = Session::new(HostOptions::default());
    let mut measure = |source: &str| {
        let envelope = session.run(source, "bracket.rhai");
        assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
        let wall = faces_at_x(&envelope, 0.0);
        let eye = faces_on_eye(&envelope, 44.0, 0.0, 3.0);
        assert!(!eye.is_empty(), "no eye faces");
        let text = setup(&wall, &eye).replace("force_n = 10", "force_n = 200");
        let out = session.simulate(&text).expect("simulate");
        assert_eq!(out.summary["status"], "solved", "{:#}", out.summary);
        out.summary["result"].clone()
    };
    let block = measure(&block_source);
    let design = measure(&source);
    eprintln!("[pk956] block: {block}");
    eprintln!("[pk956] design: {design}");
}
