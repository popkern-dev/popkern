//! PK-410: the wasm bundle answers as the native engine — the node
//! smoke drives the wasm-bindgen module the way the page's worker will
//! and this side replays the same scripts natively, comparing classes,
//! counts, mesh topology (bit for bit), geometry (value for value: the
//! two platforms' libm differ in the last ulp of transcendental
//! functions, measured 1.5e-16 on the boss plate), the document bytes
//! and the canonical script (bit for bit).
//!
//! Needs the nodejs bundle from `tools/build-live.sh --node`
//! (CI builds it first). Set `POPKERN_LIVE_BUILD=1` to have this
//! test run the build itself. Without the bundle the test records why
//! it skipped and passes — the build rot gate is the CI job, which
//! never skips.

use std::path::{Path, PathBuf};
use std::process::Command;

#[path = "support/node_bundle.rs"]
mod node_bundle;

use node_bundle::{bundle, root};
use popkern_live::mesh_buffers::{fnv1a_f32, fnv1a_u32};
use popkern_live::{HostOptions, Session};
use serde_json::Value;

fn smoke(pkg: &Path, script: &Path, max_ops: Option<u64>, timeout_ms: Option<f64>) -> Value {
    smoke_dumping(pkg, script, max_ops, timeout_ms, None)
}

fn smoke_dumping(
    pkg: &Path,
    script: &Path,
    max_ops: Option<u64>,
    timeout_ms: Option<f64>,
    dump: Option<&Path>,
) -> Value {
    let mut cmd = Command::new("node");
    cmd.arg(root().join("crates/popkern-live/tests/node/smoke.mjs"))
        .arg(pkg)
        .arg(script);
    cmd.arg(max_ops.map(|v| v.to_string()).unwrap_or_default());
    cmd.arg(timeout_ms.map(|v| v.to_string()).unwrap_or_default());
    if let Some(dir) = dump {
        cmd.arg(dir);
    }
    let output = cmd.output().expect("node");
    assert!(
        output.status.success(),
        "smoke.mjs failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("smoke JSON")
}

fn write(name: &str, body: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("popkern_live_wasm");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(name);
    std::fs::write(&path, body).expect("write");
    path
}

#[test]
fn wasm_bundle_answers_as_the_native_engine_does() {
    let Some(pkg) = bundle() else { return };
    let script = root().join("crates/popkern-rhai/examples/boss_plate.rhai");
    let dump = Path::new(env!("CARGO_TARGET_TMPDIR")).join("popkern_live_wasm/dump");
    let wasm = smoke_dumping(&pkg, &script, None, None, Some(&dump));
    assert_eq!(wasm["protocol_version"], "1");
    assert_eq!(wasm["clock_available"], true, "performance.now registered");
    assert_eq!(wasm["outcome"], "served", "{wasm}");
    assert_eq!(wasm["features"], 6);
    assert_eq!(wasm["report_failures"], 0);
    assert!(
        wasm["elapsed_ms"].as_f64().is_some(),
        "timing live under the clock"
    );
    let mut session = Session::new(HostOptions::default());
    let native = session.run(
        &std::fs::read_to_string(&script).expect("read"),
        "boss_plate.rhai",
    );
    assert!(native.outcome.is_served());
    assert_eq!(
        wasm["printed"],
        serde_json::to_value(&native.printed).expect("json")
    );
    let bodies = wasm["bodies"].as_array().expect("bodies");
    assert_eq!(bodies.len(), native.bodies.len());
    for (i, (w, n)) in bodies.iter().zip(&native.bodies).enumerate() {
        assert_eq!(w["typed"], true, "typed arrays cross the boundary");
        assert_eq!(w["vertex_count"], n.vertex_count);
        assert_eq!(w["triangle_count"], n.triangle_count);
        assert_eq!(w["face_count"], n.face_count);
        assert_eq!(w["open_edges"], n.open_edges);
        assert_eq!(
            w["face_labels"],
            serde_json::to_value(&n.face_labels).expect("json")
        );
        // Mesh topology crosses bit for bit; positions and normals are
        // compared value for value below (the two platforms' libm
        // differ in the last ulp of transcendental functions).
        assert_eq!(
            w["indices_fnv"],
            fnv1a_u32(&n.indices).to_string(),
            "indices"
        );
        assert_eq!(
            w["face_ids_fnv"],
            fnv1a_u32(&n.face_ids).to_string(),
            "face ids"
        );
        let positions = read_f32(&dump.join(format!("body{i}.positions.f32")));
        let normals = read_f32(&dump.join(format!("body{i}.normals.f32")));
        let dp = max_deviation(&positions, &n.positions);
        let dn = max_deviation(&normals, &n.normals);
        let identical_positions = same_hash(&w["positions_fnv"], fnv1a_f32(&n.positions));
        let identical_normals = same_hash(&w["normals_fnv"], fnv1a_f32(&n.normals));
        eprintln!(
            "[pk410] body {i}: positions max |wasm - native| = {dp:e} (bit-identical: \
             {identical_positions}), normals {dn:e} (bit-identical: {identical_normals})"
        );
        assert!(dp <= 1e-6, "positions deviate by {dp:e}");
        assert!(dn <= 1e-5, "normals deviate by {dn:e}");
    }
    // The document and its canonical script depend only on the tree:
    // bit-identical. The STEP text prints replayed geometry, so its
    // floats carry the same ulp noise as the mesh: entity for entity
    // the same file, and the wasm bytes pass the OCC oracle themselves.
    let step = session.export_step("smoke").expect("step");
    let pkn = session.export_pkn().expect("pkn");
    let rhai = session.to_rhai().expect("rhai");
    assert_eq!(wasm["exports"]["pkn_length"], pkn.len());
    assert_eq!(
        wasm["exports"]["pkn_fnv"],
        fnv1a_bytes(&pkn).to_string(),
        ".pkn bytes"
    );
    assert_eq!(
        wasm["exports"]["rhai_fnv"],
        fnv1a_bytes(rhai.as_bytes()).to_string(),
        "canonical Rhai"
    );
    let wasm_step = std::fs::read(dump.join("smoke.step")).expect("wasm STEP");
    let entities = |bytes: &[u8]| {
        String::from_utf8_lossy(bytes)
            .lines()
            .filter(|l| l.starts_with('#'))
            .map(|l| l.split('=').nth(1).unwrap_or("").trim().to_string())
            .map(|body| body.split('(').next().unwrap_or("").to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        entities(&wasm_step),
        entities(&step),
        "STEP entity sequence"
    );
    let identical_step = same_hash(&wasm["exports"]["step_fnv"], fnv1a_bytes(&step));
    eprintln!("[pk410] STEP bytes bit-identical across platforms: {identical_step}");
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = root().join("tools/occ-oracle/occ_check.py");
    if oracle.exists() && check.exists() {
        let output = Command::new(&oracle)
            .arg(&check)
            .arg(dump.join("smoke.step"))
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(
            report.contains("\"passed\": true"),
            "OCC oracle on wasm STEP: {report}"
        );
    } else {
        eprintln!("[pk410] OCC oracle absent; wasm STEP checked structurally only");
    }
    eprintln!(
        "[pk410] wasm boss plate: {} ms in wasm vs {:?} ms native",
        wasm["elapsed_ms"], native.elapsed_ms
    );
}

#[test]
fn wasm_bundle_stops_runaway_scripts_under_both_budgets() {
    let Some(pkg) = bundle() else { return };
    let spin = write(
        "spin.rhai",
        "let m = model();\nlet i = 0;\nloop { i += 1; }\n",
    );
    let ops = smoke(&pkg, &spin, Some(200_000), None);
    assert_eq!(ops["outcome"], "script_error", "{ops}");
    assert_eq!(ops["kind"], "BudgetExceeded");
    assert_eq!(ops["payload_budget"], "operations");
    assert_eq!(ops["wall_budget"], "none");
    let pattern = write(
        "pattern.rhai",
        "let m = model();\nlet i = 0;\nloop { m.profile(xy_at(i * 1.0)).rect(1.0, 1.0).sketch(\"s\" + i).extrude(blind(0.5)); i += 1; }\n",
    );
    let wall = smoke(&pkg, &pattern, None, Some(500.0));
    assert_eq!(wall["outcome"], "script_error", "{wall}");
    assert_eq!(wall["kind"], "BudgetExceeded");
    assert_eq!(wall["payload_budget"], "wall_clock");
    assert_eq!(
        wall["wall_budget"], "armed",
        "the registered clock enforces it"
    );
}

#[test]
fn wasm_bundle_classifies_declines_like_the_native_engine() {
    let Some(pkg) = bundle() else { return };
    let decline = write(
        "decline.rhai",
        "let m = model();\nlet cube = m.profile(xy()).rect(1.0, 1.0).sketch(\"cube\").extrude(blind(1.0));\ncube.fillet(edges(), 0.9);\nm\n",
    );
    let wasm = smoke(&pkg, &decline, None, None);
    assert_eq!(wasm["outcome"], "kernel_decline", "{wasm}");
    assert_eq!(wasm["kind"], "FeatureFailed");
    let syntax = write("syntax.rhai", "let m = model(\n");
    let wasm = smoke(&pkg, &syntax, None, None);
    assert_eq!(wasm["outcome"], "script_error");
    assert_eq!(wasm["kind"], "SyntaxError");
}

fn same_hash(reported: &Value, hash: u64) -> bool {
    reported.as_str() == Some(hash.to_string().as_str())
}

fn read_f32(path: &Path) -> Vec<f32> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let (words, _) = bytes.as_chunks::<4>();
    words.iter().map(|w| f32::from_le_bytes(*w)).collect()
}

fn max_deviation(a: &[f32], b: &[f32]) -> f64 {
    assert_eq!(a.len(), b.len(), "buffer lengths");
    a.iter()
        .zip(b)
        .map(|(x, y)| (f64::from(*x) - f64::from(*y)).abs())
        .fold(0.0, f64::max)
}

fn fnv1a_bytes(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// PK-412: every gallery script runs in the wasm bundle with the
/// feature count the manifest promises.
#[test]
fn wasm_bundle_serves_every_gallery_script() {
    let Some(pkg) = bundle() else { return };
    let manifest: Value = serde_json::from_str(
        &std::fs::read_to_string(root().join("web/live/examples.json")).expect("manifest"),
    )
    .expect("manifest json");
    let entries = manifest["examples"].as_array().expect("examples");
    assert!(entries.len() >= 8);
    for entry in entries {
        // Project examples build through run_project (project_parity.rs).
        let Some(file) = entry["file"].as_str() else {
            continue;
        };
        let script = root().join("crates/popkern-rhai/examples").join(file);
        let out = smoke(&pkg, &script, None, None);
        assert_eq!(out["outcome"], "served", "[{file}] {out}");
        assert_eq!(out["features"], entry["features"], "[{file}] feature count");
        assert_eq!(out["report_failures"], 0, "[{file}]");
        for body in out["bodies"].as_array().expect("bodies") {
            assert_eq!(body["open_edges"], 0, "[{file}] watertight in wasm");
        }
    }
}
