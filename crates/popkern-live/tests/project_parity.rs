//! PK-660 parity: a project built by the wasm bundle (the node smoke
//! drives `LiveSession::run_project` with the folder's files, the way
//! the page hands a project over) answers as the native session —
//! outcome, feature count, every migration's row, and each body's mesh
//! volume (to 1e-9 relative: the platforms' libm differ in the last
//! ulp of transcendental functions). Skips like the wasm smoke when the
//! nodejs bundle is missing or stale.

#[path = "support/node_bundle.rs"]
mod node_bundle;

use std::path::{Path, PathBuf};
use std::process::Command;

use node_bundle::{bundle, root};
use popkern_live::{HostOptions, Session};
use serde_json::Value;

fn fixture(name: &str) -> PathBuf {
    root()
        .join("crates/popkern-project/tests/fixtures")
        .join(name)
}

fn files(dir: &Path, from: &Path, out: &mut Vec<(String, Vec<u8>)>) {
    for entry in std::fs::read_dir(dir).expect("fixture") {
        let path = entry.expect("entry").path();
        if path.is_dir() {
            files(&path, from, out);
        } else {
            let key = path
                .strip_prefix(from)
                .expect("inside")
                .to_string_lossy()
                .replace('\\', "/");
            out.push((key, std::fs::read(&path).expect("read")));
        }
    }
}

fn mesh_volume(positions: &[f32], indices: &[u32]) -> f64 {
    let p = |i: u32| {
        let i = 3 * i as usize;
        [
            positions[i] as f64,
            positions[i + 1] as f64,
            positions[i + 2] as f64,
        ]
    };
    indices
        .chunks_exact(3)
        .map(|t| {
            let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        })
        .sum::<f64>()
        / 6.0
}

#[test]
fn projects_build_the_same_in_wasm_and_natively() {
    let Some(pkg) = bundle() else {
        return;
    };
    for name in ["boss_plate", "labelled_rim", "step_base"] {
        let dir = fixture(name);
        let output = Command::new("node")
            .arg(root().join("crates/popkern-live/tests/node/project_smoke.mjs"))
            .arg(&pkg)
            .arg(&dir)
            .output()
            .expect("node");
        assert!(
            output.status.success(),
            "project_smoke.mjs failed:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let wasm: Value = serde_json::from_slice(&output.stdout).expect("smoke JSON");

        let mut all = Vec::new();
        files(&dir, &dir, &mut all);
        let mut session = Session::new(HostOptions::default());
        let native = session.run_project(all);
        let native_json = native.to_json();
        assert_eq!(wasm["outcome"], native_json["outcome"], "{name}: outcome");
        assert_eq!(
            wasm["features"], native_json["features"],
            "{name}: features"
        );
        let rows = |v: &Value| -> Vec<(Value, Value, Value, Value)> {
            v["project"]["versions"]
                .as_array()
                .expect("versions")
                .iter()
                .map(|r| {
                    (
                        r["name"].clone(),
                        r["status"].clone(),
                        r["features"].clone(),
                        r["edited"].clone(),
                    )
                })
                .collect()
        };
        assert_eq!(rows(&wasm), rows(&native_json), "{name}: migration rows");
        let volumes: Vec<f64> = native
            .bodies
            .iter()
            .map(|b| mesh_volume(&b.positions, &b.indices))
            .collect();
        let wasm_volumes: Vec<f64> = wasm["volumes"]
            .as_array()
            .expect("volumes")
            .iter()
            .map(|v| v.as_f64().expect("number"))
            .collect();
        assert_eq!(wasm_volumes.len(), volumes.len(), "{name}: body count");
        for (w, n) in wasm_volumes.iter().zip(&volumes) {
            assert!(
                (w - n).abs() <= 1e-9 * n.abs().max(1e-12),
                "{name}: volume {w} (wasm) vs {n} (native)"
            );
        }
        assert!(volumes.iter().all(|v| *v > 0.0), "{name}: solid bodies");
    }
}

#[test]
fn a_session_reuses_checkpoints_across_project_runs() {
    let dir = fixture("boss_plate");
    let mut all = Vec::new();
    files(&dir, &dir, &mut all);
    let mut session = Session::new(HostOptions::default());
    let first = session.run_project(all.clone());
    assert!(first.outcome.is_served(), "{:?}", first.outcome);
    for (path, bytes) in &mut all {
        if path.ends_with("003_rounds.rhai") {
            *bytes = String::from_utf8(bytes.clone())
                .expect("utf8")
                .replace("0.08", "0.06")
                .into_bytes();
        }
    }
    let second = session.run_project(all).to_json();
    let statuses: Vec<&str> = second["project"]["versions"]
        .as_array()
        .expect("versions")
        .iter()
        .map(|r| r["status"].as_str().expect("status"))
        .collect();
    assert_eq!(statuses, ["cached", "cached", "built"]);
}
