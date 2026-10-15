//! The clipping view's geometry (`web/live/section-geometry.js`;
//! docs/architecture/clipping-view.md §9): the pure checks in
//! `test/section.test.mjs`, and the cut outline on the kernel's own
//! display mesh — the boss plate cut at several heights and tilts
//! outlines a CLOSED polyline (every endpoint shared by an even number
//! of segments). Skips when node is not installed.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use popkern_live::{HostOptions, Session};
use serde_json::{Value, json};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn node(script: &str, input: Option<&Value>) -> Option<Vec<u8>> {
    let mut child = match Command::new("node")
        .arg(root().join("web/live/test").join(script))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
    {
        Ok(child) => child,
        Err(_) => {
            eprintln!("[page_section] SKIP: node is not installed");
            return None;
        }
    };
    let bytes = input
        .map(|v| serde_json::to_vec(v).expect("json"))
        .unwrap_or_default();
    std::io::Write::write_all(&mut child.stdin.take().expect("stdin"), &bytes).expect("write");
    let output = child.wait_with_output().expect("node");
    assert!(
        output.status.success(),
        "{script} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    Some(output.stdout)
}

#[test]
fn section_geometry_checks_pass() {
    let Some(out) = node("section.test.mjs", None) else {
        return;
    };
    let row: Value = serde_json::from_slice(&out).expect("json");
    assert!(row["checks"].as_u64().unwrap_or(0) >= 14, "{row}");
}

#[test]
fn a_kernel_display_mesh_cuts_to_a_closed_outline() {
    let script = root().join("crates/popkern-rhai/examples/boss_plate.rhai");
    let mut session = Session::new(HostOptions::default());
    session.run(
        &std::fs::read_to_string(&script).expect("script"),
        "boss_plate.rhai",
    );
    let steps = session.timeline().expect("timeline")[0].steps.len();
    let rows = session.step_bodies(None, steps, &[]).expect("final rows");
    let mesh = rows[0].mesh.as_ref().expect("mesh");
    let planes = json!([
        { "normal": [0.0, 0.0, 1.0], "offset": 0.2 },
        { "normal": [0.0, 0.0, 1.0], "offset": 0.55 },
        { "normal": [1.0, 0.0, 0.0], "offset": 0.0 },
        { "normal": [0.6, 0.0, 0.8], "offset": 0.3 },
        { "normal": [0.48, 0.6, 0.64], "offset": 0.1 },
    ]);
    let input = json!({
        "meshes": [{ "positions": mesh.positions, "indices": mesh.indices }],
        "planes": planes,
    });
    let Some(out) = node("section-mesh.test.mjs", Some(&input)) else {
        return;
    };
    let rows: Vec<Value> = serde_json::from_slice(&out).expect("rows");
    assert_eq!(rows.len(), 5);
    for (k, row) in rows.iter().enumerate() {
        assert!(
            row["segments"].as_u64().unwrap_or(0) > 4,
            "plane {k}: {row}"
        );
        assert_eq!(row["open"], 0, "plane {k}: outline not closed: {row}");
    }
}
