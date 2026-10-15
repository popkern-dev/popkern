//! The viewer's face-boundary edges (`web/live/edges.js`) on the
//! kernel's own step meshes: every drawn segment lies on a model edge.
//!
//! A sketch step's body is a lamina: front and back faces share every
//! triangle edge, and the earlier "two different faces meet" rule drew
//! the region's whole triangulation (the boss plate's history step 1).
//! Edges now bound a face when that face uses the mesh edge an odd
//! number of times. The same lamina is drawn front-sided (its two
//! coincident faces depth-fight double-sided), so the test also pins
//! `enclosesNoVolume`.

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

/// The boss plate's meshes after `steps` (the sketch, the extrude).
fn boss_plate_steps(steps: &[usize]) -> Vec<Value> {
    let script = root().join("crates/popkern-rhai/examples/boss_plate.rhai");
    let mut session = Session::new(HostOptions::default());
    session.run(
        &std::fs::read_to_string(&script).expect("script"),
        "boss_plate.rhai",
    );
    steps
        .iter()
        .map(|&index| {
            let rows = session.step_bodies(None, index, &[]).expect("step rows");
            assert_eq!(rows.len(), 1, "step {index}: one body");
            let mesh = rows[0].mesh.as_ref().expect("mesh");
            json!({
                "positions": mesh.positions,
                "indices": mesh.indices,
                "face_ids": mesh.face_ids,
            })
        })
        .collect()
}

fn boundary_rows(meshes: &[Value]) -> Option<Vec<Value>> {
    let mut child = match Command::new("node")
        .arg(root().join("web/live/test/edges.test.mjs"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
    {
        Ok(child) => child,
        Err(_) => {
            eprintln!("[page_edges] SKIP: node is not installed");
            return None;
        }
    };
    let input = serde_json::to_vec(meshes).expect("json");
    std::io::Write::write_all(&mut child.stdin.take().expect("stdin"), &input).expect("write");
    let output = child.wait_with_output().expect("node");
    assert!(output.status.success(), "edges.test.mjs failed");
    Some(serde_json::from_slice(&output.stdout).expect("rows"))
}

/// Points `[x, y, z]` of a row's segments.
fn points(row: &Value) -> Vec<[f64; 3]> {
    let flat: Vec<f64> = row["segments"]
        .as_array()
        .expect("segments")
        .iter()
        .map(|v| v.as_f64().expect("number"))
        .collect();
    flat.as_chunks::<3>().0.to_vec()
}

#[test]
fn sketch_laminae_and_solids_draw_only_their_model_edges() {
    let Some(rows) = boundary_rows(&boss_plate_steps(&[1, 2])) else {
        return;
    };
    // The plate: rect(2.4, 1.6) centred on the origin, extruded 0.4.
    let on = |v: f64, target: f64| (v.abs() - target).abs() < 1e-5;
    let tol = 1e-5;
    // Step 1, the plate sketch: one front/back face pair, drawn only
    // around the rectangle.
    let sketch = points(&rows[0]);
    assert!(!sketch.is_empty(), "the sketch outline is drawn");
    assert_eq!(rows[0]["edges"], 1, "one front/back face pair");
    assert_eq!(rows[0]["sheet"], true, "a lamina encloses no volume");
    assert_eq!(rows[1]["sheet"], false, "the extruded plate does");
    for p in &sketch {
        assert!(p[2].abs() < tol, "{p:?} off the sketch plane");
        assert!(
            on(p[0], 1.2) || on(p[1], 0.8),
            "{p:?} is inside the rectangle: a triangulation edge was drawn"
        );
    }
    // Step 2, the extruded plate: a box's twelve edges and nothing else.
    let solid = points(&rows[1]);
    assert_eq!(rows[1]["edges"], 12, "a box has twelve edges");
    for p in &solid {
        let on_count = [
            on(p[0], 1.2),
            on(p[1], 0.8),
            p[2].abs() < tol || (p[2] - 0.4).abs() < tol,
        ]
        .iter()
        .filter(|b| **b)
        .count();
        assert!(on_count >= 2, "{p:?} is not on a box edge");
    }
}
