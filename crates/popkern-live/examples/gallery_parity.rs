//! Native vs wasm parity over the gallery (live-projects-and-history.md
//! §11 E-parity, PK-665): every gallery script runs in the nodejs wasm
//! bundle (tests/node/smoke.mjs, buffers dumped) and natively; per body
//! the topology (counts, triangle indices, face ids) must be identical,
//! the mesh volume must agree to 1e-9 relative, and the share of
//! bit-identical vertex buffers is reported (the two platforms' libm
//! may differ in the last ulp of transcendental functions).
//!
//! `cargo run --profile live -p popkern-live --example gallery_parity`
//! (after `tools/build-live.sh --node`). Exit status 1 on a topology or
//! volume disagreement.

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_live::mesh_buffers::{fnv1a_f32, fnv1a_u32};
use popkern_live::{HostOptions, Session};
use serde_json::Value;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn volume(positions: &[f32], indices: &[u32]) -> f64 {
    let p = |i: u32| {
        let i = 3 * i as usize;
        [
            f64::from(positions[i]),
            f64::from(positions[i + 1]),
            f64::from(positions[i + 2]),
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

fn read_words(path: &Path) -> Vec<[u8; 4]> {
    std::fs::read(path)
        .unwrap_or_default()
        .chunks_exact(4)
        .map(|b| [b[0], b[1], b[2], b[3]])
        .collect()
}

fn main() -> Result<std::process::ExitCode, Box<dyn std::error::Error>> {
    let pkg = root().join("target/live/node");
    if !pkg.join("popkern_live.js").exists() {
        eprintln!("no nodejs bundle: run tools/build-live.sh --node");
        return Ok(std::process::ExitCode::from(2));
    }
    let manifest: Value = serde_json::from_str(&std::fs::read_to_string(
        root().join("web/live/examples.json"),
    )?)?;
    let dump_root = root().join("target/live/parity");
    let (mut bodies, mut identical, mut failures) = (0usize, 0usize, Vec::new());
    let mut worst_volume = 0.0f64;
    println!(
        "{:<28} {:>6} {:>12} {:>14}",
        "script", "bodies", "bit-equal", "max |dV|/V"
    );
    for entry in manifest["examples"].as_array().cloned().unwrap_or_default() {
        let Some(file) = entry["file"].as_str() else {
            continue;
        };
        let script = root().join("crates/popkern-rhai/examples").join(file);
        let dump = dump_root.join(file);
        let _ = std::fs::remove_dir_all(&dump);
        let output = Command::new("node")
            .arg(root().join("crates/popkern-live/tests/node/smoke.mjs"))
            .arg(&pkg)
            .arg(&script)
            .arg("")
            .arg("")
            .arg(&dump)
            .output()?;
        let wasm: Value = serde_json::from_slice(&output.stdout)?;
        let mut session = Session::new(HostOptions::default());
        let source = std::fs::read_to_string(&script)?;
        let native = session.run(&source, file);
        let native_json = native.to_json();
        if wasm["outcome"] != native_json["outcome"] || wasm["features"] != native_json["features"]
        {
            failures.push(format!("{file}: outcome/features differ"));
            continue;
        }
        let wasm_bodies = wasm["bodies"].as_array().cloned().unwrap_or_default();
        let (mut equal, mut worst) = (0, 0.0f64);
        for (i, (w, n)) in wasm_bodies.iter().zip(&native.bodies).enumerate() {
            let topology = w["vertex_count"] == n.vertex_count
                && w["triangle_count"] == n.triangle_count
                && w["face_count"] == n.face_count
                && w["indices_fnv"].as_str() == Some(fnv1a_u32(&n.indices).to_string().as_str())
                && w["face_ids_fnv"].as_str() == Some(fnv1a_u32(&n.face_ids).to_string().as_str());
            if !topology {
                failures.push(format!("{file} body {i}: topology differs"));
            }
            if w["positions_fnv"].as_str() == Some(fnv1a_f32(&n.positions).to_string().as_str()) {
                equal += 1;
            }
            // Each platform's volume from its own buffers: with a
            // differing triangulation, crossing them measures nothing.
            let dumped = |kind: &str| read_words(&dump.join(format!("body{i}.{kind}")));
            let wasm_positions: Vec<f32> = dumped("positions.f32")
                .into_iter()
                .map(f32::from_le_bytes)
                .collect();
            let wasm_indices: Vec<u32> = dumped("indices.u32")
                .into_iter()
                .map(u32::from_le_bytes)
                .collect();
            let (vw, vn) = (
                volume(&wasm_positions, &wasm_indices),
                volume(&n.positions, &n.indices),
            );
            let rel = (vw - vn).abs() / vn.abs().max(1e-12);
            worst = worst.max(rel);
            if rel > 1e-9 {
                failures.push(format!("{file} body {i}: volume {vw} vs {vn}"));
            }
        }
        if wasm_bodies.len() != native.bodies.len() {
            failures.push(format!("{file}: body count differs"));
        }
        bodies += native.bodies.len();
        identical += equal;
        worst_volume = worst_volume.max(worst);
        println!(
            "{file:<28} {:>6} {:>12} {worst:>14.2e}",
            native.bodies.len(),
            format!("{equal}/{}", native.bodies.len())
        );
    }
    println!(
        "\n{bodies} bodies: {identical} bit-identical vertex buffers ({:.0}%), worst volume disagreement {worst_volume:.2e}",
        100.0 * identical as f64 / bodies.max(1) as f64
    );
    if failures.is_empty() {
        println!("topology and outcomes identical; volumes within 1e-9");
        Ok(std::process::ExitCode::SUCCESS)
    } else {
        for f in &failures {
            println!("FAIL {f}");
        }
        Ok(std::process::ExitCode::from(1))
    }
}
