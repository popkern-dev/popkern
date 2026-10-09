//! The six-axis arm's sliders in the wasm engine (design
//! `six-axis-arm.md` §8: a joint re-solve within 30 ms, an `upper_arm`
//! edit — one tube replay and a re-solve — within 6 s; §9.8): the
//! nodejs bundle builds the arm, then a drag of `j2` answers with
//! placement-only envelopes (the engine keeps its body buffers) with a
//! step within the budget, and the tube edit rebuilds within
//! its own. Skips like the wasm smoke when the bundle is missing or
//! stale.

#[path = "support/node_bundle.rs"]
mod node_bundle;

use std::process::Command;

use node_bundle::{bundle, root};
use serde_json::Value;

#[test]
fn the_arm_s_joint_slider_steps_within_the_wasm_budget() {
    let Some(pkg) = bundle() else {
        return;
    };
    let output = Command::new("node")
        .arg(root().join("crates/popkern-live/tests/node/joint_steps.mjs"))
        .arg(&pkg)
        .arg(root().join("crates/popkern-rhai/examples/six_axis_arm.rhai"))
        .args(["j2", "30", "21", "upper_arm", "5"])
        .output()
        .expect("node");
    let out: Value = serde_json::from_slice(&output.stdout).unwrap_or_else(|e| {
        panic!(
            "joint_steps output: {e}\n{}",
            String::from_utf8_lossy(&output.stderr)
        )
    });
    eprintln!("[pk725] wasm arm: {out}");
    assert_eq!(out["outcome"], "served", "{out}");
    assert_eq!(out["served"], true, "{out}");
    assert_eq!(out["bodies_unchanged"], true, "{out}");
    assert_eq!(out["bodies_before"], out["bodies_after"], "{out}");
    // Every step does the same work: the fastest is its cost (a median
    // also carries the machine's load).
    let fastest = out["fastest_ms"].as_f64().unwrap_or(f64::MAX);
    assert!(
        fastest < 30.0,
        "j2 step {fastest} ms over the 30 ms wasm budget"
    );
    let tube = &out["rebuilt"];
    assert_eq!(tube["outcome"], "served", "{out}");
    assert_eq!(
        tube["bodies_unchanged"], false,
        "a tube edit re-meshes: {out}"
    );
    let ms = tube["ms"].as_f64().unwrap_or(f64::MAX);
    assert!(
        ms < 6000.0,
        "upper_arm edit {ms} ms over the 6 s wasm budget"
    );
}
