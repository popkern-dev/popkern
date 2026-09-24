//! Multi-body documents export and render in a deterministic body
//! order (PK-525).
//!
//! The history's body manager keys bodies in a `HashMap`, and
//! `bodies()` once returned them in hash iteration order: the
//! whole-model STEP and thumbnail concatenated bodies in a random
//! order per map instance (measured on the countersunk-joint gallery
//! script: STEP data and thumbnail pixels changed run to run, the
//! screw head z-fighting its countersink in draw order), and
//! first-match body scans (auto-merge host detection, sweep path
//! lookup, the current-body fallback) inherited the same randomness.
//! Bodies now come back in ascending id — creation — order.

use popkern_rhai::ScriptHost;

/// Eight independent bodies (each base-plane extrude is its own
/// body), enough that a hash-ordered concatenation is a different
/// permutation on practically every run.
const EIGHT_BLOCKS: &str = r#"
let m = model();
for k in 0..8 {
    m.profile(xy()).rect(1.0, 1.0).at(3.0 * k, 0.0).sketch("block_" + k).extrude(blind(1.0 + 0.25 * k));
}
m
"#;

fn whole_model_step(script: &str) -> Vec<u8> {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(script, "blocks.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    assert_eq!(model.alive_bodies().len(), 8, "eight live bodies");
    model
        .to_part_all()
        .expect("whole-model part")
        .step_bytes("blocks")
        .expect("STEP bytes")
}

#[test]
fn whole_model_step_is_identical_across_runs() {
    let first = whole_model_step(EIGHT_BLOCKS);
    for run in 1..5 {
        assert!(
            whole_model_step(EIGHT_BLOCKS) == first,
            "run {run}: the whole-model STEP must be byte-identical to run 0"
        );
    }
}

#[test]
fn bodies_come_back_in_creation_order() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(EIGHT_BLOCKS, "blocks.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let ids: Vec<u32> = model.alive_bodies().iter().map(|body| body.id.0).collect();
    let mut sorted = ids.clone();
    sorted.sort_unstable();
    assert_eq!(ids, sorted, "alive bodies in ascending (creation) order");
}
