//! PK-583: a step whose filleted reflex edge ends in two saddle
//! vertices (six edges, equal radius 0.06) serves both saddles through
//! the saddle setback fan.
//!
//! It used to achieve G0, with five seams creased at g1 0.42-0.49 at
//! both requests. Four defects stacked:
//! - The sphere pre-screen claimed the equal-radius saddles for a
//!   rolling ball, so the G1 request never reached the setback rung.
//! - The x=0 saddle's ring carve registered the support plane as a
//!   strip: the chord `kef` kept the strip's handle. Its bridge then
//!   found no support.
//! - Plane pcurves projected from a split NURBS rail spanned the whole
//!   rail (G2 request: 0.237).
//! - With the fans serving, the C1 pair-corner flap arrived split in
//!   two by a mitre chord, and the settlement skipped it.
//!
//! A saddle cannot meet G2, so the G2 request downgrades honestly to a
//! certified G1.

use popkern_rhai::{ScriptHost, report_value};

fn step_rounds(continuity: &str) -> serde_json::Value {
    let script = format!(
        r#"let m = model();
let step = m.profile(yz()).polygon([[0.0, 3.0], [1.0, 3.0], [1.0, 3.6], [0.6, 3.6], [0.6, 4.0], [0.0, 4.0]]).sketch("step").extrude(blind(1.0));
step.fillet(edges().straight().direction("x").at_height(3.6)
    .or(edges().straight().direction("y").at_height(3.6).chord_length_between(0.39, 0.41))
    .or(edges().straight().direction("z").chord_length_between(0.39, 0.41)),
    0.06, #{{ continuity: "{continuity}", label: "rounds" }});
m
"#
    );
    let mut host = ScriptHost::new();
    let outcome = host.run(&script, "step.rhai").expect("the step serves");
    let value = report_value(&outcome.model.expect("model"));
    value["features"]
        .as_array()
        .expect("rows")
        .iter()
        .find(|r| r["label"] == "rounds")
        .cloned()
        .expect("rounds row")
}

fn assert_certified_g1(rounds: &serde_json::Value, request: &str) {
    let blend = &rounds["blend"];
    let c = &blend["continuity"];
    // Two setback fans of 16 pieces each.
    assert_eq!(blend["corner_patches"], 32, "{request}: {blend}");
    assert_eq!(c["achieved"], "G1", "{request}: {c}");
    assert_eq!(c["certified_achieved"], "G1", "{request}");
    assert_eq!(c["uncertified_seams"], 0, "{request}");
    assert_eq!(c["unmeasured_seams"], 0, "{request}");
    // Measured 8.0e-8 (G1 request) and 3.3e-8 (G2 request); the old
    // creases read 0.42-0.49.
    let g1 = c["worst_g1_normal"].as_f64().expect("g1");
    assert!(g1 < 1e-6, "{request}: worst g1 {g1:.3e}");
    // The commercial G0 gate. Measured 1.9e-10 (G1 request) and 3.4e-6
    // (G2 request). The whole-rail plane pcurves read 0.225.
    let g0 = c["worst_g0_gap"].as_f64().expect("g0");
    assert!(g0 < 1e-5, "{request}: worst g0 {g0:.3e}");
}

#[test]
fn a_g1_request_serves_both_step_saddles_certified_g1() {
    let rounds = step_rounds("g1");
    assert_certified_g1(&rounds, "G1");
    assert_eq!(rounds["blend"]["downgraded"], false);
}

#[test]
fn a_g2_request_on_the_step_downgrades_honestly_to_certified_g1() {
    let rounds = step_rounds("g2");
    assert_certified_g1(&rounds, "G2");
    assert_eq!(rounds["blend"]["downgraded"], true);
    assert_eq!(rounds["blend"]["continuity"]["requested"], "G2");
}
