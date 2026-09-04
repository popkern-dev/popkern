//! PK-436 / PK-438 regression: a keyed bore that protrudes below a
//! union-derived body serves watertight. This is the sketch-bridge
//! reproducer for BOTH stacked defects — the analytic `BRep::cylinder`
//! path does not reproduce them, so the acceptance lives here.
//!
//! A carrier disc (union of a disc, three pins and an output boss) is
//! bored by `hole ∪ key`, a tool that overshoots the disc bottom by
//! 0.06. Defect 1 (PK-436) opened the boss-top keyway loop; defect 2
//! (PK-438) glued the below-body protrusion into the boss-region wall
//! fragment and stranded the flush rim. With both fixed the body is
//! watertight — so `assert_watertight` inside the script serving at all
//! is the gate (it raises `NotWatertight` otherwise), and we also check
//! the recorded report carries no failures.

use popkern_rhai::ScriptHost;

const SCRIPT: &str = r#"
let m = model();
let pi = 3.141592653589793;
let planets = 3;
let carrier_r = 0.5;
let carrier = m.profile(xy_at(0.44)).circle(0.62).sketch("carrier_disc").extrude(blind(0.12));
for k in 0..planets {
    let phi = k * 2.0 * pi / planets;
    let pin = m.tool(|m| m.profile(xy_at(0.02)).circle(0.07).at(carrier_r * cos(phi), carrier_r * sin(phi)).sketch("pin_" + k).extrude(blind(0.48)));
    carrier = carrier.union(pin);
}
let boss = m.tool(|m| m.profile(xy_at(0.52)).circle(0.16).sketch("output_boss").extrude(blind(0.28)));
carrier = carrier.union(boss);
let out_bore = m.tool(|m| {
    let hole = m.profile(xy_at(0.40)).circle(0.1).sketch("output_bore").extrude(blind(0.45));
    let key = m.profile(xy_at(0.38)).rect(0.06, 0.05).at(0.0, 0.115).sketch("output_key").extrude(blind(0.51));
    hole.union(key)
});
carrier = carrier.subtract(out_bore);
m.assert_watertight();
m
"#;

#[test]
fn keyed_bore_over_union_serves_watertight() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(SCRIPT, "keyed_bore_over_union.rhai")
        .unwrap_or_else(|e| {
            panic!(
                "keyed bore over union must serve watertight: {e}\n{}",
                e.payload
            )
        });
    let scripted = outcome.model.expect("the script returns its model");
    assert_eq!(
        scripted.report().failures().count(),
        0,
        "no feature may fail"
    );
}
