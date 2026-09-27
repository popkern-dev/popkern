//! PK-411: the report JSON every host shows carries per-row derived
//! facts — wall time in milliseconds and the blend's downgrade
//! verdict — so a panel never has to re-derive them (and never shows
//! a bare success mark for a blend whose measured continuity fell
//! below the request). The fixture pairs a blend that meets its G2
//! request (a cube's top edges) with one that serves a valid body but
//! measures below it: a step whose reflex edge ends in two saddle
//! vertices, rounded at G2. A saddle can never meet G2 (at each bridge
//! foot the plane asks for curvature 0 and the strip for 1/r), so the
//! downgrade row does not depend on construction quality. It reads a
//! certified G1 since PK-583 (it had been a G0 crease; see
//! `step_saddle_setback.rs`). The triangular prism this used before
//! meets G2 since PK-575.
//!
//! PK-530: the fixture was a blend on a blend, whose second G2 fillet
//! served corner residue off its faces. The territory gate now
//! declines that body (PK-534), so the downgrade row comes from a
//! served, valid blend instead.

use popkern_rhai::{ScriptHost, report_value};

const MET_AND_DOWNGRADED: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges().at_height(1.0), 0.15, #{ continuity: "g2", label: "top" });
let step = m.profile(yz()).polygon([[0.0, 3.0], [1.0, 3.0], [1.0, 3.6], [0.6, 3.6], [0.6, 4.0], [0.0, 4.0]]).sketch("step").extrude(blind(1.0));
step.fillet(edges().straight().direction("x").at_height(3.6)
    .or(edges().straight().direction("y").at_height(3.6).chord_length_between(0.39, 0.41))
    .or(edges().straight().direction("z").chord_length_between(0.39, 0.41)),
    0.06, #{ continuity: "g2", label: "rounds" });
m
"#;

#[test]
fn report_rows_carry_elapsed_ms_and_the_blend_downgrade_verdict() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(MET_AND_DOWNGRADED, "met_and_downgraded.rhai")
        .expect("served");
    let model = outcome.model.expect("model");
    let value = report_value(&model);
    assert_eq!(value["failures"], 0);
    assert_eq!(value["downgraded_blends"], 1, "{value}");
    let rows = value["features"].as_array().expect("rows");
    assert_eq!(rows.len(), 6, "two sketches, two extrudes, two fillets");
    for row in rows {
        assert!(
            row["elapsed_ms"].as_f64().is_some(),
            "every executed row carries elapsed_ms: {row}"
        );
    }
    let top = rows.iter().find(|r| r["label"] == "top").expect("top row");
    assert_eq!(top["blend"]["downgraded"], false);
    assert_eq!(top["blend"]["continuity"]["requested"], "G2");
    assert_eq!(top["blend"]["continuity"]["achieved"], "G2");
    let rounds = rows
        .iter()
        .find(|r| r["label"] == "rounds")
        .expect("rounds row");
    assert_eq!(rounds["blend"]["downgraded"], true, "{rounds}");
    assert_eq!(rounds["blend"]["continuity"]["requested"], "G2");
    assert_eq!(rounds["blend"]["continuity"]["achieved"], "G1");
    // The recorder agrees with the projection.
    assert_eq!(model.report().downgraded_blends().count(), 1);
}

/// Blend on blend at G2: a G2 fillet on the vertical edges of a cube
/// whose top edges already carry a G2 fillet. It once served residue
/// off its faces (8 edges, worst 3.6e-2). The vertical strips' ends now
/// settle against the top strips (PK-534), and the fitted G2 verticals
/// reach past their end stations to do so (PK-531). The body serves
/// with every edge on its faces.
#[test]
fn blend_on_blend_g2_serves_consistent() {
    const BLEND_ON_BLEND: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube = cube.fillet(edges().at_height(1.0), 0.15, #{ continuity: "g2", label: "top" });
cube.fillet(edges().straight().direction("z"), 0.1, #{ continuity: "g2", label: "verticals" });
m
"#;
    let mut host = ScriptHost::new();
    let outcome = host
        .run(BLEND_ON_BLEND, "blend_on_blend.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let errors = model.to_part().expect("part").brep().validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
}
