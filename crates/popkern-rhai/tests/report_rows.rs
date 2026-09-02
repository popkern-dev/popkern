//! PK-411: the report JSON every host shows carries per-row derived
//! facts — wall time in milliseconds and the blend's downgrade
//! verdict — so a panel never has to re-derive them (and never shows
//! a bare success mark for a blend whose measured continuity fell
//! below the request). The fixture is a blend on a blend: the second
//! G2 fillet lands at G0 today, which is exactly the row a user must
//! see flagged.

use popkern_rhai::{ScriptHost, report_value};

const BLEND_ON_BLEND: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube = cube.fillet(edges().at_height(1.0), 0.15, #{ continuity: "g2", label: "top" });
cube.fillet(edges().straight().direction("z"), 0.1, #{ continuity: "g2", label: "verticals" });
m
"#;

#[test]
fn report_rows_carry_elapsed_ms_and_the_blend_downgrade_verdict() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(BLEND_ON_BLEND, "blend_on_blend.rhai")
        .expect("served");
    let model = outcome.model.expect("model");
    let value = report_value(&model);
    assert_eq!(value["failures"], 0);
    assert_eq!(value["downgraded_blends"], 1, "{value}");
    let rows = value["features"].as_array().expect("rows");
    assert_eq!(rows.len(), 4);
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
    let verticals = rows
        .iter()
        .find(|r| r["label"] == "verticals")
        .expect("verticals row");
    assert_eq!(verticals["blend"]["downgraded"], true, "{verticals}");
    assert_eq!(verticals["blend"]["continuity"]["requested"], "G2");
    assert_ne!(verticals["blend"]["continuity"]["achieved"], "G2");
    // The recorder agrees with the projection.
    assert_eq!(model.report().downgraded_blends().count(), 1);
}
