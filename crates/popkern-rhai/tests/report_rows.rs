//! PK-411: the report JSON every host shows carries per-row derived
//! facts — wall time in milliseconds and the blend's downgrade
//! verdict — so a panel never has to re-derive them (and never shows
//! a bare success mark for a blend whose measured continuity fell
//! below the request). The fixture pairs a blend that meets its G2
//! request (a cube's top edges) with one that serves a valid body but
//! measures below it: every edge of a triangular prism at G2, whose
//! corners land at G1.
//!
//! PK-530: the fixture was a blend on a blend, whose second G2 fillet
//! served corner residue off its faces. The territory gate now
//! declines that body (PK-534), so the downgrade row comes from a
//! served, valid blend instead.

use popkern_rhai::{ScriptHost, report_value};

const MET_AND_DOWNGRADED: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges().at_height(1.0), 0.15, #{ continuity: "g2", label: "top" });
let prism = m.profile(xy()).polygon([[3.0, 0.0], [4.0, 0.0], [3.5, 0.8]]).sketch("prism").extrude(blind(1.0));
prism.fillet(edges(), 0.08, #{ continuity: "g2", label: "rounds" });
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
    assert_ne!(rounds["blend"]["continuity"]["achieved"], "G2");
    // The recorder agrees with the projection.
    assert_eq!(model.report().downgraded_blends().count(), 1);
}

/// PK-531 pin: a G2 fillet on the vertical edges of a cube whose top
/// edges already carry a G2 fillet once served residue off its faces
/// (8 edges, worst 3.6e-2). G1 verticals on a G1 or G2 top now serve
/// (PK-534). G2 verticals are fitted strips whose surface stops at the
/// edge's end station, short of the crossing with the top strips' seam,
/// so their ends cannot settle until fitted strips extend (PK-531). The
/// territory gate declines them typed, as `InconsistentGeometry` in the
/// failure payload.
#[test]
fn blend_on_blend_g2_residue_declines_typed() {
    const BLEND_ON_BLEND: &str = r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube = cube.fillet(edges().at_height(1.0), 0.15, #{ continuity: "g2", label: "top" });
cube.fillet(edges().straight().direction("z"), 0.1, #{ continuity: "g2", label: "verticals" });
m
"#;
    let mut host = ScriptHost::new();
    let Err(error) = host.run(BLEND_ON_BLEND, "blend_on_blend.rhai") else {
        panic!("PK-531 served the G2 blend-on-blend corners: restore it as the report fixture");
    };
    assert_eq!(error.kind, "FeatureFailed");
    assert!(
        error.payload.to_string().contains("InconsistentGeometry"),
        "typed territory decline: {}",
        error.payload
    );
}
