//! PK-603: the `corner_cap` / `corner_cap_strict` fillet options. The
//! preference parses, round-trips through "show as code", is refused
//! where it cannot apply, and the report names what served each corner.

use popkern_rhai::{ScriptHost, report_value, to_rhai};

/// The unit block's (+x, +y, top) corner, three radii, at G2.
fn mixed_corner(options: &str) -> String {
    format!(
        r#"let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
let corner = [0.5, 0.5, 1.0];
block.fillet_each([
    [edges().at_vertex(corner).direction("x"), 0.05],
    [edges().at_vertex(corner).direction("y"), 0.2],
    [edges().at_vertex(corner).direction("z"), 0.3],
], #{{ continuity: "g2", label: "corner"{options} }});
m
"#
    )
}

fn corners(model: &popkern_script::Model) -> Vec<serde_json::Value> {
    let value = report_value(model);
    let row = value["features"]
        .as_array()
        .expect("rows")
        .iter()
        .find(|r| r["label"] == "corner")
        .cloned()
        .expect("corner row");
    row["blend"]["corners"].as_array().expect("corners").clone()
}

fn refused(script: &str) -> String {
    match ScriptHost::new().run(script, "refused.rhai") {
        Err(error) => error.message,
        Ok(_) => panic!("must be refused:\n{script}"),
    }
}

#[test]
fn the_report_names_the_fitted_cap_and_the_option_round_trips() {
    let outcome = ScriptHost::new()
        .run(
            &mixed_corner(r#", corner_cap: "fitted", corner_cap_strict: true"#),
            "fitted.rhai",
        )
        .unwrap_or_else(|e| panic!("{e}"));
    let mut model = outcome.model.expect("model");
    let served = corners(&model);
    assert_eq!(served.len(), 1, "{served:?}");
    assert_eq!(served[0]["construction"], "fitted", "{served:?}");
    assert_eq!(served[0]["sides"], 3);
    let at = served[0]["centroid"].as_array().expect("centroid");
    assert!(at[2].as_f64().expect("z") > 0.6, "{at:?}");

    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        emitted.contains(r#"corner_cap: "fitted""#) && emitted.contains("corner_cap_strict: true"),
        "{emitted}"
    );
}

#[test]
fn the_default_ladder_is_not_spelled_back() {
    let mut model = ScriptHost::new()
        .run(&mixed_corner(""), "auto.rhai")
        .unwrap_or_else(|e| panic!("{e}"))
        .model
        .expect("model");
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(!emitted.contains("corner_cap"), "{emitted}");
}

#[test]
fn a_mitred_corner_has_no_cap_to_choose() {
    let message = refused(&mixed_corner(r#", corner: "mitre", corner_cap: "fitted""#));
    assert!(message.contains("no cap to choose"), "{message}");
}

#[test]
fn strictness_needs_a_construction_that_can_decline() {
    let message = refused(&mixed_corner(", corner_cap_strict: true"));
    assert!(message.contains("\"setback\" or \"fitted\""), "{message}");
}

#[test]
fn an_unknown_corner_cap_is_refused_with_the_spellings() {
    let message = refused(&mixed_corner(r#", corner_cap: "sphere""#));
    assert!(
        message.contains("auto | setback | fitted | classic"),
        "{message}"
    );
}

#[test]
fn a_strict_setback_cap_that_folds_fails_the_feature() {
    let script = mixed_corner(r#", corner_cap: "setback", corner_cap_strict: true"#);
    match ScriptHost::new().run(&script, "strict.rhai") {
        Err(error) => assert!(
            error.message.contains("setback") && error.message.contains("fold"),
            "{}",
            error.message
        ),
        Ok(outcome) => {
            let value = report_value(&outcome.model.expect("model"));
            panic!("a strict setback decline must fail the feature: {value}");
        }
    }
}
