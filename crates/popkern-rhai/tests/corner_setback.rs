//! PK-623: the `setback` fillet option. It parses, round-trips through
//! "show as code", is refused where it cannot apply, and a setback the
//! corner cannot take fails the feature naming the range.

use popkern_rhai::{ScriptHost, to_rhai};

/// The unit block's (+x, +y, top) corner, three radii.
fn mixed_corner(options: &str) -> String {
    format!(
        r#"let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
let corner = [0.5, 0.5, 1.0];
block.fillet_each([
    [edges().at_vertex(corner).direction("x"), 0.05],
    [edges().at_vertex(corner).direction("y"), 0.2],
    [edges().at_vertex(corner).direction("z"), 0.3],
], #{{ label: "corner"{options} }});
m.assert_watertight();
m
"#
    )
}

fn refused(script: &str) -> String {
    match ScriptHost::new().run(script, "refused.rhai") {
        Err(error) => error.message,
        Ok(_) => panic!("must be refused:\n{script}"),
    }
}

#[test]
fn a_setback_serves_and_round_trips() {
    let outcome = ScriptHost::new()
        .run(&mixed_corner(", setback: 0.45"), "setback.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    let mut model = outcome.model.expect("model");
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(emitted.contains("setback: 0.45"), "{emitted}");
    // The automatic setback is not spelled back.
    let mut auto = ScriptHost::new()
        .run(&mixed_corner(""), "auto.rhai")
        .unwrap_or_else(|e| panic!("{e}"))
        .model
        .expect("model");
    assert!(!to_rhai(&mut auto).expect("emit").contains("setback"));
}

#[test]
fn a_setback_the_corner_cannot_take_names_the_range() {
    let message = refused(&mixed_corner(", setback: 0.001"));
    assert!(
        message.contains("setback 0.0010") && message.contains("more than"),
        "{message}"
    );
}

#[test]
fn a_setback_is_refused_where_it_cannot_apply() {
    let mitre = refused(&mixed_corner(r#", corner: "mitre", setback: 0.3"#));
    assert!(mitre.contains("no setback"), "{mitre}");
    let classic = refused(&mixed_corner(r#", corner_cap: "classic", setback: 0.3"#));
    assert!(classic.contains("do not set rounds back"), "{classic}");
    let negative = refused(&mixed_corner(", setback: -0.1"));
    assert!(negative.contains("positive distance"), "{negative}");
}
