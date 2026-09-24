//! The fillet `corner` option (PK-529): `"smooth"` (default) or
//! `"mitre"`, which trims every adjacent strip pair at a corner to its
//! mutual intersection (PK-522, PK-528). The option is document state:
//! it is recorded in the feature tree, survives "show as code", and a
//! mitre the kernel cannot build fails with its typed reason in the
//! payload.

use popkern_rhai::{ScriptHost, to_rhai};

const MITRED_BOX: &str = r#"let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
block.fillet(edges(), 0.1, #{ corner: "mitre", label: "rounds" });
m
"#;

#[test]
fn mitred_box_serves_and_round_trips() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(MITRED_BOX, "mitred_box.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let errors = part.brep().validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");

    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        emitted.contains("corner: \"mitre\""),
        "the corner mode is spelled back:\n{emitted}"
    );
    let mut again = host
        .run(&emitted, "emitted.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload))
        .model
        .expect("model");
    assert_eq!(
        again.to_part().expect("part").brep().faces.len(),
        part.brep().faces.len(),
        "the emitted script rebuilds the same mitred body"
    );
}

#[test]
fn smooth_is_the_default_and_is_not_spelled() {
    let script = MITRED_BOX.replace("corner: \"mitre\", ", "");
    let mut host = ScriptHost::new();
    let mut model = host
        .run(&script, "smooth_box.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload))
        .model
        .expect("model");
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        !emitted.contains("corner:"),
        "default stays implicit:\n{emitted}"
    );
}

#[test]
fn an_unknown_corner_mode_is_a_typed_argument_error() {
    let script = MITRED_BOX.replace("\"mitre\"", "\"bevel\"");
    let mut host = ScriptHost::new();
    let Err(error) = host.run(&script, "bad_corner.rhai") else {
        panic!("an unknown corner mode must not run");
    };
    assert!(
        error.message.contains("smooth") && error.message.contains("mitre"),
        "the error names the spellings: {}",
        error.message
    );
}

/// A mitre the kernel cannot build (G2 strips on a hexagonal prism:
/// crossings off the fitted strips, PK-531) fails the feature with the
/// mitre's typed reason in the payload.
#[test]
fn a_declined_mitre_surfaces_its_typed_reason() {
    let script = r#"let m = model();
let hex = m.profile(xy()).polygon([[1.0, 0.0], [0.5, 0.866], [-0.5, 0.866], [-1.0, 0.0], [-0.5, -0.866], [0.5, -0.866]]).sketch("hex").extrude(blind(1.0));
hex.fillet(edges(), 0.07, #{ corner: "mitre", continuity: "g2", label: "rounds" });
m
"#;
    let mut host = ScriptHost::new();
    let Err(error) = host.run(script, "hex_g2_mitre.rhai") else {
        panic!("PK-531 serves the G2 hexagonal mitre now: pick another declining fixture");
    };
    assert_eq!(error.kind, "FeatureFailed");
    assert!(
        error.payload.to_string().contains("mitre declined"),
        "typed mitre decline in the payload: {}",
        error.payload
    );
}
