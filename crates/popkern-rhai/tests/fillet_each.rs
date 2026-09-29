//! PK-602: `fillet_each` — one fillet feature with several radius
//! groups — and the `at_vertex` edge selector. A corner whose edges
//! carry different radii is built as one vertex blend; groups never
//! share an edge; the groups are document state that round-trips
//! through "show as code" and stays parametric.

use popkern_rhai::{ScriptHost, to_rhai};

/// A unit block and its (+x, +y, top) corner filleted with three radii.
const MIXED_CORNER: &str = r#"let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
let corner = [0.5, 0.5, 1.0];
block.fillet_each([
    [edges().at_vertex(corner).direction("x"), 0.05],
    [edges().at_vertex(corner).direction("y"), 0.2],
    [edges().at_vertex(corner).direction("z"), 0.3],
], #{ continuity: "g2", label: "mixed_corner" });
m
"#;

fn run(script: &str) -> popkern_script::Model {
    ScriptHost::new()
        .run(script, "fillet_each.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload))
        .model
        .expect("model")
}

#[test]
fn a_mixed_radius_corner_serves_as_one_feature_and_round_trips() {
    let mut model = run(MIXED_CORNER);
    let part = model.to_part().expect("part");
    let brep = part.brep();
    assert!(brep.validate().is_ok(), "{:?}", brep.validate());
    // Three strips and the corner cap beyond the block's six faces.
    assert!(brep.faces.len() > 9, "faces: {}", brep.faces.len());

    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        emitted.contains(".fillet_each(["),
        "spelled back as fillet_each:\n{emitted}"
    );
    let mut again = run(&emitted);
    assert_eq!(
        again.to_part().expect("part").brep().faces.len(),
        brep.faces.len(),
        "the emitted script rebuilds the same body:\n{emitted}"
    );
}

#[test]
fn at_vertex_picks_exactly_the_three_corner_edges() {
    let one = r#"let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
block.fillet(edges().at_vertex([0.5, 0.5, 1.0]), 0.1, #{ label: "corner" });
m
"#;
    let model = run(one);
    let text = model
        .describe(&model.feature("corner").expect("feature"))
        .expect("describe");
    assert!(text.contains("on 3 edges"), "{text}");
}

#[test]
fn overlapping_groups_are_refused_with_the_groups_named() {
    let script = MIXED_CORNER.replace(
        "[edges().at_vertex(corner).direction(\"z\"), 0.3]",
        "[edges().at_vertex(corner), 0.3]",
    );
    let Err(error) = ScriptHost::new().run(&script, "overlap.rhai") else {
        panic!("overlapping groups must not run");
    };
    assert!(
        error.message.contains("groups 0 and 2 (1 edge)")
            && error.message.contains("groups 1 and 2 (1 edge)")
            && error.message.contains("one radius"),
        "the message names the groups: {}",
        error.message
    );
}

#[test]
fn a_radius_law_is_refused_on_groups() {
    let script = MIXED_CORNER.replace("#{ continuity", "#{ law: linear(0.1), continuity");
    let Err(error) = ScriptHost::new().run(&script, "law.rhai") else {
        panic!("a law with groups must not run");
    };
    assert!(
        error.message.contains("constant radii"),
        "{}",
        error.message
    );
}

#[test]
fn a_malformed_group_is_refused() {
    let script = MIXED_CORNER.replace("[edges().at_vertex(corner).direction(\"y\"), 0.2]", "[0.2]");
    let Err(error) = ScriptHost::new().run(&script, "shape.rhai") else {
        panic!("a malformed group must not run");
    };
    assert!(
        error
            .message
            .contains("group 1 must be [edge selector, radius]"),
        "{}",
        error.message
    );
}

#[test]
fn a_group_radius_is_a_parameter() {
    let mut model = run(MIXED_CORNER);
    let describe = model
        .describe(&model.feature("mixed_corner").expect("feature"))
        .expect("describe");
    assert!(
        describe.contains("3 radius groups (0.05, 0.2, 0.3)"),
        "{describe}"
    );
    model
        .set_param("mixed_corner.radius_1", 0.15)
        .expect("set group radius");
    let describe = model
        .describe(&model.feature("mixed_corner").expect("feature"))
        .expect("describe");
    assert!(describe.contains("(0.05, 0.15, 0.3)"), "{describe}");
    let part = model.to_part().expect("replays with the new group radius");
    assert!(part.brep().validate().is_ok());
}
