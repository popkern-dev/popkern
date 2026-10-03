//! PK-633: labels and spelled names follow an entity through later
//! booleans — a face the boolean leaves whole (renamed `Modified`), a
//! tool face the boolean trims (minted `Generated`), and the same after
//! a second boolean — while no stable name changes.

use popkern_rhai::ScriptHost;

const SCRIPT: &str = r#"
let m = model();
let body = m.profile(xy_at(0.44)).rect(1.0, 1.0).sketch("plate").extrude(blind(0.12));
m.label_face("plate_side", faces().facing("x"));
let spelled_side = m.entity_spelling("plate_side");
let pin = m.tool(|m| {
    let p = m.profile(xy_at(0.02)).circle(0.07).at(0.3, 0.0).sketch("pin_0").extrude(blind(0.48));
    m.label_face("pin_0", faces().kind("cylinder"));
    p
});
body = body.union(pin);
let other = m.tool(|m| m.profile(xy_at(0.02)).circle(0.07).at(-0.3, 0.0).sketch("pin_1").extrude(blind(0.48)));
body = body.union(other);
for l in ["plate_side", "pin_0"] {
    m.label_face("now_" + l, faces().labeled(l));
    print(l + " -> " + m.entity_spelling("now_" + l));
}
// The spelling minted before both booleans still selects its face.
m.label_face("by_spelling", faces().named(spelled_side));
print("spelled -> " + m.entity_spelling("by_spelling"));
m
"#;

#[test]
fn labels_follow_faces_through_two_unions() {
    let outcome = ScriptHost::new()
        .run(SCRIPT, "lineage.rhai")
        .expect("labels resolve");
    let lines = outcome.printed;
    assert_eq!(lines.len(), 3, "{lines:?}");
    // The plate side survived both unions as one face; the pin's trimmed
    // cylinder is one face; the spelling resolves to the side.
    assert!(
        lines[0].starts_with("plate_side -> boolean-2/Modified/"),
        "{lines:?}"
    );
    assert!(lines[1].starts_with("pin_0 -> boolean-"), "{lines:?}");
    assert_eq!(
        lines[2].replace("spelled", "plate_side"),
        lines[0],
        "{lines:?}"
    );
}

#[test]
fn a_labelled_edge_the_boolean_splits_resolves_to_its_pieces() {
    let script = r#"
let m = model();
let plate = m.profile(xy()).rect(2.0, 1.0).sketch("plate").extrude(blind(0.2));
m.label_edge("front_top", edges().straight().direction("x").at_height(0.2).at_vertex([1.0, -0.5, 0.2]));
// A block standing over the middle of the front edge splits it in two.
let block = m.tool(|m| m.profile(xy_at(0.1)).rect(0.4, 0.4).at(0.0, -0.5).sketch("block").extrude(blind(0.5)));
plate.union(block);
m
"#;
    let mut outcome = ScriptHost::new().run(script, "split.rhai").expect("runs");
    let Some(popkern_assembly::Document::Part(mut model)) = outcome.take_document() else {
        panic!("a model")
    };
    let brep = model.to_brep().expect("body");
    let pieces = model
        .resolve_edges_in(
            &brep,
            &popkern_script::selection::edges().labeled("front_top"),
        )
        .expect("the label resolves");
    assert_eq!(
        pieces.names.len(),
        2,
        "both pieces of the split edge: {:?}",
        pieces.names
    );
}
