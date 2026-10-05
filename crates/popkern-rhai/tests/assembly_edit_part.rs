//! PK-678: `a.edit_part(name, |m| …)` edits a part definition with the
//! part verbs — every instance follows — and a closure that raises or
//! a part that does not exist leaves the assembly as it was.

use popkern_assembly::{Assembly, Document};
use popkern_rhai::ScriptHost;

const PINS: &str = r#"let a = assembly("pins");
let pin = a.part("pin", |m| { m.profile(xy()).circle(0.1).sketch("pin").extrude(blind(1.0)); });
a.place(pin, "left");
a.place(pin, "right").translate(0.5, 0.0, 0.0);
"#;

fn run(tail: &str) -> Result<Assembly, String> {
    let source = format!("{PINS}{tail}\na\n");
    let mut outcome = ScriptHost::new()
        .run(&source, "pins.rhai")
        .map_err(|e| e.to_string())?;
    match outcome.take_document() {
        Some(Document::Assembly(assembly)) => Ok(*assembly),
        _ => Err("not an assembly".into()),
    }
}

fn pin_labels(asm: &mut Assembly) -> Vec<String> {
    let id = asm.definition_named("pin").expect("pin");
    asm.part_model_mut(id)
        .expect("part")
        .features()
        .into_iter()
        .map(|f| f.label)
        .collect()
}

#[test]
fn an_inserted_feature_reaches_every_instance() {
    let mut asm = run(r#"a.edit_part("pin", |m| {
    m.insert_after(m.feature("extrude-1"), |m, body| body.chamfer(edges().circular().at_height(1.0), 0.02));
});"#)
    .expect("edit");
    assert_eq!(pin_labels(&mut asm), ["pin", "extrude-1", "chamfer-1"]);
    assert_eq!(asm.leaves().len(), 2, "both instances still placed");
    let flat = asm.flatten().expect("flatten");
    assert_eq!(flat.len(), 2);
    for body in &flat {
        assert!(
            body.brep.faces.len() > 3,
            "{}: the chamfer face is on every instance",
            body.path
        );
    }
}

#[test]
fn a_raising_closure_or_an_unknown_part_changes_nothing() {
    let mut caught = run(r#"try {
    a.edit_part("pin", |m| { m.active().chamfer(edges().circular().at_height(1.0), 0.02); throw "stop"; });
} catch (e) { }"#)
    .expect("caught");
    assert_eq!(
        pin_labels(&mut caught),
        ["pin", "extrude-1"],
        "the half-done edit is dropped"
    );
    let Err(unknown) = run(r#"a.edit_part("nope", |m| { });"#) else {
        panic!("an unknown part must raise");
    };
    assert!(unknown.contains("no part 'nope'"), "{unknown}");
}
