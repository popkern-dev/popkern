//! Part appearance (six-axis-arm §4.5): `m.appearance` / `a.appearance`
//! set a definition's colour and finish; STEP carries the colour each
//! way (styled items owned by their own product), `.pkn` and
//! show-as-code carry all of it.

use popkern_assembly::{Assembly, Document};
use popkern_rhai::{ScriptHost, assembly_to_rhai};

const COLOURS: &str = r##"let a = assembly("colours");
let blue = a.part("housing", |m| { m.profile(xy()).circle(10.0).sketch("s").extrude(blind(20.0)); m.appearance(#{ color: "#2b5cab", metal: 0.6, rough: 0.35 }); });
let steel = a.part("bolt", |m| { m.profile(xy()).rect(4.0, 4.0).sketch("s").extrude(blind(12.0)); });
let plain = a.part("plain", |m| { m.profile(xy()).rect(2.0, 2.0).sketch("s").extrude(blind(2.0)); });
a.appearance(steel, #{ color: "steel", metal: 0.9, rough: 0.3 });
a.place(blue, "housing_1");
a.place(steel, "bolt_1").translate(20.0, 0.0, 0.0);
a.place(plain, "plain_1").translate(-20.0, 0.0, 0.0);
a
"##;

fn run(source: &str) -> Assembly {
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, "colours.rhai")
        .unwrap_or_else(|e| panic!("{}\n{source}", e.render()));
    match outcome.take_document() {
        Some(Document::Assembly(a)) => *a,
        _ => panic!("an assembly"),
    }
}

fn colour_of(asm: &mut Assembly, part: &str) -> Option<[f64; 3]> {
    let id = asm.definition_named(part).expect("definition");
    asm.part_model_mut(id)
        .expect("part")
        .appearance()
        .map(|a| a.color)
}

fn near(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() < 1e-6)
}

#[test]
fn appearance_sets_colour_and_finish() {
    let mut asm = run(COLOURS);
    let blue = [43.0 / 255.0, 92.0 / 255.0, 171.0 / 255.0];
    assert!(near(colour_of(&mut asm, "housing").expect("housing"), blue));
    let id = asm.definition_named("bolt").expect("bolt");
    let a = asm
        .part_model_mut(id)
        .expect("part")
        .appearance()
        .expect("bolt");
    assert!(near(a.color, [0.55, 0.56, 0.58]) && a.metal == 0.9 && a.rough == 0.3);
    assert!(colour_of(&mut asm, "plain").is_none());
}

#[test]
fn step_carries_each_part_its_own_colour() {
    let mut asm = run(COLOURS);
    let bytes = asm.to_step_bytes("colours").expect("step");
    let text = String::from_utf8_lossy(&bytes);
    assert_eq!(
        text.matches("STYLED_ITEM(").count(),
        2,
        "one styled solid per coloured part"
    );
    let (mut back, _) = Assembly::from_step_bytes(&bytes, "colours").expect("import");
    let blue = [43.0 / 255.0, 92.0 / 255.0, 171.0 / 255.0];
    // f32 in the document: compare to its precision.
    let close = |a: [f64; 3], b: [f64; 3]| (0..3).all(|i| (a[i] - b[i]).abs() < 1e-6);
    assert!(close(
        colour_of(&mut back, "housing").expect("housing colour"),
        blue
    ));
    assert!(close(
        colour_of(&mut back, "bolt").expect("bolt colour"),
        [0.55, 0.56, 0.58]
    ));
    assert!(
        colour_of(&mut back, "plain").is_none(),
        "an uncoloured part stays uncoloured"
    );
}

#[test]
fn pkn_and_show_as_code_keep_the_appearance() {
    let mut asm = run(COLOURS);
    let bytes = asm.to_pkn_bytes().expect("save");
    let Document::Assembly(mut back) = Document::from_pkn_bytes(&bytes).expect("open") else {
        panic!("assembly");
    };
    let id = back.definition_named("bolt").expect("bolt");
    let a = back
        .part_model_mut(id)
        .expect("part")
        .appearance()
        .expect("kept");
    assert_eq!((a.metal, a.rough), (0.9, 0.3));
    let code = assembly_to_rhai(&mut asm).expect("emit");
    assert_eq!(code.matches("m.appearance(").count(), 2, "{code}");
    let mut again = run(&code);
    assert_eq!(
        assembly_to_rhai(&mut again).expect("emit again"),
        code,
        "fixed point"
    );
}

#[test]
fn bad_appearances_are_argument_errors() {
    for (spelling, needle) in [
        (r##"#{ color: "#12345" }"##, "color"),
        (r##"#{ color: "chartreuse" }"##, "aluminium"),
        (r##"#{ color: [1.5, 0.0, 0.0] }"##, "red"),
        (r##"#{ color: "steel", metal: 2.0 }"##, "metal"),
        (r##"#{ metal: 0.5 }"##, "needs a color"),
    ] {
        let source = format!(
            "let m = model(); m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(1.0)); m.appearance({spelling}); m"
        );
        let err = ScriptHost::new()
            .run(&source, "bad.rhai")
            .expect_err("refused");
        assert!(
            err.to_string().contains(needle) || err.payload.to_string().contains(needle),
            "{spelling}: {err} / {}",
            err.payload
        );
    }
}
