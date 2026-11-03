//! The partition's operations from Rhai (PK-1052, design: operations on
//! the cell complex, 5.4): recorded from a script, emitted as code,
//! re-run, the two documents replay to the same body with the same face
//! names; set expressions and piece rules spell back as written.

use popkern_rhai::{ScriptHost, to_rhai};
use popkern_script::prelude::*;

fn run(host: &mut ScriptHost, source: &str, name: &str) -> Model {
    host.run(source, name)
        .unwrap_or_else(|e| panic!("{name}: {}\n{}", e.render(), source))
        .model
        .expect("model")
}

/// Face count, every face named, the sorted face names as text, and
/// the volume to 1e-9.
fn signature(m: &mut Model) -> (usize, Vec<String>, i64) {
    let part = m.to_part().expect("part");
    let brep = part.brep();
    assert!(brep.validate().is_ok(), "{:?}", brep.validate());
    let mut names: Vec<String> = brep
        .faces
        .values()
        .map(|f| {
            let id = f.name.expect("every face is named");
            format!("{:?}", brep.name_pool.get(id))
        })
        .collect();
    names.sort();
    let tol = popkern_core::tolerance::Tolerance::default();
    let volume = popkern_measure::body_volume(brep, &tol).expect("volume");
    (brep.faces.len(), names, (volume * 1e9).round() as i64)
}

const SCRIPTS: &[(&str, &str)] = &[
    (
        "split",
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate").extrude(blind(0.4));
let a = m.tool(|m| m.profile(xy_at(-0.5)).rect(0.4, 3.0).at(-0.5, 0.0).sketch("a").extrude(blind(2.0)));
let b = m.tool(|m| m.profile(xy_at(-0.5)).rect(0.4, 3.0).at(0.5, 0.0).sketch("b").extrude(blind(2.0)));
plate.split([a, b], inside_tools(body(0).or(body(1))));
m"#,
    ),
    (
        "split_largest",
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate").extrude(blind(0.4));
let a = m.tool(|m| m.profile(xy_at(-0.5)).rect(0.4, 3.0).at(-0.5, 0.0).sketch("a").extrude(blind(2.0)));
plate.split([a], largest(1), #{ label: "keep_big" });
m"#,
    ),
    (
        "nary",
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate").extrude(blind(0.4));
let a = m.tool(|m| m.profile(xy_at(-0.5)).rect(0.8, 3.0).at(-0.3, 0.0).sketch("a").extrude(blind(2.0)));
let b = m.tool(|m| m.profile(xy_at(-0.5)).rect(0.8, 3.0).at(0.3, 0.0).sketch("b").extrude(blind(2.0)));
plate.nary_boolean([a, b], at_least(2));
m"#,
    ),
    (
        "xor",
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate").extrude(blind(0.4));
let s = m.tool(|m| m.profile(xy_at(-0.5)).circle(0.3).sketch("s").extrude(blind(2.0)));
plate.symmetric_difference(s);
m"#,
    ),
    (
        "imprint",
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate").extrude(blind(0.4));
let s = m.tool(|m| m.profile(xy_at(-0.5)).circle(0.3).sketch("s").extrude(blind(2.0)));
plate.imprint(s);
m"#,
    ),
];

#[test]
fn partition_features_emit_and_replay_the_same_body() {
    let mut host = ScriptHost::new();
    for (name, source) in SCRIPTS {
        let mut original = run(&mut host, source, name);
        let emitted = to_rhai(&mut original).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let mut again = run(&mut host, &emitted, name);
        let twice = to_rhai(&mut again).expect("emit again");
        assert_eq!(emitted, twice, "{name}: emit(run(emit(m))) == emit(m)");
        assert_eq!(
            signature(&mut original),
            signature(&mut again),
            "{name}:\n{emitted}"
        );
    }
}

#[test]
fn set_expressions_and_piece_rules_spell_back_as_written() {
    let mut host = ScriptHost::new();
    let mut m = run(&mut host, SCRIPTS[0].1, "split");
    let emitted = to_rhai(&mut m).expect("emit");
    assert!(
        emitted.contains("inside_tools(body(0).or(body(1)))"),
        "{emitted}"
    );
    let mut m = run(&mut host, SCRIPTS[3].1, "xor");
    let emitted = to_rhai(&mut m).expect("emit");
    assert!(emitted.contains(".symmetric_difference("), "{emitted}");
}
