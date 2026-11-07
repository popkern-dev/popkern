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
    // Every alive body: a separated piece is a body of its own.
    let part = m.to_part_all().expect("part");
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
    (
        "union_edge_contact",
        r#"let m = model();
let a = m.profile(xy()).rect(1.0, 1.0).sketch("a").extrude(blind(1.0));
let b = m.tool(|m| m.profile(xy()).rect(1.0, 1.0).at(1.0, 1.0).sketch("b").extrude(blind(1.0)));
a.union(b, #{ non_manifold: "split_to_manifold" });
m"#,
    ),
    (
        "pieces",
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate").extrude(blind(0.4));
let s = m.tool(|m| m.profile(xy_at(-0.5)).rect(0.4, 3.0).sketch("s").extrude(blind(2.0)));
plate.split([s]);
let right = plate.separate_pieces(containing([1.0, 0.0, 0.2]), #{ label: "right" });
plate.keep_pieces(largest(1));
right.fillet(edges().straight().direction("z"), 0.05);
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

#[test]
fn a_two_body_boolean_records_its_non_manifold_policy() {
    let mut host = ScriptHost::new();
    let source = SCRIPTS[5].1;
    let mut m = run(&mut host, source, "union_edge_contact");
    let emitted = to_rhai(&mut m).expect("emit");
    assert!(
        emitted.contains(r#"non_manifold: "split_to_manifold""#),
        "{emitted}"
    );
    let part = m.to_part().expect("part");
    assert!(
        !part.brep().topology_contacts.is_empty(),
        "contact recorded"
    );
    // The default declines the same edge contact.
    let declined = source.replace(r#", #{ non_manifold: "split_to_manifold" }"#, "");
    let error = host
        .run(&declined, "union_edge_contact_default")
        .expect_err("an edge contact declines by default");
    assert!(
        error.render().contains("non-manifold"),
        "{}",
        error.render()
    );
}

#[test]
fn piece_features_spell_their_rules_and_the_separated_body() {
    let mut host = ScriptHost::new();
    let mut m = run(&mut host, SCRIPTS[6].1, "pieces");
    let emitted = to_rhai(&mut m).expect("emit");
    assert!(
        emitted.contains(
            "let b2 = b0.separate_pieces(containing([1.0, 0.0, 0.2]), #{ label: \"right\" });"
        ),
        "{emitted}"
    );
    assert!(
        emitted.contains("b0 = b0.keep_pieces(largest(1)"),
        "{emitted}"
    );
    assert!(
        emitted.contains("b2 = b2.fillet("),
        "the fillet stays on the piece"
    );
    let alive = m.bodies().into_iter().filter(|b| b.is_alive(&m)).count();
    assert_eq!(alive, 2, "the plate and the separated piece");
}

#[test]
fn imprint_classes_answer_by_spelled_name() {
    let mut host = ScriptHost::new();
    let source = SCRIPTS[4].1;
    let mut m = run(&mut host, source, "imprint");
    let brep = m.to_part().expect("part").brep().clone();
    let mut asked = Vec::new();
    for face in brep.faces.values() {
        let name = brep
            .name_pool
            .get(face.name.expect("named"))
            .expect("pooled")
            .clone();
        let class = m.imprint_class(&name).expect("query").expect("classified");
        let spelled = m.spell_entity(&name).expect("spell");
        let want = match class {
            FaceClass::Inside => "inside",
            FaceClass::Outside => "outside",
            FaceClass::OnAligned => "on_aligned",
            FaceClass::OnOpposed => "on_opposed",
        };
        asked.push(format!(
            "if m.imprint_class({spelled:?}) != {want:?} {{ throw \"{spelled}\"; }}"
        ));
    }
    assert!(asked.iter().any(|a| a.contains("\"inside\"")), "{asked:?}");
    let script = source.replace("\nm", &format!("\n{}\nm", asked.join("\n")));
    let script = script.replacen(
        "\nm",
        "\nif m.imprint_class(\"extrude-1/Cap/9\") != () { throw \"unclassified\"; }\nm",
        1,
    );
    run(&mut host, &script, "imprint_classes");
}

/// Under `non_manifold: "general"` two boxes meeting along an edge
/// union into one general body with no contact left recorded; the
/// policy spells back as written and the emitted script replays to the
/// same body.
#[test]
fn a_general_union_spells_back_and_replays() {
    let mut host = ScriptHost::new();
    let source = SCRIPTS[5].1.replace("split_to_manifold", "general");
    let mut m = run(&mut host, &source, "union_general");
    let emitted = to_rhai(&mut m).expect("emit");
    assert!(emitted.contains(r#"non_manifold: "general""#), "{emitted}");
    let counts = |m: &mut Model| {
        let part = m.to_part().expect("part");
        let brep = part.brep();
        assert_eq!(brep.body_kind(), popkern_core::topology::BodyKind::General);
        assert!(brep.topology_contacts.is_empty());
        let errors = brep.validate_all();
        assert!(errors.is_empty(), "{errors:?}");
        (brep.vertices.len(), brep.edges.len(), brep.faces.len())
    };
    let first = counts(&mut m);
    let mut replayed = run(&mut host, &emitted, "union_general_replayed");
    assert_eq!(counts(&mut replayed), first);
}
