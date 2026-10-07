//! PK-700: the hinge verbs from Rhai — the returned map, option
//! errors, and document-as-code: a hinge recorded whole, one recorded
//! piecewise (`only` + `hinge_component`) and a living hinge emit, run
//! and replay to the same bodies, and the emission is a fixed point.

use popkern_rhai::{ScriptHost, to_rhai};
use popkern_script::prelude::*;

fn run(source: &str) -> Model {
    ScriptHost::new()
        .run(source, "hinge.rhai")
        .unwrap_or_else(|e| panic!("{}\n{source}", e.render()))
        .model
        .expect("model")
}

/// Each body's vertex positions, sorted: geometry, independent of
/// arena order.
fn geometry(m: &mut Model) -> Vec<Vec<[f64; 3]>> {
    m.bodies()
        .into_iter()
        .map(|b| {
            let mut v: Vec<[f64; 3]> = m
                .brep_of(b)
                .expect("brep")
                .vertices
                .values()
                .map(|v| [v.position.x, v.position.y, v.position.z])
                .collect();
            v.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
            v
        })
        .collect()
}

const SOURCE: &str = r#"let m = model();
let h = m.hinge("butt", 60.0, #{ knuckles: 7, clearance: 0.2, screw: "M3", label: "door" });
let lid = m.hinge("print_in_place", 40.0, #{ only: "leaf_b", angle: 120.0, label: "lid" });
let lid_a = m.hinge_component("lid", "leaf_a");
let flex = m.living_hinge("slotted", 12.0, 40.0, #{ slot_width: 0.4, label: "flex" });
let web = m.living_hinge("straight", 10.0, 30.0, #{ web: 0.3, label: "web" });
m
"#;

#[test]
fn hinges_emit_as_code_that_rebuilds_the_same_bodies() {
    let mut original = run(SOURCE);
    assert_eq!(original.bodies().len(), 3 + 2 + 2);
    let emitted = to_rhai(&mut original).expect("emit");
    assert!(emitted.contains("m.hinge(\"butt\", 60"), "{emitted}");
    assert!(
        emitted.contains("m.hinge_component(\"door\", \"leaf_b\""),
        "{emitted}"
    );
    assert!(
        emitted.contains("m.hinge_component(\"lid\", \"leaf_a\""),
        "{emitted}"
    );
    // Defaults that follow the length are not spelled out.
    assert!(!emitted.contains("leaf_width"), "{emitted}");
    let mut again = run(&emitted);
    assert_eq!(geometry(&mut original), geometry(&mut again), "{emitted}");
    // The emitted document emits itself (fixed point).
    assert_eq!(to_rhai(&mut again).expect("emit again"), emitted);
    // The parameter tables match too: one set per hinge.
    let names = |m: &Model| {
        let mut v: Vec<(String, Option<String>)> = m
            .params()
            .into_iter()
            .map(|p| (p.name, p.expression))
            .collect();
        v.sort();
        v
    };
    assert_eq!(names(&original), names(&again));
}

#[test]
fn the_map_and_its_errors() {
    let m = run(r#"let m = model();
let h = m.hinge("print_in_place", 50.0);
if h.pin != () { throw "print-in-place has no pin"; }
if !h.leaf_a.is_alive() || !h.leaf_b.is_alive() { throw "leaves"; }
m
"#);
    assert_eq!(m.bodies().len(), 2);
    for (source, needle) in [
        ("let m = model(); m.hinge(\"hinged\", 50.0); m", "\"butt\""),
        (
            "let m = model(); m.hinge(\"butt\", 50.0, #{ holes: 2 }); m",
            "need a screw",
        ),
        (
            "let m = model(); m.hinge(\"butt\", 50.0, #{ only: \"lid\" }); m",
            "\"leaf_a\"",
        ),
        (
            "let m = model(); m.hinge(\"butt\", 50.0, #{ wobble: 1 }); m",
            "no option 'wobble'",
        ),
        (
            "let m = model(); m.living_hinge(\"straight\", 10.0, 30.0, #{ slot_pitch: 1.0 }); m",
            "takes web, root_radius",
        ),
        (
            "let m = model(); m.hinge_component(\"nobody\", \"leaf_b\"); m",
            "no hinge labelled",
        ),
    ] {
        let err = ScriptHost::new()
            .run(source, "bad.rhai")
            .err()
            .unwrap_or_else(|| panic!("{source} should fail"));
        let text = err.render();
        assert!(text.contains(needle), "{source}: {text}");
    }
}
