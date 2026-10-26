//! PK-956: generative design through the Rhai skin, the kernel suite's
//! coarse case — a 12 × 4 × 4 mm cantilever carved to half its volume
//! at 1 mm cells (seconds): the body serves valid and every face named,
//! the kept wall keeps its label's face as an ancestor, the report row
//! records the design, show-as-code spells `.generative(..)` and the
//! emitted script rebuilds the same body; bad maps decline by argument.
//! (The map selects geometrically: an emitted script declares its
//! labels after the features, so a `labeled()` selection does not
//! rebuild from show-as-code, whatever the feature — PK-962.)

use popkern_rhai::{ScriptHost, to_rhai};
use popkern_script::Model;
use serde_json::Value;

const CANTILEVER: &str = r#"let m = model();
let b = m.profile(yz()).rect(4.0, 4.0).sketch("block").extrude(blind(12.0));
m.label_face("wall", faces().facing_neg("x"));
m.label_face("tip", faces().facing("x"));
b = b.generative(#{
    supports: [#{ faces: faces().facing_neg("x"), kind: "fixed" }],
    loads: [#{ faces: faces().facing("x"), force: 10.0, direction: [0.0, 0.0, -1.0] }],
    keep: [faces().facing_neg("x"), faces().facing("x")],
    volume_fraction: 0.5,
    filter: 1.5,
    keep_depth: 1.0,
    cell: 1.0,
    build_direction: "+z",
    max_iterations: 150,
    label: "carve",
});
m.assert_watertight();
m
"#;

fn run(host: &mut ScriptHost, source: &str) -> Model {
    host.run(source, "cantilever.rhai")
        .unwrap_or_else(|e| panic!("{}\n{source}", e.render()))
        .model
        .unwrap_or_else(|| panic!("no model"))
}

/// The generative row of a model's report.
fn row(m: &mut Model) -> Value {
    let report = serde_json::to_value(m.report()).unwrap_or(Value::Null);
    report["features"]
        .as_array()
        .and_then(|f| f.iter().find(|r| r["operation"] == "generative"))
        .cloned()
        .unwrap_or(Value::Null)
}

#[test]
fn a_carved_cantilever_serves_named_and_round_trips_as_code() {
    let mut host = ScriptHost::new();
    let mut m = run(&mut host, CANTILEVER);
    let brep = m.to_brep().unwrap_or_else(|e| panic!("{e}"));
    brep.validate().unwrap_or_else(|e| panic!("{e}"));
    assert!(brep.faces.values().all(|f| f.name.is_some()));
    // The wall's label still resolves on the design: its patches are
    // the design space's wall, modified.
    let wall = m.entity_by_label("wall").unwrap_or_else(|e| panic!("{e}"));
    let modified = brep
        .faces
        .values()
        .filter_map(|f| brep.name_pool.get(f.name?))
        .filter(|n| {
            n.role.category == popkern_core::naming::RoleCategory::Modified
                && n.ancestors
                    .iter()
                    .any(|a| brep.name_pool.get(*a) == Some(&wall))
        })
        .count();
    assert!(modified >= 1, "the wall's patches, named from it");
    let r = row(&mut m);
    assert_eq!(r["outcome"], "Served", "{r}");
    let design = &r["generative"];
    assert!(design["faces"].as_u64().unwrap_or(0) > 6, "{r}");
    let kept = design["body_fraction"].as_f64().unwrap_or(0.0);
    assert!((0.35..0.65).contains(&kept), "kept {kept}");
    assert_eq!(
        design["unsupported_fraction"], 0.0,
        "printable along +z: {r}"
    );

    // Show-as-code spells the feature, and the emitted script rebuilds
    // the same body (the second build is the cache's).
    let emitted = to_rhai(&mut m).unwrap_or_else(|e| panic!("{e}"));
    assert!(emitted.contains(".generative(#{"), "{emitted}");
    assert!(emitted.contains("build_direction: \"+z\""), "{emitted}");
    let mut again = run(&mut host, &emitted);
    let rebuilt = again.to_brep().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(rebuilt.faces.len(), brep.faces.len());
    assert_eq!(row(&mut again)["generative"]["cached"], true);
}

#[test]
fn bad_maps_decline_naming_the_argument() {
    let mut host = ScriptHost::new();
    for (bad, argument) in [
        ("supports: [],", "supports"),
        ("volume_fraction: 0.5, volume: 1.0,", "generative"),
        ("build_direction: \"up\",", "build_direction"),
    ] {
        let source = CANTILEVER
            .replace(
                r#"supports: [#{ faces: faces().facing_neg("x"), kind: "fixed" }],"#,
                if bad.starts_with("supports") {
                    bad
                } else {
                    r#"supports: [#{ faces: faces().facing_neg("x"), kind: "fixed" }],"#
                },
            )
            .replace(
                "volume_fraction: 0.5,",
                if bad.starts_with("volume_fraction") {
                    bad
                } else {
                    "volume_fraction: 0.5,"
                },
            )
            .replace(
                "build_direction: \"+z\",",
                if bad.starts_with("build_direction") {
                    bad
                } else {
                    "build_direction: \"+z\","
                },
            );
        let error = match host.run(&source, "bad.rhai") {
            Ok(_) => panic!("{bad} served"),
            Err(e) => e,
        };
        let text = error.render();
        assert!(text.contains(argument), "{bad}: {text}");
    }
}
