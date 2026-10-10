//! PK-404: the canonical emitter is a fixed point and the round trip
//! is bit-identical — boss plate recorded → emitted → run → saved, the
//! two documents replay to the same body; migration verbs and
//! configurations spelled from Rhai match the Rust recorder; features
//! the skin cannot spell decline by name.

use popkern_part::prelude::*;
use popkern_rhai::{ScriptHost, to_rhai};
use popkern_script::prelude::*;

const BOSS_PLATE: &str = include_str!("../examples/boss_plate.rhai");

fn signature(part: &Part) -> (usize, usize, Vec<[i64; 3]>) {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut verts: Vec<[i64; 3]> = mesh
        .vertices
        .iter()
        .map(|v| {
            let p = v.position;
            [
                (p.x * 1e9).round() as i64,
                (p.y * 1e9).round() as i64,
                (p.z * 1e9).round() as i64,
            ]
        })
        .collect();
    verts.sort();
    (part.face_count(), part.edge_count(), verts)
}

fn positions(part: &Part) -> Vec<[f64; 3]> {
    let mut v: Vec<[f64; 3]> = part
        .brep()
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    v.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    v
}

fn run(host: &mut ScriptHost, source: &str, name: &str) -> Model {
    host.run(source, name)
        .unwrap_or_else(|e| panic!("{name}: {}\n{}", e.render(), source))
        .model
        .expect("model")
}

fn outcomes(m: &Model) -> Vec<(String, String, bool)> {
    m.features()
        .iter()
        .map(|f| {
            let suppressed = m
                .tree()
                .feature(f.id)
                .map(|n| n.suppressed)
                .unwrap_or(false);
            (f.label.clone(), f.operation.clone(), suppressed)
        })
        .collect()
}

#[test]
fn boss_plate_emission_is_a_fixed_point_and_replays_bit_identical() {
    let mut host = ScriptHost::new();
    let mut original = run(&mut host, BOSS_PLATE, "boss_plate.rhai");
    let emitted = to_rhai(&mut original).expect("emit");
    eprintln!("{emitted}");
    let mut again = run(&mut host, &emitted, "emitted.rhai");
    let twice = to_rhai(&mut again).expect("emit again");
    assert_eq!(emitted, twice, "emit(run(emit(m))) == emit(m)");
    assert_eq!(outcomes(&again), outcomes(&original));
    let a = original.to_part().expect("part");
    let b = again.to_part().expect("part");
    assert_eq!(signature(&a), signature(&b));
    assert_eq!(
        positions(&a),
        positions(&b),
        "vertex positions bit-identical"
    );
    // Both documents save to the same wire form modulo nothing: labels,
    // parameters and selections all survived.
    let names = |m: &Model| {
        let mut v: Vec<(String, f64)> = m
            .params()
            .iter()
            .map(|p| (p.name.clone(), p.value))
            .collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    };
    assert_eq!(names(&original), names(&again));
    for line in [
        "m.profile(xy())",
        ".extrude(blind(0.4), #{ label: \"extrude-1\" })",
        "edges().between(\"cylinder\", \"plane\")",
        "#{ label: \"rounds\" }",
    ] {
        assert!(emitted.contains(line), "missing {line:?} in\n{emitted}");
    }
}

#[test]
fn migration_script_matches_the_rust_recorder_and_survives_emission() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("emit_round_trip");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("boss_plate.pkn");
    let mut host = ScriptHost::new();
    let original = run(&mut host, BOSS_PLATE, "boss_plate.rhai");
    original.save(&path).expect("save");
    // §12 migration, as a GUI would emit it.
    let migration = format!(
        r#"let m = open({path:?});
m.feature("rounds").modify(#{{ radius: 0.1 }});
m.set_param("extrude-1.distance", 0.5);
let report = m.replay();
if report.failures.len() > 0 {{ throw "failed features"; }}
m.save({path:?});
m"#,
        path = path.display()
    );
    let mut migrated = run(&mut host, &migration, "migrate.rhai");
    // The same migration through the Rust skin.
    let mut rust = Model::open(&path).expect("open");
    rust.replay().expect("replay");
    assert_eq!(
        signature(&migrated.to_part().expect("part")),
        signature(&rust.to_part().expect("part"))
    );
    let mut direct = Model::new();
    {
        let plate = direct
            .profile(popkern_script::Plane::XY)
            .rect(2.4, 1.6)
            .sketch("plate_sketch")
            .expect("s")
            .extrude(&mut direct, End::Blind(0.5.into()))
            .expect("e");
        let boss = direct
            .tool(|m| {
                m.profile(popkern_script::Plane::xy_at(0.2))
                    .circle(0.35)
                    .sketch("boss_sketch")?
                    .extrude(m, End::Blind(0.8.into()))
            })
            .expect("boss");
        let plate = plate.union(&mut direct, boss).expect("union");
        plate
            .fillet_continuity(
                &mut direct,
                edges()
                    .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                    .at_height(0.5)
                    .or(edges().straight().at_height(0.5))
                    .or(edges().straight().direction(Dir::Z)),
                0.1,
                ContinuityKind::G1,
                "rounds",
            )
            .expect("fillet");
    }
    assert_eq!(
        migrated.to_part().expect("part").face_count(),
        direct.to_part().expect("part").face_count()
    );
    // Emission of the migrated document carries the new values and
    // still round-trips.
    let emitted = to_rhai(&mut migrated).expect("emit");
    assert!(emitted.contains("blind(0.5)"), "{emitted}");
    assert!(
        emitted.contains(", 0.1, #{ label: \"rounds\" }"),
        "{emitted}"
    );
    let mut again = run(&mut host, &emitted, "emitted.rhai");
    assert_eq!(
        positions(&again.to_part().expect("part")),
        positions(&migrated.to_part().expect("part"))
    );
}

#[test]
fn every_timeline_verb_is_callable_and_replay_equivalent() {
    let mut host = ScriptHost::new();
    let body = BOSS_PLATE
        .trim_end()
        .strip_suffix("\nm")
        .expect("returns m");
    let source = format!(
        r#"{body}
let rounds = m.feature("rounds");
let chamfer = m.insert_before(rounds, |m, b| b.chamfer(edges().circular().at_height(1.0), 0.05, #{{ label: "top_bevel" }}));
print(chamfer.label);
print(m.features().len());
chamfer.suppress();
print(m.report().failures.len() == 0);
chamfer.unsuppress();
m.reorder(chamfer, m.features().len() - 1);
print(m.features()[m.features().len() - 1].label);
let deleted = chamfer.delete();
print(deleted);
print(m.feature("boolean-1").dependents().len());
rounds.suppress();
print(m.undo());
print(m.redo());
m.undo();
let before = m.rollback_to(rounds);
print(before.failures.len() == 0);
m.rollback_to_end();
m"#
    );
    let outcome = host
        .run(&source, "verbs.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    // The boss-plate script prints the fillet's describe line first.
    let mut printed = outcome.printed.into_iter().skip(1);
    assert_eq!(printed.next().as_deref(), Some("top_bevel"));
    assert_eq!(
        printed.next().as_deref(),
        Some("7"),
        "six boss-plate features plus the insert"
    );
    assert_eq!(printed.next().as_deref(), Some("true"));
    assert_eq!(printed.next().as_deref(), Some("top_bevel"));
    assert_eq!(printed.next().as_deref(), Some("1"));
    assert_eq!(printed.next().as_deref(), Some("1"));
    assert_eq!(printed.next().as_deref(), Some("true"));
    assert_eq!(printed.next().as_deref(), Some("true"));
    assert_eq!(printed.next().as_deref(), Some("true"));
    let mut model = outcome.model.expect("model");
    let mut reference = run(&mut ScriptHost::new(), BOSS_PLATE, "ref.rhai");
    assert_eq!(
        signature(&model.to_part().expect("part")),
        signature(&reference.to_part().expect("part")),
        "insert + delete + undo/redo net to the original document"
    );
}

#[test]
fn configurations_labels_and_suppression_survive_emission() {
    let mut host = ScriptHost::new();
    let source = r#"
let m = model();
m.param("size", 1.0, "m");
m.param_expr("height", "size * 0.75", "m");
let cube = m.profile(xy()).rect(p("size"), p("size")).sketch("base").extrude(blind(p("height")));
cube = cube.fillet(edges().straight().direction("z"), 0.1, #{ label: "posts" });
cube = cube.translate(0.0, 0.0, 0.25, #{ label: "lift" });
m.feature("lift").suppress();
m.label_body("block", cube);
m.label_edge("post_edge", edges().straight().direction("z").chord_length_between(0.7, 0.8).allow_empty());
m.add_config("tall", "double size");
m.set_override("tall", "size", 2.0);
m.set_suppression("tall", m.feature("posts"), true);
m
"#;
    let source = source.replace(
        "m.label_edge(\"post_edge\", edges().straight().direction(\"z\").chord_length_between(0.7, 0.8).allow_empty());",
        "m.label_face(\"top\", faces().facing(\"z\"));",
    );
    let mut original = run(&mut host, &source, "config.rhai");
    let emitted = to_rhai(&mut original).expect("emit");
    eprintln!("{emitted}");
    for line in [
        "m.param(\"size\", 1.0, \"m\");",
        "m.param_expr(\"height\", \"size * 0.75\", \"m\");",
        ".rect(p(\"size\"), p(\"size\"))",
        ".extrude(blind(p(\"height\"))",
        "m.feature(\"lift\").suppress();",
        "m.label_body(\"block\", b",
        "m.label_entity(\"top\", ",
        "m.add_config(\"tall\", \"double size\");",
        "m.set_override(\"tall\", \"size\", 2.0);",
        "m.set_suppression(\"tall\", m.feature(\"posts\"), true);",
    ] {
        assert!(emitted.contains(line), "missing {line:?} in\n{emitted}");
    }
    let mut again = run(&mut host, &emitted, "emitted.rhai");
    assert_eq!(to_rhai(&mut again).expect("emit"), emitted, "fixed point");
    assert_eq!(outcomes(&again), outcomes(&original));
    assert_eq!(again.labels(), original.labels(), "entity labels");
    assert_eq!(again.configs(), original.configs());
    assert_eq!(
        positions(&again.to_part().expect("part")),
        positions(&original.to_part().expect("part"))
    );
    again.set_active_config(Some("tall")).expect("tall");
    original.set_active_config(Some("tall")).expect("tall");
    assert_eq!(
        positions(&again.to_part().expect("part")),
        positions(&original.to_part().expect("part")),
        "configuration replays identically"
    );
}

#[test]
fn unspellable_features_decline_by_name() {
    let mut m = Model::new();
    let cube = m
        .profile(popkern_script::Plane::XY)
        .rect(1.0, 1.0)
        .sketch("cube")
        .expect("s")
        .extrude(&mut m, End::Blind(1.0.into()))
        .expect("e");
    let brep = m.to_brep().expect("brep");
    let _ = cube;
    m.raw_body(brep, SnapshotOrigin::Polyhedron)
        .expect("raw body");
    let err = to_rhai(&mut m).expect_err("raw_body has no spelling");
    match err {
        popkern_rhai::EmitError::NotApplicable {
            feature, operation, ..
        } => {
            assert_eq!(operation, "raw_body");
            assert_eq!(feature, "raw_body-1");
        }
    }
}

/// PK-754: a unit-suffixed literal in a parameter expression keeps its
/// spelling and dimension through emission; a bare number added to a
/// length is refused with a unit-suffix hint instead of silently
/// meaning metres.
#[test]
fn suffixed_literals_survive_emission_and_bare_sums_are_refused() {
    let mut host = ScriptHost::new();
    let source = r#"
let m = model();
m.param("size", 1.0, "m");
m.param_expr("height", "size - 250mm", "m");
let cube = m.profile(xy()).rect(p("size"), p("size")).sketch("base").extrude(blind(p("height")));
m
"#;
    let mut original = run(&mut host, source, "suffixed.rhai");
    let emitted = to_rhai(&mut original).expect("emit");
    assert!(
        emitted.contains("m.param_expr(\"height\", \"size - 250mm\", \"m\");"),
        "{emitted}"
    );
    let mut again = run(&mut host, &emitted, "emitted.rhai");
    assert_eq!(to_rhai(&mut again).expect("emit"), emitted, "fixed point");
    assert_eq!(
        positions(&again.to_part().expect("part")),
        positions(&original.to_part().expect("part"))
    );
    let bare = source.replace("size - 250mm", "size - 250");
    let err = host
        .run(&bare, "bare.rhai")
        .expect_err("bare literal in a sum");
    assert!(err.render().contains("unit suffix"), "{}", err.render());
}
