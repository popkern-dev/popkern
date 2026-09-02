//! Phase 3a acceptance (design §10): the boss plate authored in Rhai
//! is mesh-identical to the Rust recorder's build, in 30 lines or
//! fewer, and the recorded document carries the same labels and
//! parameters.

use popkern_part::prelude::*;
use popkern_rhai::ScriptHost;
use popkern_script::Plane;
use popkern_script::prelude::*;

const SCRIPT: &str = include_str!("../examples/boss_plate.rhai");

fn rust_recorded() -> Model {
    let mut m = Model::new();
    let plate = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate_sketch")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.4.into()))
        .expect("extrude");
    let boss = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.2))
                .circle(0.35)
                .sketch("boss_sketch")?
                .extrude(m, End::Blind(0.8.into()))
        })
        .expect("boss");
    let plate = plate.union(&mut m, boss).expect("union");
    plate
        .fillet_continuity(
            &mut m,
            edges()
                .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(0.4)
                .or(edges().straight().at_height(0.4))
                .or(edges().straight().direction(Dir::Z)),
            0.08,
            ContinuityKind::G1,
            "rounds",
        )
        .expect("fillet");
    m
}

fn mesh_signature(part: &Part) -> (usize, usize, Vec<[i64; 3]>) {
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
    (mesh.vertices.len(), mesh.triangles.len(), verts)
}

#[test]
fn boss_plate_in_rhai_is_mesh_identical_to_the_rust_recorder() {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(SCRIPT, "boss_plate.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut scripted = outcome.model.expect("the script returns its model");
    let mut reference = rust_recorded();
    let a = scripted.to_part().expect("part");
    let b = reference.to_part().expect("part");
    assert_eq!(a.face_count(), b.face_count());
    assert_eq!(a.edge_count(), b.edge_count());
    assert_eq!(mesh_signature(&a), mesh_signature(&b), "mesh identical");
    assert_eq!(scripted.report().failures().count(), 0);
    assert_eq!(scripted.report().downgraded_blends().count(), 0);
    let labels: Vec<String> = scripted
        .features()
        .iter()
        .map(|f| f.label.clone())
        .collect();
    let expected: Vec<String> = reference
        .features()
        .iter()
        .map(|f| f.label.clone())
        .collect();
    assert_eq!(labels, expected, "same timeline labels");
    let names = |m: &Model| {
        let mut names: Vec<String> = m.params().iter().map(|p| p.name.clone()).collect();
        names.sort();
        names
    };
    assert_eq!(names(&scripted), names(&reference), "same parameters");
    assert_eq!(
        outcome.printed,
        vec![
            reference
                .describe(&reference.feature("rounds").expect("label"))
                .expect("describe")
        ],
        "print reaches the host sink"
    );
    if let Some(elapsed) = outcome.elapsed {
        assert!(
            elapsed.as_secs_f64() <= 2.0,
            "boss plate must run in <= 2 s (spec §9), took {elapsed:?}"
        );
        eprintln!("[pk402] boss plate in Rhai ran in {elapsed:?}");
    }
}

#[test]
fn boss_plate_script_fits_the_thirty_line_budget() {
    let lines = SCRIPT
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with("//"))
        .count();
    assert!(lines <= 30, "boss plate script is {lines} lines");
}

#[test]
fn saved_document_reopens_and_replays_identically() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("boss_plate_rhai");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("boss_plate.pkn");
    let body = SCRIPT
        .trim_end()
        .strip_suffix("\nm")
        .expect("the example ends by returning its model");
    let source = format!(
        "{body}\nm.save({:?});\nm.to_step({:?});\nm.assert_step_roundtrip({:?});\nm\n",
        path.display(),
        dir.join("boss_plate.step").display(),
        dir.join("boss_plate_roundtrip.step").display()
    );
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&source, "boss_plate.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    let mut scripted = outcome.model.expect("model");
    let mut again = Model::open(&path).expect("open");
    again.replay().expect("replay");
    assert_eq!(
        mesh_signature(&again.to_part().expect("part")),
        mesh_signature(&scripted.to_part().expect("part"))
    );
    // The saved document is editable from Rhai too (campaign sentence).
    let migration = format!(
        "let m = open({:?}); m.set_param(\"extrude-1.distance\", 0.5); m.replay(); m.assert_watertight(); m",
        path.display()
    );
    let outcome = host
        .run(&migration, "migrate.rhai")
        .unwrap_or_else(|e| panic!("{e}"));
    let mut taller = outcome.model.expect("model");
    let z_max = taller
        .to_part()
        .expect("part")
        .brep()
        .vertices
        .values()
        .map(|v| v.position.z)
        .fold(f64::MIN, f64::max);
    assert!(
        (z_max - 1.0).abs() < 1e-9,
        "boss tops out at 1.0, got {z_max}"
    );
}
