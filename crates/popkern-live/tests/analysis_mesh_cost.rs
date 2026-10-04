//! PK-577: the zebra view's 6° analysis mesh against the 12° display
//! mesh on representative gallery parts (the whole gallery is the
//! `analysis_mesh_cost` example): triangles, mesh time and buffer bytes
//! (positions + normals + indices + face ids). The zebra view pays this
//! once per served model, on entering zebra mode. The analysis mesh is
//! never coarser than the display mesh.

use popkern_live::{HostOptions, Session};
use std::time::Instant;

fn bytes(bodies: &[popkern_live::envelope::BodyMesh]) -> usize {
    bodies
        .iter()
        .map(|b| 4 * (b.positions.len() + b.normals.len() + b.indices.len() + b.face_ids.len()))
        .sum()
}

/// Run gallery part `name`, then build its analysis mesh; print the
/// cost row and assert the analysis mesh is never coarser.
fn cost(dir: &std::path::Path, name: &str) {
    let source = std::fs::read_to_string(dir.join(format!("{name}.rhai"))).expect("script");
    let mut session = Session::new(HostOptions::default());
    let t0 = Instant::now();
    let envelope = session.run(&source, name);
    let run_ms = t0.elapsed().as_secs_f64() * 1e3;
    assert!(envelope.outcome.is_served(), "{name}");
    let display_tris: usize = envelope.bodies.iter().map(|b| b.indices.len() / 3).sum();
    let t1 = Instant::now();
    let analysis = match session.analysis_meshes() {
        Ok(analysis) => analysis,
        Err(outcome) => {
            eprintln!("COST {name}: display {display_tris} tris; analysis declined: {outcome:?}");
            return;
        }
    };
    let analysis_ms = t1.elapsed().as_secs_f64() * 1e3;
    let analysis_tris: usize = analysis.iter().map(|b| b.indices.len() / 3).sum();
    eprintln!(
        "COST {name}: display {display_tris} tris {:.1} KiB | analysis {analysis_tris} tris {:.1} KiB in {analysis_ms:.0} ms (x{:.2} tris; run incl. display mesh {run_ms:.0} ms)",
        bytes(&envelope.bodies) as f64 / 1024.0,
        bytes(&analysis) as f64 / 1024.0,
        analysis_tris as f64 / display_tris.max(1) as f64
    );
    assert_eq!(analysis.len(), envelope.bodies.len(), "{name}");
    assert!(
        analysis_tris >= display_tris,
        "{name}: analysis mesh coarser than display"
    );
}

fn gallery() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/live/examples")
}

/// Five fast, representative parts (curved NURBS, sphere, torus rounds,
/// setback apex).
#[test]
fn analysis_mesh_is_finer_on_representative_parts() {
    for name in [
        "pyramid_apex",
        "revolved_knob",
        "ball_bearing",
        "boss_plate",
        "bored_block",
    ] {
        cost(&gallery(), name);
    }
}

/// A model whose estimated analysis mesh exceeds the budget declines
/// typed, and the zebra view keeps the display mesh. Four 6.5-turn
/// springs carry over 100k display triangles, so their estimated
/// analysis mesh is over the 400k budget. (coil_spring's three springs,
/// 94.7k display triangles, fit since PK-589: 361k analysis triangles
/// in about 7 s native.)
#[test]
fn an_over_budget_analysis_mesh_declines_typed() {
    let source = r#"let m = model();
for i in 0..4 {
    m.profile(xy()).circle(0.04).sketch(`wire_${i}`)
        .helical_sweep(axis_at([i.to_float(), 0.0, 0.0], [0.0, 0.0, 1.0]),
                       0.3, 0.15, 6.5, #{ label: `coil_${i}` });
}
m
"#;
    let mut session = Session::new(HostOptions::default());
    let envelope = session.run(source, "springs");
    assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
    let display: usize = envelope.bodies.iter().map(|b| b.indices.len() / 3).sum();
    assert!(
        4 * display > popkern_live::mesh_buffers::ANALYSIS_TRIANGLE_BUDGET,
        "fixture must be over budget: {display} display triangles"
    );
    match session.analysis_meshes() {
        Err(popkern_live::Outcome::GateFailure { step, detail }) => {
            assert_eq!(step, "analysis_mesh");
            assert!(detail.contains("budget"), "{detail}");
        }
        other => panic!("expected a typed decline, got {:?}", other.map(|b| b.len())),
    }
}
