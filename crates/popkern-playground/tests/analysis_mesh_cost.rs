//! PK-577: the zebra view's 6° analysis mesh against the 12° display
//! mesh, over every gallery part: triangles, mesh time and buffer bytes
//! (positions + normals + indices + face ids). The zebra view pays this
//! once per served model, on entering zebra mode. The analysis mesh is
//! never coarser than the display mesh.

use popkern_playground::{HostOptions, Session};
use std::time::Instant;

fn bytes(bodies: &[popkern_playground::envelope::BodyMesh]) -> usize {
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
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/playground/examples")
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

/// The whole gallery. About 20 minutes, most of it the scripts
/// themselves (fastener_kit alone runs for minutes), so it is a report
/// run on demand. The representative test above is the gate.
#[test]
#[ignore = "full-gallery cost report (~20 min): cargo test --release -p popkern-playground --test analysis_mesh_cost -- --ignored --nocapture"]
fn analysis_mesh_cost_over_the_gallery() {
    let dir = gallery();
    let mut parts: Vec<String> = std::fs::read_dir(&dir)
        .expect("gallery")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "rhai"))
        .filter_map(|p| p.file_stem().and_then(|s| s.to_str()).map(str::to_string))
        .collect();
    parts.sort();
    assert!(parts.len() >= 20, "gallery parts: {}", parts.len());
    for name in parts {
        cost(&dir, &name);
    }
}

/// A model whose estimated analysis mesh exceeds the budget declines
/// typed, and the zebra view keeps the display mesh (coil_spring: 94.6k
/// display triangles, an estimated 378k analysis triangles, and 50-107 s
/// to mesh natively).
#[test]
fn an_over_budget_analysis_mesh_declines_typed() {
    let source = std::fs::read_to_string(gallery().join("coil_spring.rhai")).expect("script");
    let mut session = Session::new(HostOptions::default());
    assert!(session.run(&source, "coil_spring").outcome.is_served());
    match session.analysis_meshes() {
        Err(popkern_playground::Outcome::GateFailure { step, detail }) => {
            assert_eq!(step, "analysis_mesh");
            assert!(detail.contains("budget"), "{detail}");
        }
        other => panic!("expected a typed decline, got {:?}", other.map(|b| b.len())),
    }
}
