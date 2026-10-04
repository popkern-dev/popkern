//! PK-577 cost report: the zebra view's 6° analysis mesh against the
//! 12° display mesh over the whole gallery (triangles, mesh time and
//! buffer bytes). About 20 minutes, most of it the scripts themselves,
//! so it is an on-demand report, not a test; `tests/analysis_mesh_cost.rs`
//! gates five representative parts.
//!
//! `cargo run --release -p popkern-live --example analysis_mesh_cost [part…]`

use popkern_core::time::Stopwatch;
use popkern_live::{HostOptions, Session};

fn bytes(bodies: &[popkern_live::envelope::BodyMesh]) -> usize {
    bodies
        .iter()
        .map(|b| 4 * (b.positions.len() + b.normals.len() + b.indices.len() + b.face_ids.len()))
        .sum()
}

fn main() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/live/examples");
    let mut parts: Vec<String> = std::fs::read_dir(&dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.extension().is_some_and(|x| x == "rhai"))
                .filter_map(|p| p.file_stem().and_then(|s| s.to_str()).map(str::to_string))
                .collect()
        })
        .unwrap_or_default();
    parts.sort();
    // Optional part names on the command line narrow the report.
    let only: Vec<String> = std::env::args().skip(1).collect();
    if !only.is_empty() {
        parts.retain(|p| only.contains(p));
    }
    for name in parts {
        let Ok(source) = std::fs::read_to_string(dir.join(format!("{name}.rhai"))) else {
            println!("{name}: unreadable");
            continue;
        };
        let mut session = Session::new(HostOptions::default());
        let t0 = Stopwatch::start();
        let envelope = session.run(&source, &name);
        let run_ms = t0.elapsed().map_or(f64::NAN, |d| d.as_secs_f64() * 1e3);
        if !envelope.outcome.is_served() {
            println!("{name}: not served");
            continue;
        }
        let display_tris: usize = envelope.bodies.iter().map(|b| b.indices.len() / 3).sum();
        let t1 = Stopwatch::start();
        match session.analysis_meshes() {
            Ok(analysis) => {
                let analysis_ms = t1.elapsed().map_or(f64::NAN, |d| d.as_secs_f64() * 1e3);
                let analysis_tris: usize = analysis.iter().map(|b| b.indices.len() / 3).sum();
                println!(
                    "{name}: display {display_tris} tris {:.1} KiB | analysis {analysis_tris} tris {:.1} KiB in {analysis_ms:.0} ms (x{:.2} tris; run incl. display mesh {run_ms:.0} ms)",
                    bytes(&envelope.bodies) as f64 / 1024.0,
                    bytes(&analysis) as f64 / 1024.0,
                    analysis_tris as f64 / display_tris.max(1) as f64
                );
            }
            Err(outcome) => {
                println!("{name}: display {display_tris} tris; analysis declined: {outcome:?}")
            }
        }
    }
}
