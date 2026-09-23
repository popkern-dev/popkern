//! Run-to-run determinism gate for multi-wrap thread serves
//! (PK-512).
//!
//! The cellular propagation planner's seed and audit choice decides
//! which fragments are ray-sampled and therefore which record-side
//! witness votes exist — when that choice rode HashMap iteration
//! order, the SAME len-110 M6 cut declined non-manifold three times
//! and then served (measured 2026-09-23): a thin, order-dependent
//! sign calibration overrode two correct wall verdicts. The planner
//! now iterates ordered structures; this gate pins run-to-run
//! identity of the served topology and volume by building the same
//! rod twice in one process (distinct HashMap seeds per instance, so
//! order dependence cannot hide).
//!
//! Release-only like the perf budgets (`--perf` tier): debug builds
//! self-skip.

use popkern_core::tolerance::Tolerance;
use popkern_part::fastener::{ThreadSpec, ThreadedRod};

fn build_counts(len: f64) -> (usize, usize, usize, usize, f64) {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let rod = ThreadedRod::new(spec, len)
        .tip_chamfer(false)
        .build(Tolerance::default())
        .expect("rod");
    assert_eq!(rod.open_edge_count().expect("mesh"), 0);
    let b = rod.brep();
    let mesh = rod
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let c = mesh.vertices[t[1]].position;
        let d = mesh.vertices[t[2]].position;
        v6 += a.x * (c.y * d.z - c.z * d.y) - a.y * (c.x * d.z - c.z * d.x)
            + a.z * (c.x * d.y - c.y * d.x);
    }
    (
        b.faces.len(),
        b.edges.len(),
        b.vertices.len(),
        b.loops.len(),
        (v6 / 6.0).abs(),
    )
}

#[test]
fn thread_serve_is_run_deterministic() {
    if cfg!(debug_assertions) {
        eprintln!("skipping thread determinism gate in debug (release tier)");
        return;
    }
    let a = build_counts(20.0);
    let b = build_counts(20.0);
    assert_eq!(
        (a.0, a.1, a.2, a.3),
        (b.0, b.1, b.2, b.3),
        "served topology differs between identical builds: {a:?} vs {b:?}"
    );
    assert!(
        (a.4 - b.4).abs() <= 1e-9 * a.4.max(1.0),
        "served volume differs between identical builds: {} vs {}",
        a.4,
        b.4
    );
}
