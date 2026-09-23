//! PK-496 serve-cost probe: attribute the superlinear thread-cut
//! serve cost to boolean pipeline stages via the §16 telemetry
//! (`boolean_with_diagnostics` `PhaseMetric` rows). Measured at
//! PK-495 close (release, dev machine, partial contention): len 60 ≈
//! 4.5 min, 90 ≈ 9.8, 110 ≈ 18.1, 150 ≈ 26.5 — a 100-turn cut at
//! ~14 min is at the edge of interactive use.
//!
//! Run by hand:
//! `PK496_LEN=<len> cargo test -p popkern-fastener --release probe_thread_serve_cost -- --ignored --nocapture`

use popkern_core::tolerance::Tolerance;
use popkern_sketch::ProfileShape;

use crate::engine::{extrude_shape, helical_tool, plan_turns};
use crate::spec::ThreadSpec;
use crate::thread_form::external_groove;

#[test]
#[ignore = "heavy: PK-496 serve-cost telemetry probe — run by hand with PK496_LEN; feeds the cost gates"]
fn probe_thread_serve_cost() {
    let len: f64 = std::env::var("PK496_LEN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(20.0);
    let tolerance = Tolerance::default();
    let spec = ThreadSpec::parse("M6").expect("M6");
    let groove = external_groove(&spec).expect("groove");
    let shank = extrude_shape(
        ProfileShape::Circle {
            radius: spec.major_diameter() / 2.0,
        },
        0.0,
        len,
        "probe shank",
        tolerance,
    )
    .expect("shank");
    let (z0, turns) = plan_turns(0.0, len, spec.lead());
    let t_tool = popkern_core::time::Stopwatch::start();
    let tool = helical_tool(&spec, &groove, z0, turns, 0, "probe tool", tolerance).expect("tool");
    println!(
        "PK496_COST len={len} turns={turns} tool_build={:.1?} tool_faces={}",
        t_tool.elapsed(),
        tool.faces.len()
    );
    let t_bool = popkern_core::time::Stopwatch::start();
    let (result, diag) = popkern_bool::boolean_with_diagnostics(
        popkern_bool::BoolOperation::Subtraction,
        &shank,
        &tool,
        tolerance,
    )
    .expect("subtract");
    println!(
        "PK496_COST len={len} boolean_total={:.1?} faces={}",
        t_bool.elapsed(),
        result.faces.len()
    );
    for m in &diag.phase_metrics {
        println!(
            "PK496_PHASE len={len} {:?} {:.1?} entities={}",
            m.stage, m.duration, m.entity_count
        );
    }
}
