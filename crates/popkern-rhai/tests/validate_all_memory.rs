//! `validate_all` retains memory per call, not per call count (PK-829).
//!
//! An arm build with `validate_all` after every boolean reached 16.5 GB
//! resident over 16 rayon workers (1.6 GB on one). Attributed with a
//! counting allocator (a probe; the workspace forbids `unsafe`): over
//! 60 calls on an M10 hex bolt the live heap returns to its baseline
//! after every call (2.7 -> 2.9 MB), each call's transient peak is
//! 27.9 MB, and resident memory plateaus at 138 MB from the tenth call
//! on — glibc's per-thread arenas keep the freed transient pages (the
//! arena cap measured 16.5 -> 7.3 GB on the arm). Retention is bounded
//! by the largest transient peak times the workers, independent of how
//! many times validation runs. This pins it: over 40 calls on a
//! four-worker pool resident memory grows in at most a couple of steps
//! per worker, never call after call (a 5 MB leak per call fails it).

use popkern_rhai::ScriptHost;

/// Validation workers: arenas grow at most about once each.
const WORKERS: usize = 4;

fn resident_mb() -> Option<f64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with("VmRSS"))?;
    Some(line.split_whitespace().nth(1)?.parse::<f64>().ok()? / 1024.0)
}

#[test]
fn validate_all_memory_plateaus_with_call_count() {
    if resident_mb().is_none() {
        eprintln!("[pk829] SKIP: no /proc/self/status");
        return;
    }
    let source = include_str!("../examples/hex_bolt.rhai");
    let mut run = ScriptHost::new()
        .run(source, "hex_bolt.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    let brep = parts[0].brep().clone();
    // A fixed pool: each worker's arena grows once, so the bound below
    // is per worker (beside the whole suite on a many-core host the
    // global pool spread the work over more workers call by call, 15
    // steps in 40).
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(WORKERS)
        .build()
        .expect("pool");
    let mut at = Vec::new();
    for _ in 0..40 {
        assert!(pool.install(|| brep.validate_all()).is_empty());
        at.push(resident_mb().unwrap_or(0.0));
    }
    // A leak grows on (nearly) every call; arenas grow in a few steps,
    // the first time a worker holds a transient peak (measured series:
    // 82, 89, 89, 109 x 35, 131 x 2 MB — steps at calls 2, 4 and 39).
    let steps = at.windows(2).filter(|w| w[1] > w[0] + 1.0).count();
    eprintln!("[pk829] {steps} growth steps over 40 calls; resident {at:.0?} MB");
    assert!(
        steps <= 2 * WORKERS,
        "resident memory grew on {steps} of {} validate_all calls: {at:.0?}",
        at.len()
    );
}
