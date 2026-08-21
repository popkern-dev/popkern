//! §15 boolean cancellation-latency acceptance (PK-276).
//!
//! `boolean_cancellable` historically polled the token only between
//! pipeline stages, so a cancel landing inside a long SSI stage waited
//! the stage out — the "< 1 s clean termination" claim held only while
//! no stage was long. PK-276 moved the entry point onto the shared
//! detach-on-cancel layer (`popkern_core::cancellation`), making the
//! latency structural. This is the PK-190-pattern gate proving it
//! mid-flight.
//!
//! Fixture: two overlapping unit spheres (offset 0.6/0.35/0.2) —
//! measured 4.89 s release best-of-1 (2026-08-21, 16-core dev
//! reference), comfortably past the 300 ms mid-flight floor. (The
//! fixture probe that selected it also surfaced two heavier failing
//! pairs, filed as PK-304 / PK-305.)
//!
//! Release-only: the < 1 s claim is a shipped-product property. CI
//! enforces via the release `perf` job; locally
//! `tools/quality-check.sh --perf`.

use popkern_bool::cancellation::CancellationToken;
use popkern_bool::{BoolError, BoolOperation, boolean, boolean_cancellable};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::{Duration, Instant};

fn heavy_operands() -> (BRep, BRep) {
    let tol = Tolerance::default();
    let a = BRep::sphere(1.0, tol).expect("sphere builds");
    let b = BRep::sphere(1.0, tol)
        .expect("sphere builds")
        .transform(&Transform3::translation(0.6, 0.35, 0.2))
        .expect("translation applies");
    (a, b)
}

/// §15 latency gate: a cancellation fired mid-flight on a
/// sphere-sphere union must terminate in < 1 s with the inputs
/// untouched.
#[test]
fn mid_flight_boolean_cancellation_terminates_within_one_second() {
    if cfg!(debug_assertions) {
        eprintln!(
            "§15 boolean cancellation latency: skipped in debug — the < 1 s claim \
             is a release-build product property, enforced by the CI perf job"
        );
        return;
    }
    let tol = Tolerance::default();
    let (a, b) = heavy_operands();
    let pre = (a.faces.len(), a.edges.len(), b.faces.len(), b.edges.len());

    // Uncancelled reference run: calibrates where "mid-flight" is.
    let t0 = Instant::now();
    let uncancelled = boolean(BoolOperation::Union, a.clone(), b.clone(), tol, None);
    let full_run = t0.elapsed();
    assert!(
        uncancelled.is_ok(),
        "reference sphere-sphere union must succeed: {uncancelled:?}"
    );
    assert!(
        full_run > Duration::from_millis(300),
        "workload too fast ({full_run:?}) for a mid-flight cancel — pick a \
         heavier fixture so the latency claim is actually exercised"
    );

    let token = CancellationToken::new();
    let worker_token = token.clone();
    let (worker_a, worker_b) = (a.clone(), b.clone());
    let worker = std::thread::spawn(move || {
        boolean_cancellable(
            BoolOperation::Union,
            &worker_a,
            &worker_b,
            tol,
            &worker_token,
        )
    });

    std::thread::sleep(full_run / 3);
    let t_cancel = Instant::now();
    token.cancel();
    let result = worker.join().expect("worker thread must not panic");
    let latency = t_cancel.elapsed();

    eprintln!(
        "§15 boolean cancellation latency: {latency:?} from cancel() to return \
         (uncancelled reference {full_run:?}, cancelled at {:?})",
        full_run / 3
    );
    assert!(
        matches!(result, Err(BoolError::Cancelled)),
        "mid-flight cancel must surface BoolError::Cancelled, got {result:?}"
    );
    assert!(
        latency < Duration::from_secs(1),
        "§15 violation: cancellation took {latency:?} to terminate (required \
         < 1 s) — the detach layer is not doing its job"
    );

    // Inputs untouched by the cancelled run.
    assert_eq!(
        (a.faces.len(), a.edges.len(), b.faces.len(), b.edges.len()),
        pre,
        "operand entity counts changed by a cancelled run"
    );
    assert!(a.validate().is_ok(), "operand A no longer validates");
    assert!(b.validate().is_ok(), "operand B no longer validates");
}
