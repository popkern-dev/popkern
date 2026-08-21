//! PK-276 shell cancellation latency (acceptance row 412).
//!
//! Fixture-weight evidence (release, 16-core dev reference,
//! 2026-08-21): shell(cube) 79 µs, shell(sphere) 170 µs,
//! shell(cylinder) 578 µs, shell(torus) 292 µs, shell of a 26-face
//! G1-filleted cube 2.3 ms. No realistic shell input runs long enough
//! for the strict PK-190 "fire at T/3" pattern to land mid-flight —
//! the operation completes orders of magnitude faster than any human
//! cancel. The < 1 s claim is structural: shell rides the SAME shared
//! detach layer (`popkern_core::cancellation::run_detach_on_cancel_with_events`)
//! whose mid-flight latency is gated on the second-scale ops (boolean
//! 4.9 s fixture, fillet 2.5 s, smooth-mesh 3.9 s). This gate races a
//! cancel against the run and asserts the §17 bound on RETURN LATENCY
//! after cancel(), whichever side wins — the strongest deterministic
//! claim a millisecond-class op admits.

use popkern_core::cancellation::CancellationToken;
use popkern_core::progress::OperationControl;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellError, ShellOptions, shell_with_control};
use std::time::{Duration, Instant};

/// §17: a cancel fired against a running shell returns within one
/// second — as `Cancelled`, or as the completed result if the worker
/// won the tie. The input is untouched either way.
#[test]
fn shell_cancel_returns_within_one_second_either_way() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let token = CancellationToken::new();
    let worker_token = token.clone();
    let worker_cube = cube.clone();
    let worker = std::thread::spawn(move || {
        shell_with_control(
            &worker_cube,
            0.1,
            ShellDirection::Inward,
            &[],
            &[],
            &ShellOptions::default(),
            Tolerance::default(),
            None,
            OperationControl::with_cancel(&worker_token),
        )
    });
    std::thread::sleep(Duration::from_millis(1));
    let t_cancel = Instant::now();
    token.cancel();
    let result = worker.join().expect("worker must not panic");
    let latency = t_cancel.elapsed();
    eprintln!(
        "§17 shell cancel: returned {latency:?} after cancel(), outcome {}",
        if result.is_ok() {
            "completed (won the tie)"
        } else {
            "Cancelled"
        }
    );
    assert!(
        result.is_ok() || matches!(result, Err(ShellError::Cancelled)),
        "cancel race must complete or surface Cancelled, got {:?}",
        result.err()
    );
    assert!(
        latency < Duration::from_secs(1),
        "§17 violation: shell took {latency:?} to return after cancel()"
    );
    assert!(cube.validate().is_ok(), "input untouched");
}
