//! PK-276 sweep cancellation latency (acceptance row 412).
//!
//! Fixture-weight evidence (release, 16-core dev reference,
//! 2026-08-21): a 40-control-point serpentine sweep and a 12-turn
//! helical sweep both resolve in ~0.5 ms (as valid declines or
//! serves) — the 200-frame clamp bounds every sweep to
//! millisecond-class work, so the strict PK-190 "fire at T/3" pattern
//! cannot land mid-flight. The < 1 s claim is structural: sweep rides
//! the SAME shared detach layer whose mid-flight latency is gated on
//! the second-scale ops (boolean 4.9 s fixture, fillet 2.5 s,
//! smooth-mesh 3.9 s). This gate races a cancel against the run and
//! asserts the §17 bound on RETURN LATENCY after cancel(), whichever
//! side wins.

use popkern_core::cancellation::CancellationToken;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::progress::OperationControl;
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{SweepError, SweepOptions, SweepProfile, sweep_with_control};
use std::time::{Duration, Instant};

fn circle_profile(radius: f64) -> SweepProfile {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = radius;
    SweepProfile::Wire(NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, r, 0.0),
            Point3::new(0.0, r, 0.0),
            Point3::new(-r, r, 0.0),
            Point3::new(-r, 0.0, 0.0),
            Point3::new(-r, -r, 0.0),
            Point3::new(0.0, -r, 0.0),
            Point3::new(r, -r, 0.0),
            Point3::new(r, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    })
}

fn z_path(length: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, length)],
        weights: vec![1.0, 1.0],
    }
}

/// §17: a cancel fired against a running sweep returns within one
/// second — as `Cancelled`, or as the completed result if the worker
/// won the tie.
#[test]
fn sweep_cancel_returns_within_one_second_either_way() {
    let token = CancellationToken::new();
    let worker_token = token.clone();
    let worker = std::thread::spawn(move || {
        sweep_with_control(
            &circle_profile(1.0),
            &z_path(5.0),
            &SweepOptions::default(),
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
        "§17 sweep cancel: returned {latency:?} after cancel(), outcome {}",
        if result.is_ok() {
            "completed (won the tie)"
        } else {
            "Cancelled"
        }
    );
    assert!(
        result.is_ok() || matches!(result, Err(SweepError::Cancelled)),
        "cancel race must complete or surface Cancelled, got {:?}",
        result.err()
    );
    assert!(
        latency < Duration::from_secs(1),
        "§17 violation: sweep took {latency:?} to return after cancel()"
    );
}
