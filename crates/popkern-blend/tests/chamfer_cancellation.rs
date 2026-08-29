//! §17 chamfer cancellation (PK-276, acceptance row 458).
//!
//! Chamfers historically had NO cancellation: `chamfer_edges` passed
//! `None` into the apply phase (`lib_chamfer.rs`, "chamfer entry
//! points expose no cancellation token"). PK-276 wired the token
//! through the chamfer pipeline (entry poll, per-edge poll, chain
//! poll, and the apply phase's stage + per-corner polls) and added
//! `chamfer_edges_cancellable` on the shared detach layer.

use popkern_blend::cancellation::CancellationToken;
use popkern_blend::{
    BlendError, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges, chamfer_edges_cancellable,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::{Duration, Instant};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_chamfer_specs(cube: &BRep, count: usize) -> Vec<ChamferSpec> {
    cube.edges
        .keys()
        .take(count)
        .map(|e| ChamferSpec {
            edge: e,
            kind: ChamferKind::Equal(0.1),
        })
        .collect()
}

/// A pre-signalled token is honored at the entry poll: `Cancelled`
/// before any work, input untouched.
#[test]
fn pre_signalled_token_cancels_chamfer_before_any_work() {
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs = cube_chamfer_specs(&cube, 12);
    let token = CancellationToken::new();
    token.cancel();
    let result = chamfer_edges_cancellable(
        &cube,
        &specs,
        &ChamferOptions::default(),
        tol(),
        None,
        &token,
    );
    assert!(
        matches!(result, Err(BlendError::Cancelled)),
        "pre-signalled token must surface Cancelled, got {result:?}"
    );
    assert!(cube.validate().is_ok(), "input untouched by cancelled run");
}

/// Detach transparency: a cancellable call whose token never fires
/// produces a result equivalent to the plain call (same face/edge/
/// vertex counts, both validate) — the worker thread and clones are
/// invisible.
#[test]
fn never_fired_token_matches_plain_chamfer() {
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs = cube_chamfer_specs(&cube, 12);
    let plain = chamfer_edges(&cube, &specs, &ChamferOptions::default(), tol(), None)
        .expect("plain 12-edge chamfer serves");
    let token = CancellationToken::new();
    let detached = chamfer_edges_cancellable(
        &cube,
        &specs,
        &ChamferOptions::default(),
        tol(),
        None,
        &token,
    )
    .expect("cancellable 12-edge chamfer serves");
    assert_eq!(plain.brep.faces.len(), detached.brep.faces.len());
    assert_eq!(plain.brep.edges.len(), detached.brep.edges.len());
    assert_eq!(plain.brep.vertices.len(), detached.brep.vertices.len());
    assert!(detached.brep.validate().is_ok());
}

/// Mid-flight cancel returns within one second of the token firing —
/// the structural detach guarantee. A chamfer may be fast enough that
/// the worker finishes first (completion wins the tie, returning the
/// full result); either outcome is correct, but the RETURN LATENCY
/// after cancel() must stay under the §17 second either way.
#[test]
fn mid_flight_chamfer_cancel_returns_within_one_second() {
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs = cube_chamfer_specs(&cube, 12);
    let token = CancellationToken::new();
    let worker_token = token.clone();
    let (worker_cube, worker_specs) = (cube.clone(), specs.clone());
    let worker = std::thread::spawn(move || {
        chamfer_edges_cancellable(
            &worker_cube,
            &worker_specs,
            &ChamferOptions::default(),
            tol(),
            None,
            &worker_token,
        )
    });
    std::thread::sleep(Duration::from_millis(5));
    let t_cancel = Instant::now();
    token.cancel();
    let result = worker.join().expect("worker must not panic");
    let latency = t_cancel.elapsed();
    let outcome = match &result {
        Ok(_) => "completed (won the tie)".to_string(),
        Err(e) => format!("{e:?}"),
    };
    eprintln!("§17 chamfer cancel: returned {latency:?} after cancel(), outcome {outcome}");
    assert!(
        result.is_ok() || matches!(result, Err(BlendError::Cancelled)),
        "mid-flight cancel must complete or surface Cancelled, got {result:?}"
    );
    assert!(
        latency < Duration::from_secs(1),
        "§17 violation: chamfer took {latency:?} to return after cancel()"
    );
    assert!(cube.validate().is_ok(), "input untouched");
}
