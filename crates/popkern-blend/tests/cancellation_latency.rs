//! §17 cancellation-latency acceptance (PK-190).
//!
//! Requirements item 17 promises: *"Cancelable (token, clean
//! termination < 1s, model unchanged)"*. The pre-signaled functional
//! test in `lib_tests.rs` proves the token is honored at the first
//! poll; this test proves the **latency** half of the claim — a token
//! fired *mid-flight* on a heavy multi-edge operation terminates the
//! pipeline in under one second, without mutating the input model.
//!
//! Method: measure the uncancelled wall-clock `T` of a 12-edge G2 cube
//! fillet (the heaviest mainstream pipeline path — the full Gregory /
//! KP construction ladder). Then run the same operation in a worker
//! thread and fire the token at `T/3`, guaranteeing the signal lands
//! while the pipeline is executing. The measured latency is the time
//! from `cancel()` to the worker observing `BlendError::Cancelled`.
//!
//! The strict 1 s assertion runs in release builds only — the claim is
//! about the shipped product, and debug-build phase timings are
//! meaningless for it. CI enforces this via the release-mode `perf`
//! job (`.github/workflows/ci.yml`); local runs use
//! `tools/quality-check.sh --perf`.

use popkern_blend::cancellation::CancellationToken;
use popkern_blend::{
    BlendError, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
    fillet_edges_cancellable,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::{Duration, Instant};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The heaviest mainstream workload: all 12 cube edges through the G2
/// construction ladder (rational Gregory rungs + KP setback/polar
/// caps). Long enough in release (~seconds) that a `T/3` cancel is
/// guaranteed to land mid-pipeline.
fn heavy_workload() -> (BRep, Vec<FilletSpec>, FilletOptions) {
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    (cube, specs, options)
}

/// §17 latency gate: a cancellation fired mid-flight on a 12-edge G2
/// cube fillet must terminate the pipeline in < 1 s with the input
/// model untouched.
#[test]
fn mid_flight_cancellation_terminates_within_one_second() {
    if cfg!(debug_assertions) {
        eprintln!(
            "§17 cancellation latency: skipped in debug — the < 1 s claim is a \
             release-build product property, enforced by the CI perf job"
        );
        return;
    }

    let (cube, specs, options) = heavy_workload();
    let pre_faces = cube.faces.len();
    let pre_edges = cube.edges.len();
    let pre_vertices = cube.vertices.len();

    // Uncancelled reference run: calibrates where "mid-flight" is and
    // proves the workload is long enough for the cancel to land inside
    // the pipeline rather than after it.
    let t0 = Instant::now();
    let uncancelled = fillet_edges(&cube, &specs, &options, tol(), None);
    let full_run = t0.elapsed();
    assert!(
        uncancelled.is_ok(),
        "reference 12-edge G2 cube fillet must succeed for the latency \
         measurement to be meaningful: {uncancelled:?}"
    );
    assert!(
        full_run > Duration::from_millis(300),
        "workload too fast ({full_run:?}) for a mid-flight cancel — pick a \
         heavier fixture so the latency claim is actually exercised"
    );

    let token = CancellationToken::new();
    let worker_token = token.clone();
    let (worker_cube, worker_specs, worker_options) =
        (cube.clone(), specs.clone(), options.clone());
    let worker = std::thread::spawn(move || {
        fillet_edges_cancellable(
            &worker_cube,
            &worker_specs,
            &worker_options,
            Tolerance::default(),
            None,
            &worker_token,
        )
    });

    // Fire the token a third of the way through the reference runtime:
    // deep inside phase 2/3 construction, far from both the entry poll
    // and natural completion.
    std::thread::sleep(full_run / 3);
    let t_cancel = Instant::now();
    token.cancel();
    let result = worker.join().expect("worker thread must not panic");
    let latency = t_cancel.elapsed();

    eprintln!(
        "§17 cancellation latency: {latency:?} from cancel() to return \
         (uncancelled reference {full_run:?}, cancelled at {:?})",
        full_run / 3
    );

    assert!(
        matches!(result, Err(BlendError::Cancelled)),
        "mid-flight cancel must surface BlendError::Cancelled, got {result:?}"
    );
    assert!(
        latency < Duration::from_secs(1),
        "§17 violation: cancellation took {latency:?} to terminate \
         (required < 1 s) — a pipeline phase is running too long between \
         cancellation polls"
    );

    // Model unchanged: the input BRep is untouched by a cancelled run.
    assert_eq!(
        cube.faces.len(),
        pre_faces,
        "face count changed by cancelled run"
    );
    assert_eq!(
        cube.edges.len(),
        pre_edges,
        "edge count changed by cancelled run"
    );
    assert_eq!(
        cube.vertices.len(),
        pre_vertices,
        "vertex count changed by cancelled run"
    );
    assert!(
        cube.validate().is_ok(),
        "input BRep no longer validates after a cancelled run"
    );
}
