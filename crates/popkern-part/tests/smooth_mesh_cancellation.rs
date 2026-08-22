//! PK-276 smooth-mesh interactive control (acceptance row 412 for
//! smooth mesh): pre-signalled functional test, never-fired identity,
//! progress shape, and the PK-190 mid-flight latency gate on the cube
//! cage (measured 3.86 s best-of-3 release at PK-190 close — long
//! enough that a T/3 cancel lands deep inside the spline pipeline).

use popkern_core::cancellation::CancellationToken;
use popkern_core::progress::{CaptureReporter, OperationControl};
use popkern_core::tolerance::Tolerance;
use popkern_part::{
    PartError, SmoothMeshOptions, smooth_solid_from_obj, smooth_solid_from_obj_with_control,
};
use std::time::{Duration, Instant};

/// Closed unit-cube control cage (8 vertices, 6 quads).
const CUBE_CAGE: &str = "\
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v  1  1  1
v -1  1  1
f 1 4 3 2
f 5 6 7 8
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";

fn tol() -> Tolerance {
    Tolerance::default()
}

/// FNV-1a over the serialized body — the workspace-standard
/// whole-body fingerprint (identity pin, PK-276 row 457).
fn body_fingerprint(brep: &popkern_core::topology::BRep) -> u64 {
    let bytes = serde_json::to_vec(brep).expect("output BRep must serialize");
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

#[test]
fn pre_signalled_token_cancels_smooth_mesh_before_any_work() {
    let token = CancellationToken::new();
    token.cancel();
    let result = smooth_solid_from_obj_with_control(
        CUBE_CAGE,
        tol(),
        SmoothMeshOptions::default(),
        OperationControl::with_cancel(&token),
    );
    assert!(
        matches!(result, Err(PartError::Cancelled)),
        "pre-signalled token must surface Cancelled, got {:?}",
        result.err()
    );
}

/// Never-fired token + capture reporter: result matches the plain
/// call, and the progress stream has the promised shape (monotone
/// 0-based indices, honest total, known labels).
#[test]
fn never_fired_token_matches_plain_smooth_mesh_and_reports_stages() {
    if !toms_tables_or_skip() {
        return;
    }
    let options = SmoothMeshOptions::default();
    let (plain_part, plain_report) =
        smooth_solid_from_obj(CUBE_CAGE, tol(), options).expect("plain cube cage serves");

    let token = CancellationToken::new();
    let mut reporter = CaptureReporter::default();
    let (part, report) = smooth_solid_from_obj_with_control(
        CUBE_CAGE,
        tol(),
        options,
        OperationControl::new(&token, &mut reporter),
    )
    .expect("cube cage serves under control");

    assert_eq!(
        body_fingerprint(part.brep()),
        body_fingerprint(plain_part.brep()),
        "controlled smooth mesh must be fingerprint-identical to the plain run"
    );
    assert_eq!(report.refinement_steps, plain_report.refinement_steps);
    assert_eq!(report.accepted_components, plain_report.accepted_components);

    // One serving attempt = 4 stages, indices monotone from 0.
    assert_eq!(
        reporter.events.len(),
        4,
        "direct-serving cage reports exactly one attempt's stages: {:?}",
        reporter.events.iter().map(|e| e.stage).collect::<Vec<_>>()
    );
    for (i, event) in reporter.events.iter().enumerate() {
        assert_eq!(event.stage_index, i);
        assert!(event.stage.starts_with("smooth_mesh."));
        assert!(event.total_stages >= reporter.events.len());
    }
    assert_eq!(reporter.events[0].stage, "smooth_mesh.spline_construction");
    assert_eq!(
        reporter.events[3].stage,
        "smooth_mesh.self_intersection_screen"
    );
}

/// §17 latency gate (PK-190 pattern): a token fired a third of the way
/// through the cube-cage construction terminates in < 1 s.
/// Release-only — enforced by the CI release perf job.
#[test]
fn mid_flight_smooth_mesh_cancellation_terminates_within_one_second() {
    if !toms_tables_or_skip() {
        return;
    }
    if cfg!(debug_assertions) {
        eprintln!(
            "§17 smooth-mesh cancellation latency: skipped in debug — release-build \
             product property, enforced by the CI perf job"
        );
        return;
    }
    let options = SmoothMeshOptions::default();

    let t0 = Instant::now();
    let uncancelled = smooth_solid_from_obj(CUBE_CAGE, tol(), options);
    let full_run = t0.elapsed();
    assert!(uncancelled.is_ok(), "reference run must serve");
    assert!(
        full_run > Duration::from_millis(300),
        "workload too fast ({full_run:?}) for a mid-flight cancel"
    );

    let token = CancellationToken::new();
    let worker_token = token.clone();
    // Boxed: PartError's largest variant trips clippy::result_large_err
    // when returned through the thread join.
    let worker = std::thread::spawn(move || {
        smooth_solid_from_obj_with_control(
            CUBE_CAGE,
            Tolerance::default(),
            options,
            OperationControl::with_cancel(&worker_token),
        )
        .map_err(Box::new)
    });
    std::thread::sleep(full_run / 3);
    let t_cancel = Instant::now();
    token.cancel();
    let result = worker.join().expect("worker must not panic");
    let latency = t_cancel.elapsed();

    eprintln!(
        "§17 smooth-mesh cancellation latency: {latency:?} from cancel() to return \
         (uncancelled reference {full_run:?})"
    );
    assert!(
        matches!(
            result.as_ref().err().map(|e| &**e),
            Some(PartError::Cancelled)
        ),
        "mid-flight cancel must surface PartError::Cancelled, got {:?}",
        result.err()
    );
    assert!(
        latency < Duration::from_secs(1),
        "§17 violation: smooth-mesh cancellation took {latency:?} (required < 1 s)"
    );
}

/// The PnS mask tables are optional (ACM ToMS package, untracked —
/// PK-314/PK-315). Tests that build smooth-mesh bodies need them at
/// runtime; on checkouts without the package they skip with a visible
/// notice instead of failing.
fn toms_tables_or_skip() -> bool {
    if popkern_core::geometry::polyhedral_spline::toms_tables_available() {
        return true;
    }
    eprintln!("[toms] ACM ToMS 1032 package absent — smooth-mesh test skipped (PK-314/PK-315)");
    false
}
