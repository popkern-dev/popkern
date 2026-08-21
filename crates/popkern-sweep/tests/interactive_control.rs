//! PK-276 sweep-family interactive control: pre-signalled functional
//! tests for every `*_with_control` entry, never-fired identity, and
//! progress-stream shape. The mid-flight latency gate lives in
//! `cancellation_latency.rs`.

use popkern_core::cancellation::CancellationToken;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::progress::{CaptureReporter, OperationControl};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{
    LoftOptions, LoftSection, PipeRadius, SweepError, SweepOptions, SweepProfile,
    guided_sweep_with_control, loft_with_control, pipe_with_control, revolve_with_control, sweep,
    sweep_with_control, swept_blend_with_control,
};

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

fn loft_sections() -> Vec<LoftSection> {
    [0.0_f64, 1.0]
        .iter()
        .map(|&z| {
            let SweepProfile::Wire(mut curve) = circle_profile(1.0) else {
                unreachable!("circle profile is a wire");
            };
            for cp in &mut curve.control_points {
                cp.z = z;
            }
            LoftSection {
                curve,
                position: None,
            }
        })
        .collect()
}

#[test]
fn pre_signalled_token_cancels_every_sweep_entry_point() {
    let token = CancellationToken::new();
    token.cancel();
    let cancelled = |r: Result<popkern_sweep::SweepResult, SweepError>, name: &str| {
        assert!(
            matches!(r, Err(SweepError::Cancelled)),
            "{name}: pre-signalled token must surface Cancelled"
        );
    };

    cancelled(
        sweep_with_control(
            &circle_profile(1.0),
            &z_path(5.0),
            &SweepOptions::default(),
            tol(),
            None,
            OperationControl::with_cancel(&token),
        ),
        "sweep",
    );
    cancelled(
        loft_with_control(
            &loft_sections(),
            &LoftOptions::default(),
            tol(),
            None,
            OperationControl::with_cancel(&token),
        ),
        "loft",
    );
    cancelled(
        pipe_with_control(
            &z_path(5.0),
            &PipeRadius::Constant(0.4),
            &SweepOptions::default(),
            tol(),
            None,
            OperationControl::with_cancel(&token),
        ),
        "pipe",
    );
    cancelled(
        guided_sweep_with_control(
            &circle_profile(1.0),
            &z_path(5.0),
            &[],
            &SweepOptions::default(),
            tol(),
            None,
            OperationControl::with_cancel(&token),
        ),
        "guided_sweep",
    );
    cancelled(
        swept_blend_with_control(
            &[(0.0, circle_profile(1.0)), (1.0, circle_profile(0.5))],
            &z_path(5.0),
            &SweepOptions::default(),
            tol(),
            None,
            OperationControl::with_cancel(&token),
        ),
        "swept_blend",
    );
    cancelled(
        revolve_with_control(
            &circle_profile(0.5),
            &Point3::new(3.0, 0.0, 0.0),
            &Vec3::new(0.0, 1.0, 0.0),
            std::f64::consts::TAU,
            tol(),
            None,
            OperationControl::with_cancel(&token),
        ),
        "revolve",
    );
}

/// Never-fired token + capture reporter: sweep result matches the
/// plain call and the progress stream has the promised shape.
#[test]
fn never_fired_token_matches_plain_sweep_and_reports_stages() {
    let profile = circle_profile(1.0);
    let path = z_path(5.0);
    let options = SweepOptions::default();
    let plain = sweep(&profile, &path, &options, tol(), None).expect("plain sweep serves");

    let token = CancellationToken::new();
    let mut reporter = CaptureReporter::default();
    let controlled = sweep_with_control(
        &profile,
        &path,
        &options,
        tol(),
        None,
        OperationControl::new(&token, &mut reporter),
    )
    .expect("controlled sweep serves");

    assert_eq!(
        body_fingerprint(&plain.brep),
        body_fingerprint(&controlled.brep),
        "controlled sweep must be fingerprint-identical to the plain run"
    );

    assert_eq!(
        reporter.events.len(),
        6,
        "sweep reports its 6 stages: {:?}",
        reporter.events.iter().map(|e| e.stage).collect::<Vec<_>>()
    );
    for (i, event) in reporter.events.iter().enumerate() {
        assert_eq!(event.stage_index, i);
        assert!(event.stage.starts_with("sweep."));
        assert_eq!(event.total_stages, 6);
    }
    assert_eq!(reporter.events[0].stage, "sweep.analyze_profile");
    assert_eq!(reporter.events[5].stage, "sweep.topology");
}
