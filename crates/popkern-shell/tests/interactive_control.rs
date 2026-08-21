//! PK-276 shell-family interactive control: pre-signalled functional
//! tests for every `*_with_control` entry, never-fired identity, and
//! progress-stream shape. The mid-flight latency gate lives in
//! `cancellation_latency.rs`.

use popkern_core::cancellation::CancellationToken;
use popkern_core::progress::{CaptureReporter, OperationControl};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{
    ShellDirection, ShellError, ShellOptions, offset_solid_with_control, shell,
    shell_variable_with_control, shell_with_control, thicken_with_control,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// FNV-1a over the serialized body — the workspace-standard
/// whole-body fingerprint (identity pin, PK-276 row 457).
fn body_fingerprint(brep: &BRep) -> u64 {
    let bytes = serde_json::to_vec(brep).expect("output BRep must serialize");
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

#[test]
fn pre_signalled_token_cancels_every_shell_entry_point() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let token = CancellationToken::new();
    token.cancel();
    let opts = ShellOptions::default();

    let shell_result = shell_with_control(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &opts,
        tol(),
        None,
        OperationControl::with_cancel(&token),
    );
    assert!(matches!(shell_result, Err(ShellError::Cancelled)));

    let variable_result = shell_variable_with_control(
        &cube,
        0.1,
        &[],
        ShellDirection::Inward,
        &[],
        &opts,
        tol(),
        None,
        OperationControl::with_cancel(&token),
    );
    assert!(matches!(variable_result, Err(ShellError::Cancelled)));

    let offset_result = offset_solid_with_control(
        &cube,
        0.1,
        &opts,
        tol(),
        OperationControl::with_cancel(&token),
    );
    assert!(matches!(offset_result, Err(ShellError::Cancelled)));

    let thicken_result = thicken_with_control(
        &cube,
        0.1,
        &opts,
        tol(),
        OperationControl::with_cancel(&token),
    );
    assert!(matches!(thicken_result, Err(ShellError::Cancelled)));

    assert!(cube.validate().is_ok(), "input untouched by cancelled runs");
}

/// Never-fired token + capture reporter: shell result matches the
/// plain call and the progress stream has the promised shape.
#[test]
fn never_fired_token_matches_plain_shell_and_reports_stages() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let opts = ShellOptions::default();
    let plain = shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &opts,
        tol(),
        None,
    )
    .expect("plain shell serves");

    let token = CancellationToken::new();
    let mut reporter = CaptureReporter::default();
    let controlled = shell_with_control(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &opts,
        tol(),
        None,
        OperationControl::new(&token, &mut reporter),
    )
    .expect("controlled shell serves");

    assert_eq!(plain.brep.faces.len(), controlled.brep.faces.len());
    assert_eq!(plain.brep.edges.len(), controlled.brep.edges.len());
    assert_eq!(plain.brep.vertices.len(), controlled.brep.vertices.len());
    // Whole-body fingerprint identity: the detached run must be
    // byte-identical to the plain call.
    assert_eq!(
        body_fingerprint(&plain.brep),
        body_fingerprint(&controlled.brep),
        "controlled shell must be fingerprint-identical to the plain run"
    );
    // The hollow body validates since PK-306 (void skin is its own
    // shell; χ = 2 per shell).
    assert!(
        controlled.brep.validate().is_ok(),
        "hollow shell output must validate: {:?}",
        controlled.brep.validate().err()
    );

    // All 8 shell stages, monotone, honest total, known labels.
    assert_eq!(
        reporter.events.len(),
        8,
        "shell reports its 8 stages: {:?}",
        reporter.events.iter().map(|e| e.stage).collect::<Vec<_>>()
    );
    for (i, event) in reporter.events.iter().enumerate() {
        assert_eq!(event.stage_index, i);
        assert!(event.stage.starts_with("shell."));
        assert_eq!(event.total_stages, 8);
    }
    assert_eq!(reporter.events[0].stage, "shell.classify_faces");
    assert_eq!(reporter.events[7].stage, "shell.pcurve_recompute");
}
