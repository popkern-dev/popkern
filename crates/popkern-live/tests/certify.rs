//! PK-561: deferred continuity certification through popkern live
//! session — a deferred run serves on sampled evidence with no
//! certified claim, `certify` slices finish with the upgraded report,
//! and the result equals an inline run's evidence.

use std::path::Path;
use std::time::Duration;

use popkern_live::{ContinuityCertification, HostOptions, Outcome, Session};
use serde_json::Value;

fn boss_plate() -> String {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../popkern-rhai/examples/boss_plate.rhai");
    std::fs::read_to_string(path).expect("boss plate example")
}

fn session(policy: ContinuityCertification) -> Session {
    Session::new(HostOptions {
        continuity_certification: policy,
        ..HostOptions::default()
    })
}

/// The report's blend continuity records with the deferral bookkeeping
/// stripped, for comparing a certified-after record with an inline one.
fn continuity(report: &Value) -> Vec<Value> {
    report["features"]
        .as_array()
        .expect("features")
        .iter()
        .filter_map(|row| row["blend"]["continuity"].as_object().cloned())
        .map(|mut record| {
            record.remove("deferred_seams");
            if let Some(Value::Array(seams)) = record.get_mut("seams") {
                for seam in seams {
                    if let Value::Object(seam) = seam {
                        seam.remove("certification_deferred");
                    }
                }
            }
            Value::Object(record)
        })
        .collect()
}

fn deferred_seams(report: &Value) -> u64 {
    report["features"]
        .as_array()
        .expect("features")
        .iter()
        .filter_map(|row| row["blend"]["continuity"]["deferred_seams"].as_u64())
        .sum()
}

#[test]
fn a_deferred_run_serves_sampled_and_certify_slices_reproduce_the_inline_report() {
    let mut inline = session(ContinuityCertification::Inline);
    let inline_envelope = inline.run(&boss_plate(), "boss_plate.rhai");
    assert!(inline_envelope.outcome.is_served());
    let inline_report = inline_envelope.report.expect("report");
    assert_eq!(deferred_seams(&inline_report), 0);

    let mut deferred = session(ContinuityCertification::Deferred);
    let envelope = deferred.run(&boss_plate(), "boss_plate.rhai");
    assert!(envelope.outcome.is_served(), "{:?}", envelope.outcome);
    let report = envelope.report.expect("report");
    let total = deferred_seams(&report);
    assert!(total > 0, "the fillet's gated seams are deferred");
    let rounds = &continuity(&report)[0];
    assert_eq!(
        rounds["certified_achieved"],
        Value::Null,
        "no certified claim yet"
    );

    let mut remaining = total as usize;
    let mut last = None;
    while remaining > 0 {
        // A zero slice certifies exactly one seam per call.
        let step = deferred.certify(Duration::ZERO);
        assert!(step.failure.is_none(), "{:?}", step.failure);
        assert_eq!(step.certified, 1);
        assert_eq!(step.remaining, remaining - 1);
        assert_eq!(step.report.is_some(), step.remaining == 0);
        remaining = step.remaining;
        last = step.report;
    }
    let certified = last.expect("final report");
    assert_eq!(deferred_seams(&certified), 0);
    assert_eq!(continuity(&certified), continuity(&inline_report));

    // Nothing left: an empty step that still carries the report.
    let idle = deferred.certify(Duration::from_millis(250));
    assert_eq!((idle.certified, idle.remaining), (0, 0));
    assert!(idle.report.is_some());
}

#[test]
fn a_parameter_edit_replaces_the_deferred_seams() {
    let mut deferred = session(ContinuityCertification::Deferred);
    let envelope = deferred.run(&boss_plate(), "boss_plate.rhai");
    let before = deferred_seams(&envelope.report.expect("report"));
    let step = deferred.certify(Duration::ZERO);
    assert_eq!(step.remaining as u64, before - 1);
    let edited = deferred.set_param("extrude-1.distance", 0.5);
    assert!(edited.outcome.is_served(), "{:?}", edited.outcome);
    assert_eq!(
        deferred_seams(&edited.report.expect("report")),
        before,
        "the replay's fresh evidence defers every gated seam again"
    );
}

#[test]
fn certify_without_a_model_is_a_gate_failure() {
    let mut empty = session(ContinuityCertification::Deferred);
    let step = empty.certify(Duration::from_millis(10));
    assert!(matches!(
        step.failure,
        Some(Outcome::GateFailure { ref step, .. }) if step == "model"
    ));
    assert_eq!((step.certified, step.remaining), (0, 0));
}
