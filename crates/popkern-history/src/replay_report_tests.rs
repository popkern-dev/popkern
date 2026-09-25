//! Unit tests for the replay report channel: warning classification,
//! aggregation, merge semantics, and the serde contract.

use super::*;
use crate::pattern::{InstanceOutcome, PatternResult};
use crate::replay_report::evidence::FeatureEvidence;
use popkern_blend::{BlendContinuityEvidence, BlendWarning, ContinuityGates, ContinuityLevel};
use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};
use popkern_draft::DraftWarning;
use popkern_shell::ShellWarning;
use popkern_sweep::SweepWarning;

fn feature_ids(n: usize) -> Vec<FeatureId> {
    let mut map: slotmap::SlotMap<FeatureId, ()> = slotmap::SlotMap::with_key();
    (0..n).map(|_| map.insert(())).collect()
}

fn served(id: FeatureId, label: &str) -> FeatureReport {
    FeatureReport::new(id, label, "fillet", ReplayOutcome::Served)
}

fn evidence_record(requested: ContinuityLevel, achieved: ContinuityLevel) -> BlendEvidenceRecord {
    BlendEvidenceRecord {
        edges_requested: 4,
        corner_patches: 0,
        requested_continuity: requested,
        continuity: Some(BlendContinuityEvidence {
            requested,
            achieved,
            gates: ContinuityGates::commercial(),
            blend_faces: 4,
            seams: Vec::new(),
            unmeasured_seams: 0,
            sharp_seams: 0,
            tolerance_widened_seams: 0,
            certified_seams: 0,
            uncertified_seams: 0,
            certified_achieved: None,
            worst_g0_gap: 0.0,
            worst_g1_normal: 0.0,
            worst_g2_residual: 0.0,
        }),
    }
}

fn classified_blend(warning: BlendWarning) -> FeatureEvidence {
    let mut evidence = FeatureEvidence::empty();
    evidence::classify_blend_warning(&mut evidence, EdgeHandle::default(), &warning);
    evidence
}

fn counts(evidence: &FeatureEvidence) -> (usize, usize, usize) {
    (
        evidence.warnings.len(),
        evidence.fallbacks.len(),
        evidence.skipped.len(),
    )
}

#[test]
fn blend_warnings_split_into_warning_fallback_and_skip() {
    let clamped = classified_blend(BlendWarning::RadiusClamped {
        requested: 0.5,
        actual: 0.4,
    });
    assert_eq!(counts(&clamped), (1, 0, 0));
    assert_eq!(clamped.warnings[0].source, RecordSource::Blend);
    assert!(clamped.warnings[0].message.contains("RadiusClamped"));

    let skipped = classified_blend(BlendWarning::CollapsedEdgeSkipped {
        edge: EdgeHandle::default(),
        length: 1e-9,
    });
    assert_eq!(counts(&skipped), (0, 0, 1));
    assert!(skipped.skipped[0].reason.contains("1e-9"));

    let fallback = classified_blend(BlendWarning::VertexModeFallback {
        requested: popkern_blend::VertexMode::Smooth,
        actual: popkern_blend::VertexMode::Mitred,
        vertex: None,
    });
    assert_eq!(counts(&fallback), (0, 1, 0));
    assert!(
        fallback.fallbacks[0]
            .mechanism
            .ends_with("VertexModeFallback"),
        "{:?}",
        fallback.fallbacks[0]
    );
    assert!(fallback.fallbacks[0].reason.contains("Mitred"));
}

#[test]
fn shell_warnings_classify_by_variant() {
    let evidence = FeatureEvidence::shell(&[
        ShellWarning::FaceAutoRemoved {
            face: FaceHandle::default(),
            reason: "offset failed".into(),
        },
        ShellWarning::JoinStyleDowngraded {
            vertex: VertexHandle::default(),
            mechanism: "mixed corner".into(),
        },
        ShellWarning::FaceOffsetApproximated {
            face: FaceHandle::default(),
            achieved_deviation: 2e-7,
        },
        ShellWarning::ConeApexSharp {
            vertex: VertexHandle::default(),
        },
    ]);
    assert_eq!(counts(&evidence), (1, 2, 1));
    assert_eq!(evidence.skipped[0].reason, "offset failed");
    assert_eq!(
        evidence.fallbacks[0].mechanism,
        "join style Arc -> Intersection"
    );
    assert!(evidence.fallbacks[0].reason.contains("mixed corner"));
    assert!(
        evidence.fallbacks[1]
            .mechanism
            .contains("fitted approximation")
    );
    assert!(
        evidence
            .warnings
            .iter()
            .all(|w| w.source == RecordSource::Shell)
    );
}

#[test]
fn draft_and_sweep_fallbacks_are_recorded_as_fallbacks() {
    let draft = FeatureEvidence::draft(&[
        DraftWarning::AdjacentFallbackToTransition {
            face: FaceHandle::default(),
        },
        DraftWarning::EdgeCurveStale {
            edge: EdgeHandle::default(),
        },
    ]);
    assert_eq!(counts(&draft), (1, 1, 0));
    assert_eq!(draft.fallbacks[0].source, RecordSource::Draft);

    let sweep = FeatureEvidence::sweep(&[
        SweepWarning::FrenetFallback { t: 0.25 },
        SweepWarning::HelixApproximated,
        SweepWarning::SelfIntersectionAllowed { t: 0.5 },
    ]);
    assert_eq!(counts(&sweep), (1, 2, 0));
    assert_eq!(sweep.fallbacks[0].mechanism, "Frenet frame -> Bishop frame");
    assert!(sweep.warnings[0].message.contains("0.5"));
}

#[test]
fn pattern_outcomes_become_skip_and_fallback_records() {
    let partial = FeatureEvidence::pattern(&PatternResult::PartialSuccess {
        succeeded: vec![
            (0, InstanceOutcome::Replayed),
            (1, InstanceOutcome::GeometryFallback),
            (2, InstanceOutcome::Skipped),
            (
                3,
                InstanceOutcome::Failed {
                    error: "no host".into(),
                },
            ),
        ],
        failed: vec![(4, "boolean declined".into())],
        warnings: vec!["spacing clamped".into()],
    });
    assert_eq!(counts(&partial), (1, 1, 3));
    assert_eq!(partial.warnings[0].source, RecordSource::Pattern);
    assert!(partial.fallbacks[0].mechanism.starts_with("instance 1"));
    assert_eq!(partial.skipped[0].what, "instance 2");
    assert!(partial.skipped[1].reason.contains("no host"));
    assert!(partial.skipped[2].reason.contains("boolean declined"));

    let all = FeatureEvidence::pattern(&PatternResult::AllSucceeded {
        count: 3,
        warnings: Vec::new(),
    });
    assert_eq!(counts(&all), (0, 0, 0));

    let none = FeatureEvidence::pattern(&PatternResult::AllFailed {
        errors: vec![(0, "x".into()), (1, "y".into())],
    });
    assert_eq!(counts(&none), (0, 0, 2));
}

#[test]
fn boolean_fallback_events_name_the_selection_boundary() {
    let mut evidence = FeatureEvidence::empty();
    evidence.record_boolean_fallbacks(vec![popkern_bool::merge_fallback::MergeFallbackEvent {
        operation: popkern_bool::BoolOperation::Union,
        reason: popkern_bool::merge_fallback::MergeFallbackReason::CellularInvalid(
            "open shell".into(),
        ),
    }]);
    assert_eq!(counts(&evidence), (0, 1, 0));
    let record = &evidence.fallbacks[0];
    assert_eq!(record.source, RecordSource::Boolean);
    assert_eq!(
        record.mechanism,
        "Union: cellular selection -> graph-merge selection"
    );
    assert!(record.reason.contains("open shell"));
}

#[test]
fn report_is_clean_only_without_findings() {
    let ids = feature_ids(3);
    let mut report = ReplayReport {
        features: vec![served(ids[0], "A"), served(ids[1], "B")],
    };
    assert!(report.is_clean());

    report.features[0].blend = Some(evidence_record(ContinuityLevel::G1, ContinuityLevel::G1));
    assert!(report.is_clean(), "meeting the request is not a finding");

    report.features[0].blend = Some(evidence_record(ContinuityLevel::G2, ContinuityLevel::G1));
    assert!(!report.is_clean());
    assert_eq!(report.downgraded_blends().count(), 1);

    let mut report = ReplayReport {
        features: vec![FeatureReport::new(
            ids[2],
            "C",
            "shell",
            ReplayOutcome::Suppressed {
                reason: SuppressionReason::Conditional,
            },
        )],
    };
    assert!(report.is_clean(), "suppression is not a finding");
    report.features[0].outcome = ReplayOutcome::Failed {
        error: "declined".into(),
    };
    assert!(!report.is_clean());
    assert_eq!(report.failures().count(), 1);

    let mut report = ReplayReport {
        features: vec![
            served(ids[0], "A")
                .with_evidence(FeatureEvidence::sweep(&[SweepWarning::HelixApproximated])),
        ],
    };
    assert!(!report.is_clean());
    assert_eq!(report.fallbacks().count(), 1);
    report.features[0].fallbacks.clear();
    report.features[0].skipped.push(SkipRecord {
        source: RecordSource::Pattern,
        what: "instance 1".into(),
        reason: "suppressed".into(),
    });
    assert!(!report.is_clean());
    assert_eq!(report.skips().count(), 1);
}

#[test]
fn absorb_replaces_covered_features_and_orders_by_timeline() {
    let ids = feature_ids(4);
    let mut older = ReplayReport {
        features: vec![
            served(ids[0], "A"),
            served(ids[1], "B"),
            served(ids[2], "C"),
        ],
    };
    let mut newer_b = served(ids[1], "B'");
    newer_b.warnings.push(KernelWarning {
        source: RecordSource::Shell,
        message: "approximated".into(),
    });
    let newer = ReplayReport {
        features: vec![newer_b, served(ids[3], "D")],
    };
    // Timeline: C was deleted, D inserted before A.
    older.absorb(newer, &[ids[3], ids[0], ids[1]]);
    let labels: Vec<&str> = older.features.iter().map(|r| r.label.as_str()).collect();
    assert_eq!(labels, ["D", "A", "B'"]);
    assert_eq!(older.feature(ids[1]).map(|r| r.warnings.len()), Some(1));
    assert!(
        older.feature(ids[2]).is_none(),
        "records off the timeline drop"
    );
    assert_eq!(older.warnings().count(), 1);
}

#[test]
fn feature_report_round_trips_through_json() {
    let ids = feature_ids(1);
    let mut report = served(ids[0], "Fillet1").with_evidence(FeatureEvidence::shell(&[
        ShellWarning::FaceAutoRemoved {
            face: FaceHandle::default(),
            reason: "offset failed".into(),
        },
    ]));
    report.blend = Some(evidence_record(ContinuityLevel::G2, ContinuityLevel::G1));
    report.elapsed = Some(Duration::from_millis(12));
    let json = serde_json::to_string(&report).expect("serialize");
    let back: FeatureReport = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back, report);
    assert!(
        json.contains("\"G2\""),
        "continuity levels serialize by name: {json}"
    );
}
