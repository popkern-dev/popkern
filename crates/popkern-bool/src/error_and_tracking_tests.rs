//! Tests for BoolError variant coverage and boolean_ref / tracking context paths.
//!
//! Covers Gap 1 (error path coverage) and Gap 4 (boolean_ref / tracking context).

#[cfg(test)]
use popkern_core::{
    geometry::Vec3,
    history::TrackingContext,
    naming::OperationId,
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::error::{BoolError, BoolOperation, Classification, PointClassification, Stage};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Translate all vertices and geometry of a BRep by a vector offset.
fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
    for (_, vertex) in &mut brep.vertices {
        vertex.position = vertex.position + offset;
    }
    for (_, face) in &mut brep.faces {
        if let Some(ref mut sof) = face.surface {
            if let SurfaceGeometry::Plane { ref mut origin, .. } = sof.geometry {
                *origin = *origin + offset;
            }
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve {
            if let CurveGeometry::Line { ref mut origin, .. } = coe.geometry {
                *origin = *origin + offset;
            }
        }
    }
    brep
}

// ── Gap 1: BoolError variant construction and matching ─────────────────

#[test]
fn intersection_failed_error_formats_correctly() {
    let err = BoolError::IntersectionFailed("surface solver diverged".to_string());
    let msg = format!("{err}");
    assert!(
        msg.contains("intersection curve computation failed"),
        "unexpected message: {msg}"
    );
    assert!(msg.contains("surface solver diverged"));
}

#[test]
fn intersection_failed_matches_pattern() {
    let err = BoolError::IntersectionFailed("test reason".to_string());
    assert!(matches!(err, BoolError::IntersectionFailed(ref s) if s == "test reason"));
}

#[test]
fn classification_failed_error_formats_correctly() {
    let err = BoolError::ClassificationFailed {
        face: FaceHandle::default(),
        reason: "all rays inconclusive".to_string(),
    };
    let msg = format!("{err}");
    assert!(
        msg.contains("classification failed"),
        "unexpected message: {msg}"
    );
    assert!(msg.contains("all rays inconclusive"));
}

#[test]
fn classification_failed_matches_pattern() {
    let fh = FaceHandle::default();
    let err = BoolError::ClassificationFailed {
        face: fh,
        reason: "no rays".to_string(),
    };
    assert!(matches!(err, BoolError::ClassificationFailed { face, .. } if face == fh));
}

#[test]
fn split_failed_error_formats_correctly() {
    let err = BoolError::SplitFailed {
        face: FaceHandle::default(),
        reason: "degenerate loop".to_string(),
    };
    let msg = format!("{err}");
    assert!(msg.contains("face split failed"), "unexpected: {msg}");
    assert!(msg.contains("degenerate loop"));
}

#[test]
fn merge_failed_error_formats_correctly() {
    let err = BoolError::MergeFailed {
        stage: "loop_walk".to_string(),
        reason: "non-manifold edge".to_string(),
    };
    let msg = format!("{err}");
    assert!(msg.contains("merge failed"), "unexpected: {msg}");
    assert!(msg.contains("non-manifold edge"));
}

#[test]
fn validation_failed_error_formats_correctly() {
    let err = BoolError::ValidationFailed("open shell".to_string());
    let msg = format!("{err}");
    assert!(
        msg.contains("result validation failed"),
        "unexpected: {msg}"
    );
    assert!(msg.contains("open shell"));
}

#[test]
fn degenerate_input_error_formats_correctly() {
    let err = BoolError::DegenerateInput("zero-volume solid".to_string());
    let msg = format!("{err}");
    assert!(msg.contains("degenerate input"), "unexpected: {msg}");
    assert!(msg.contains("zero-volume solid"));
}

#[test]
fn coincident_faces_error_formats_correctly() {
    let err = BoolError::CoincidentFaces {
        face_a: FaceHandle::default(),
        face_b: FaceHandle::default(),
    };
    let msg = format!("{err}");
    assert!(msg.contains("coincident faces"), "unexpected: {msg}");
}

#[test]
fn no_intersection_error_formats_correctly() {
    let err = BoolError::NoIntersection;
    let msg = format!("{err}");
    assert!(msg.contains("solids do not intersect"), "unexpected: {msg}");
}

#[test]
fn bool_error_is_send_sync() {
    // Verify BoolError can be sent across threads (important for async callers).
    fn assert_send<T: Send>() {}
    fn assert_sync<T: Sync>() {}
    assert_send::<BoolError>();
    assert_sync::<BoolError>();
}

// ── Auxiliary enum coverage ────────────────────────────────────────────

#[test]
fn stage_variants_debug_format() {
    // Verify all Stage variants are representable and debuggable.
    let stages = [
        Stage::Intersection,
        Stage::Split,
        Stage::Classification,
        Stage::Selection,
        Stage::Merge,
        Stage::Validation,
    ];
    for s in &stages {
        let dbg = format!("{s:?}");
        assert!(!dbg.is_empty());
    }
}

#[test]
fn classification_variants_are_distinct() {
    assert_ne!(Classification::Inside, Classification::Outside);
    assert_ne!(
        Classification::OnBoundaryCooriented,
        Classification::OnBoundaryCrossing
    );
}

#[test]
fn point_classification_variants_are_distinct() {
    assert_ne!(PointClassification::Inside, PointClassification::Outside);
    assert_ne!(PointClassification::Inside, PointClassification::OnBoundary);
}

#[test]
fn bool_operation_variants_are_distinct() {
    assert_ne!(BoolOperation::Union, BoolOperation::Intersection);
    assert_ne!(BoolOperation::Union, BoolOperation::Subtraction);
    assert_ne!(BoolOperation::Intersection, BoolOperation::Subtraction);
}

// ── Gap 1 continued: trigger error paths through the public API ────────

#[test]
fn empty_brep_a_returns_degenerate_input() {
    let a = BRep::new(tol());
    let b = BRep::unit_cube(tol()).expect("unit cube");
    let result = crate::boolean(BoolOperation::Union, a, b, tol(), None);
    assert!(result.is_err());
    assert!(
        matches!(result.unwrap_err(), BoolError::DegenerateInput(_)),
        "expected DegenerateInput for empty A"
    );
}

#[test]
fn empty_brep_b_returns_degenerate_input() {
    let a = BRep::unit_cube(tol()).expect("unit cube");
    let b = BRep::new(tol());
    let result = crate::boolean(BoolOperation::Union, a, b, tol(), None);
    assert!(result.is_err());
    assert!(
        matches!(result.unwrap_err(), BoolError::DegenerateInput(_)),
        "expected DegenerateInput for empty B"
    );
}

#[test]
fn both_empty_breps_returns_degenerate_input() {
    let a = BRep::new(tol());
    let b = BRep::new(tol());
    let result = crate::boolean(BoolOperation::Intersection, a, b, tol(), None);
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), BoolError::DegenerateInput(_)));
}

// ── Gap 4: boolean_ref and tracking context paths ──────────────────────

#[test]
fn boolean_ref_union_matches_boolean() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );

    let result_ref = crate::boolean_ref(BoolOperation::Union, &a, &b, tol(), None);
    assert!(
        result_ref.is_ok(),
        "boolean_ref union failed: {:?}",
        result_ref.err()
    );
    let brep_ref = result_ref.unwrap();

    // Also run the consuming version to compare face count.
    let result_owned = crate::boolean(BoolOperation::Union, a, b, tol(), None);
    assert!(result_owned.is_ok());
    let brep_owned = result_owned.unwrap();

    assert_eq!(
        brep_ref.faces.len(),
        brep_owned.faces.len(),
        "boolean_ref and boolean should produce same face count"
    );
}

#[test]
fn boolean_ref_subtraction_works() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::boolean_ref(BoolOperation::Subtraction, &a, &b, tol(), None);
    assert!(
        result.is_ok(),
        "boolean_ref subtraction failed: {:?}",
        result.err()
    );
}

#[test]
fn boolean_ref_intersection_works() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::boolean_ref(BoolOperation::Intersection, &a, &b, tol(), None);
    assert!(
        result.is_ok(),
        "boolean_ref intersection failed: {:?}",
        result.err()
    );
}

#[test]
fn boolean_ref_preserves_inputs() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let a_faces_before = a.faces.len();
    let b_faces_before = b.faces.len();

    let _ = crate::boolean_ref(BoolOperation::Union, &a, &b, tol(), None);

    // Inputs must be unchanged after boolean_ref.
    assert_eq!(a.faces.len(), a_faces_before);
    assert_eq!(b.faces.len(), b_faces_before);
}

#[test]
fn boolean_ref_error_preserves_inputs() {
    let a = BRep::new(tol());
    let b = BRep::unit_cube(tol()).expect("cube b");
    let a_verts = a.vertices.len();

    let result = crate::boolean_ref(BoolOperation::Union, &a, &b, tol(), None);
    assert!(result.is_err());

    // Input a is still intact.
    assert_eq!(a.vertices.len(), a_verts);
}

#[test]
fn boolean_with_tracking_none_succeeds() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::boolean(BoolOperation::Union, a, b, tol(), None);
    assert!(result.is_ok());
}

#[test]
fn boolean_with_tracking_context_assigns_names() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let mut ctx = TrackingContext::new(OperationId(42));

    let result = crate::boolean(BoolOperation::Union, a, b, tol(), Some(&mut ctx));
    assert!(
        result.is_ok(),
        "boolean with tracking failed: {:?}",
        result.err()
    );
    let brep = result.unwrap();

    // Every face in the result should have a name assigned.
    let named_count = brep.faces.iter().filter(|(_, f)| f.name.is_some()).count();
    assert_eq!(
        named_count,
        brep.faces.len(),
        "all result faces should be named when tracking is active"
    );
}

#[test]
fn boolean_ref_with_tracking_context_assigns_names() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let mut ctx = TrackingContext::new(OperationId(99));

    let result = crate::boolean_ref(BoolOperation::Union, &a, &b, tol(), Some(&mut ctx));
    assert!(result.is_ok());
    let brep = result.unwrap();

    let named_count = brep.faces.iter().filter(|(_, f)| f.name.is_some()).count();
    assert_eq!(named_count, brep.faces.len());
}

#[test]
fn tracking_context_records_history() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let mut ctx = TrackingContext::new(OperationId(7));

    let result = crate::boolean(BoolOperation::Union, a, b, tol(), Some(&mut ctx));
    assert!(result.is_ok());

    let history = ctx.finish();
    // A union of two overlapping cubes should produce both modified and
    // generated faces. At minimum, some faces are modified from the inputs.
    let total_records = history.modified.len() + history.generated.len();
    assert!(
        total_records > 0,
        "tracking should record at least one modified or generated entity"
    );
}

#[test]
fn tracking_disabled_does_not_record() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let mut ctx = TrackingContext::new(OperationId(8));
    ctx.enabled = false;

    let result = crate::boolean(BoolOperation::Union, a, b, tol(), Some(&mut ctx));
    assert!(result.is_ok());

    let history = ctx.finish();
    assert_eq!(history.modified.len(), 0);
    assert_eq!(history.generated.len(), 0);
}

#[test]
fn boolean_with_tracking_disjoint_cubes_still_names_faces() {
    // Disjoint cubes follow the no-intersection shortcut, so tracking
    // still runs on the result.
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let mut ctx = TrackingContext::new(OperationId(10));

    let result = crate::boolean(BoolOperation::Union, a, b, tol(), Some(&mut ctx));
    assert!(result.is_ok());
    let brep = result.unwrap();

    // All faces should be named.
    let named_count = brep.faces.iter().filter(|(_, f)| f.name.is_some()).count();
    assert_eq!(named_count, brep.faces.len());
}

// ── Convenience wrapper tests ──────────────────────────────────────────

#[test]
fn ref_convenience_wrappers_work() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );
    assert!(crate::union_ref(&a, &b, tol(), None).is_ok());
    assert!(crate::intersection_ref(&a, &b, tol(), None).is_ok());
    assert!(crate::subtraction_ref(&a, &b, tol(), None).is_ok());
}

#[test]
fn owned_convenience_wrappers_work() {
    let mk = || {
        (
            BRep::unit_cube(tol()).expect("a"),
            translate(BRep::unit_cube(tol()).expect("b"), Vec3::new(0.5, 0.3, 0.2)),
        )
    };
    let (a, b) = mk();
    assert!(crate::union(a, b, tol(), None).is_ok());
    let (a, b) = mk();
    assert!(crate::intersection(a, b, tol(), None).is_ok());
    let (a, b) = mk();
    assert!(crate::subtraction(a, b, tol(), None).is_ok());
}

// ── §13.4 Remediation hint coverage ─────────────────────────

#[test]
fn every_bool_error_variant_has_a_remediation_hint() {
    use crate::error::RemediationAction;
    use popkern_core::error::CoreError;

    // Construct a representative instance of every variant and
    // verify that `remediation_hint()` returns a non-empty
    // rationale string. Pins §13.4: every BoolError must carry
    // actionable remediation data.
    let variants = vec![
        BoolError::NoIntersection,
        BoolError::IntersectionFailed("solver failed".into()),
        BoolError::SplitFailed {
            face: FaceHandle::default(),
            reason: "test".into(),
        },
        BoolError::ClassificationFailed {
            face: FaceHandle::default(),
            reason: "test".into(),
        },
        BoolError::MergeFailed {
            stage: "Test".into(),
            reason: "test".into(),
        },
        BoolError::ValidationFailed("test".into()),
        BoolError::DegenerateInput("zero volume".into()),
        BoolError::CoincidentFaces {
            face_a: FaceHandle::default(),
            face_b: FaceHandle::default(),
        },
        BoolError::CoreError(CoreError::InvalidHandle("test".into())),
        BoolError::Cancelled,
        BoolError::OpenShellNotSupported { a: true, b: false },
        BoolError::InterferenceFailed("test".into()),
        BoolError::SystemError("OOM".into()),
    ];

    for err in &variants {
        let hint = err.remediation_hint();
        assert!(
            !hint.rationale.is_empty(),
            "variant {err:?} returned an empty rationale"
        );
        // Every action must be a known variant.
        match hint.action {
            RemediationAction::IncreaseTolerance
            | RemediationAction::DecreaseTolerance
            | RemediationAction::HealInput
            | RemediationAction::UseCoincidenceMode
            | RemediationAction::SimplifyGeometry
            | RemediationAction::RetryWithPerturbation
            | RemediationAction::NoAction
            | RemediationAction::Precheck => {}
        }
    }
}

#[test]
fn remediation_hint_for_noin_intersection_is_no_action() {
    // Specific check: NoIntersection shouldn't suggest a retry or
    // tolerance bump — it's a semantic signal, not a failure.
    use crate::error::RemediationAction;
    let hint = BoolError::NoIntersection.remediation_hint();
    assert_eq!(hint.action, RemediationAction::NoAction);
}

#[test]
fn remediation_hint_for_cancelled_is_no_action() {
    use crate::error::RemediationAction;
    let hint = BoolError::Cancelled.remediation_hint();
    assert_eq!(hint.action, RemediationAction::NoAction);
}

#[test]
fn remediation_hint_for_coincident_faces_suggests_coincidence_mode() {
    use crate::error::RemediationAction;
    let hint = BoolError::CoincidentFaces {
        face_a: FaceHandle::default(),
        face_b: FaceHandle::default(),
    }
    .remediation_hint();
    assert_eq!(hint.action, RemediationAction::UseCoincidenceMode);
}

#[test]
fn remediation_hint_for_open_shell_suggests_precheck() {
    use crate::error::RemediationAction;
    let hint = BoolError::OpenShellNotSupported { a: true, b: false }.remediation_hint();
    assert_eq!(hint.action, RemediationAction::Precheck);
    assert!(
        hint.rationale.contains("boolean_with_sheets"),
        "rationale should point callers at the sheet-aware entry \
         point, got: {}",
        hint.rationale
    );
}

#[test]
fn fill_missing_remediations_populates_error_entries() {
    use crate::error::{BoolDiagnostics, DiagnosticEntry, Severity, Stage};

    let mut diag = BoolDiagnostics {
        error: Some(BoolError::IntersectionFailed("test".into())),
        failed_stage: Some(Stage::Intersection),
        entries: vec![
            DiagnosticEntry {
                severity: Severity::Error,
                phase: Stage::Intersection,
                message: "test".into(),
                entities: Vec::new(),
                location: None,
                measured_value: None,
                remediation: None,
            },
            DiagnosticEntry {
                severity: Severity::Warning,
                phase: Stage::Intersection,
                message: "warning".into(),
                entities: Vec::new(),
                location: None,
                measured_value: None,
                remediation: None,
            },
        ],
        ..Default::default()
    };

    diag.fill_missing_remediations();

    // Error entry gets filled.
    assert!(diag.entries[0].remediation.is_some());
    // Warning entry is left alone (fill is only for Error/Critical).
    assert!(diag.entries[1].remediation.is_none());
}

#[test]
fn from_error_builder_populates_entry() {
    use crate::error::{BoolDiagnostics, Severity, Stage};

    let diag = BoolDiagnostics::from_error(
        BoolError::ValidationFailed("gap of 0.01".into()),
        Stage::Validation,
    );
    assert_eq!(diag.entries.len(), 1);
    assert_eq!(diag.entries[0].severity, Severity::Error);
    assert_eq!(diag.entries[0].phase, Stage::Validation);
    assert!(diag.entries[0].remediation.is_some());
    assert_eq!(diag.failed_stage, Some(Stage::Validation));
    assert!(diag.error.is_some());
}

#[test]
fn context_aware_hint_fills_suggested_value_from_tolerance() {
    use crate::error::{RemediationAction, RemediationContext};
    use popkern_core::tolerance::Tolerance;

    let tol = Tolerance::default();
    let ctx = RemediationContext {
        tolerance: Some(&tol),
        perturbation_attempted: false,
    };

    // SplitFailed is a tolerance-bump case; the suggested value
    // should be the current linear tolerance × 10.
    let hint = BoolError::SplitFailed {
        face: FaceHandle::default(),
        reason: "test".into(),
    }
    .remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::IncreaseTolerance);
    assert_eq!(hint.parameter.as_deref(), Some("linear"));
    let expected = tol.linear() * 10.0;
    assert_eq!(hint.value, Some(expected));
}

#[test]
fn context_aware_hint_without_tolerance_has_no_value() {
    use crate::error::{RemediationAction, RemediationContext};

    // No tolerance in context — the hint has the action and
    // parameter but value stays None.
    let ctx = RemediationContext::default();
    let hint = BoolError::SplitFailed {
        face: FaceHandle::default(),
        reason: "test".into(),
    }
    .remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::IncreaseTolerance);
    assert_eq!(hint.parameter.as_deref(), Some("linear"));
    assert!(hint.value.is_none());
}

#[test]
fn context_aware_hint_avoids_perturbation_loop() {
    use crate::error::{RemediationAction, RemediationContext};

    // After perturbation was attempted, IntersectionFailed should
    // recommend IncreaseTolerance instead of RetryWithPerturbation.
    let ctx = RemediationContext {
        tolerance: None,
        perturbation_attempted: true,
    };
    let hint =
        BoolError::IntersectionFailed("solver stuck".into()).remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::IncreaseTolerance);
    assert!(
        !hint.rationale.contains("boolean_with_perturbation"),
        "rationale should not loop back to perturbation after it's been tried"
    );
}

#[test]
fn context_aware_hint_suggests_perturbation_when_untried() {
    use crate::error::{RemediationAction, RemediationContext};

    let ctx = RemediationContext::default();
    let hint =
        BoolError::IntersectionFailed("solver stuck".into()).remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::RetryWithPerturbation);
}

#[test]
fn context_aware_hint_coincidence_tolerance_scales_correctly() {
    use crate::error::{RemediationAction, RemediationContext};
    use popkern_core::tolerance::Tolerance;

    let tol = Tolerance::default();
    let ctx = RemediationContext {
        tolerance: Some(&tol),
        perturbation_attempted: false,
    };

    // ClassificationFailed suggests coincidence × 2. In the helper
    // the coincidence parameter maps to linear × 10 × factor.
    let hint = BoolError::ClassificationFailed {
        face: FaceHandle::default(),
        reason: "test".into(),
    }
    .remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::IncreaseTolerance);
    assert_eq!(hint.parameter.as_deref(), Some("coincidence"));
    assert!(hint.value.is_some());
    let expected = tol.linear() * 10.0 * 2.0;
    assert_eq!(hint.value, Some(expected));
}

#[test]
fn boolean_with_remediation_returns_ok_for_valid_inputs() {
    // Smoke test the new entry point: a successful boolean on
    // overlapping cubes should return Ok with empty diagnostics
    // (no Remediation hint — there's no error).
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let result = crate::boolean_with_remediation(BoolOperation::Union, &a, &b, tol());
    assert!(result.is_ok(), "union of overlapping cubes should succeed");
}
