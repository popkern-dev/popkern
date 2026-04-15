//! §13.4 Remediation hint coverage tests. Split out of
//! `error_and_tracking_tests.rs` to keep that file under the
//! 500-line structural cap.

#[cfg(test)]
use popkern_core::{
    geometry::Vec3,
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::error::{BoolError, BoolOperation};

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
        if let Some(ref mut sof) = face.surface
            && let SurfaceGeometry::Plane { ref mut origin, .. } = sof.geometry
        {
            *origin = *origin + offset;
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve
            && let CurveGeometry::Line { ref mut origin, .. } = coe.geometry
        {
            *origin = *origin + offset;
        }
    }
    brep
}

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
