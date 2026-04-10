use super::*;
use crate::topology::{EdgeHandle, FaceHandle, VertexHandle};

// ---------------------------------------------------------------------------
// EntityHandle
// ---------------------------------------------------------------------------

#[test]
fn entity_handle_from_impls() {
    let vh = VertexHandle::default();
    let eh_raw = EdgeHandle::default();
    let fh = FaceHandle::default();

    assert!(matches!(EntityHandle::from(vh), EntityHandle::Vertex(_)));
    assert!(matches!(EntityHandle::from(eh_raw), EntityHandle::Edge(_)));
    assert!(matches!(EntityHandle::from(fh), EntityHandle::Face(_)));
}

#[test]
fn entity_handle_display() {
    let vh = VertexHandle::default();
    let eh_raw = EdgeHandle::default();
    let fh = FaceHandle::default();

    assert!(format!("{}", EntityHandle::Vertex(vh)).starts_with("Vertex("));
    assert!(format!("{}", EntityHandle::Edge(eh_raw)).starts_with("Edge("));
    assert!(format!("{}", EntityHandle::Face(fh)).starts_with("Face("));
}

#[test]
fn entity_handle_clone_copy_eq_hash() {
    use std::collections::HashSet;
    let vh = VertexHandle::default();
    let eh = EntityHandle::Vertex(vh);
    let copied = eh;
    assert_eq!(eh, copied);

    let mut set = HashSet::new();
    set.insert(eh);
    assert!(set.contains(&EntityHandle::Vertex(vh)));

    // Different variants are not equal
    assert_ne!(eh, EntityHandle::Edge(EdgeHandle::default()));
}

// ---------------------------------------------------------------------------
// Containment
// ---------------------------------------------------------------------------

#[test]
fn containment_variants_and_equality() {
    let c = Containment::Inside;
    let c2 = c; // Copy
    assert_eq!(c, c2);
    assert_eq!(Containment::Outside, Containment::Outside);
    assert_eq!(Containment::OnBoundary, Containment::OnBoundary);
    assert_ne!(Containment::Inside, Containment::Outside);
    assert_ne!(Containment::Inside, Containment::OnBoundary);
    assert_ne!(Containment::Outside, Containment::OnBoundary);
}

// ---------------------------------------------------------------------------
// DiagnosticSeverity
// ---------------------------------------------------------------------------

#[test]
fn severity_variants_and_copy() {
    let s = DiagnosticSeverity::Warning;
    let s2 = s; // Copy
    assert_eq!(s, s2);
    assert_ne!(DiagnosticSeverity::Error, DiagnosticSeverity::Warning);
    assert_ne!(DiagnosticSeverity::Warning, DiagnosticSeverity::Info);
}

// ---------------------------------------------------------------------------
// PrecisionDiagnosticKind
// ---------------------------------------------------------------------------

#[test]
fn diagnostic_kind_clone_and_fields() {
    let kind = PrecisionDiagnosticKind::ExceedsMaxBound;
    let _cloned = kind.clone();

    let kind = PrecisionDiagnosticKind::NumericalFailure {
        reason: "division by zero".into(),
    };
    if let PrecisionDiagnosticKind::NumericalFailure { reason } = &kind {
        assert_eq!(reason, "division by zero");
    }

    let kind = PrecisionDiagnosticKind::PrecisionBudgetWarning {
        estimated_remaining: 42,
    };
    if let PrecisionDiagnosticKind::PrecisionBudgetWarning {
        estimated_remaining,
    } = &kind
    {
        assert_eq!(*estimated_remaining, 42);
    }
}

// ---------------------------------------------------------------------------
// Remediation
// ---------------------------------------------------------------------------

#[test]
fn remediation_display() {
    let vh = VertexHandle::default();
    let eh = EdgeHandle::default();

    let cases: Vec<(Remediation, &str)> = vec![
        (
            Remediation::RecoverVertexPrecision(vh),
            "recover vertex precision",
        ),
        (
            Remediation::RecoverEdgePrecision(eh),
            "recover edge precision",
        ),
        (
            Remediation::RecoverEdgeConsistency(eh),
            "recover edge consistency",
        ),
        (
            Remediation::Manual("re-mesh the face".into()),
            "re-mesh the face",
        ),
        (Remediation::None, "no remediation"),
    ];

    for (r, expected) in cases {
        assert!(
            format!("{r}").contains(expected),
            "Remediation display for {r:?} should contain '{expected}'"
        );
    }
}

// ---------------------------------------------------------------------------
// PrecisionDiagnostic
// ---------------------------------------------------------------------------

#[test]
fn precision_diagnostic_construction_and_clone() {
    let vh = VertexHandle::default();
    let diag = PrecisionDiagnostic {
        severity: DiagnosticSeverity::Error,
        kind: PrecisionDiagnosticKind::ExceedsMaxBound,
        entities: vec![EntityHandle::Vertex(vh)],
        measured: 1e-4,
        threshold: 1e-6,
        remediation: Remediation::RecoverVertexPrecision(vh),
    };
    assert_eq!(diag.severity, DiagnosticSeverity::Error);
    assert_eq!(diag.entities.len(), 1);
    assert!(diag.measured > diag.threshold);

    let _cloned = diag.clone();
}

// ---------------------------------------------------------------------------
// PrecisionHealth
// ---------------------------------------------------------------------------

#[test]
fn health_from_budget_ratio() {
    // Good: ratio >= 0.5
    assert_eq!(
        PrecisionHealth::from_budget_ratio(0.8),
        PrecisionHealth::Good
    );
    assert_eq!(
        PrecisionHealth::from_budget_ratio(0.5),
        PrecisionHealth::Good
    );
    // Degraded: 0.1 <= ratio < 0.5
    assert_eq!(
        PrecisionHealth::from_budget_ratio(0.3),
        PrecisionHealth::Degraded
    );
    assert_eq!(
        PrecisionHealth::from_budget_ratio(0.1),
        PrecisionHealth::Degraded
    );
    // Critical: ratio < 0.1
    assert_eq!(
        PrecisionHealth::from_budget_ratio(0.09),
        PrecisionHealth::Critical
    );
    assert_eq!(
        PrecisionHealth::from_budget_ratio(0.0),
        PrecisionHealth::Critical
    );
}

#[test]
fn health_clone_copy() {
    let h = PrecisionHealth::Good;
    let h2 = h; // Copy
    assert_eq!(h, h2);
}

// ---------------------------------------------------------------------------
// PrecisionSummary
// ---------------------------------------------------------------------------

#[test]
fn precision_summary_construction() {
    let summary = PrecisionSummary {
        health_score: 0.95,
        health: PrecisionHealth::Good,
        worst_entity: None,
        degraded_count: 0,
        critical_count: 0,
        estimated_remaining: Some(100),
    };
    assert_eq!(summary.health, PrecisionHealth::Good);
    assert!(summary.worst_entity.is_none());

    let vh = VertexHandle::default();
    let summary2 = PrecisionSummary {
        health_score: 0.3,
        health: PrecisionHealth::Degraded,
        worst_entity: Some((EntityHandle::Vertex(vh), 1e-4)),
        degraded_count: 3,
        critical_count: 0,
        estimated_remaining: Some(10),
    };
    assert_eq!(summary2.degraded_count, 3);
    assert!(summary2.worst_entity.is_some());
}

// ---------------------------------------------------------------------------
// OperationPrecisionReport
// ---------------------------------------------------------------------------

#[test]
fn operation_report_empty_and_with_data() {
    let report = OperationPrecisionReport {
        entity_deltas: vec![],
        max_imprecision: 0.0,
        exceeded_bound: vec![],
        diagnostics: vec![],
    };
    assert_eq!(report.max_imprecision, 0.0);
    assert!(report.exceeded_bound.is_empty());

    let vh = VertexHandle::default();
    let entity = EntityHandle::Vertex(vh);
    let report2 = OperationPrecisionReport {
        entity_deltas: vec![(entity, 1e-8, 1e-5)],
        max_imprecision: 1e-5,
        exceeded_bound: vec![entity],
        diagnostics: vec![PrecisionDiagnostic {
            severity: DiagnosticSeverity::Error,
            kind: PrecisionDiagnosticKind::ExceedsMaxBound,
            entities: vec![entity],
            measured: 1e-5,
            threshold: 1e-6,
            remediation: Remediation::RecoverVertexPrecision(vh),
        }],
    };
    assert_eq!(report2.entity_deltas.len(), 1);
    assert_eq!(report2.diagnostics.len(), 1);
}

// ---------------------------------------------------------------------------
// NearCoincidenceConfig
// ---------------------------------------------------------------------------

#[test]
fn near_coincidence_default_and_copy() {
    let cfg = NearCoincidenceConfig::default();
    assert!((cfg.range_multiplier - 3.0).abs() < f64::EPSILON);
    assert!((cfg.max_multiplier - 10.0).abs() < f64::EPSILON);
    let cfg2 = cfg; // Copy
    assert!((cfg.range_multiplier - cfg2.range_multiplier).abs() < f64::EPSILON);
}

// ---------------------------------------------------------------------------
// DegradationTracker
// ---------------------------------------------------------------------------

fn make_report(max_imprecision: f64) -> OperationPrecisionReport {
    OperationPrecisionReport {
        entity_deltas: vec![],
        max_imprecision,
        exceeded_bound: vec![],
        diagnostics: vec![],
    }
}

#[test]
fn degradation_tracker_empty_and_single() {
    let tracker = DegradationTracker::new(10);
    assert_eq!(tracker.degradation_rate(), 0.0);

    let mut tracker = DegradationTracker::new(10);
    tracker.record("op1", &make_report(1e-6));
    assert_eq!(tracker.degradation_rate(), 0.0); // <2 entries
}

#[test]
fn degradation_tracker_rate() {
    let mut tracker = DegradationTracker::new(10);
    for (i, name) in ["op1", "op2", "op3"].iter().enumerate() {
        tracker.record(name, &make_report((i + 1) as f64 * 1e-6));
    }
    // first=1e-6, last=3e-6, steps=2 -> rate = 1e-6
    assert!((tracker.degradation_rate() - 1e-6).abs() < 1e-15);
}

#[test]
fn degradation_tracker_window_eviction() {
    let mut tracker = DegradationTracker::new(3);
    for i in 0..5 {
        tracker.record("op", &make_report((i + 1) as f64 * 1e-6));
    }
    // Window: 3e-6, 4e-6, 5e-6 -> rate = (5e-6 - 3e-6)/2 = 1e-6
    assert!((tracker.degradation_rate() - 1e-6).abs() < 1e-15);
}

#[test]
fn degradation_tracker_estimated_remaining() {
    let mut tracker = DegradationTracker::new(10);
    for i in 0..3 {
        tracker.record("op", &make_report((i + 1) as f64 * 1e-6));
    }
    // rate=1e-6, remaining=(1e-3 - 3e-6)/1e-6 = 997
    assert_eq!(
        tracker.estimated_operations_remaining(3e-6, 1e-3),
        Some(997)
    );
}

#[test]
fn degradation_tracker_remaining_edge_cases() {
    // Non-degrading -> None
    let mut tracker = DegradationTracker::new(10);
    for _ in 0..3 {
        tracker.record("op", &make_report(1e-6));
    }
    assert!(tracker.estimated_operations_remaining(1e-6, 1e-3).is_none());

    // Already exceeded -> Some(0)
    let mut tracker = DegradationTracker::new(10);
    for i in 0..3 {
        tracker.record("op", &make_report((i + 1) as f64 * 1e-3));
    }
    assert_eq!(tracker.estimated_operations_remaining(5e-3, 1e-3), Some(0));
}

// ---------------------------------------------------------------------------
// SolverConfig
// ---------------------------------------------------------------------------

#[test]
fn solver_config_default_and_clone() {
    let cfg = SolverConfig::default();
    assert_eq!(cfg.max_iterations, 50);
    assert_eq!(cfg.max_fallback_iterations, 100);
    assert!(cfg.allow_escalation);
    assert!(cfg.allow_fallback);
    let _cloned = cfg.clone();
}

// ---------------------------------------------------------------------------
// SolverDiagnostics and ConvergenceOutcome
// ---------------------------------------------------------------------------

#[test]
fn solver_diagnostics_converged() {
    let diag = SolverDiagnostics {
        iterations: 10,
        final_residual: 1e-12,
        achieved_precision: 1e-10,
        convergence: ConvergenceOutcome::Converged,
    };
    assert_eq!(diag.iterations, 10);
    assert!(matches!(diag.convergence, ConvergenceOutcome::Converged));
    let _cloned = diag.clone();
}

#[test]
fn solver_diagnostics_relaxed() {
    let diag = SolverDiagnostics {
        iterations: 50,
        final_residual: 1e-8,
        achieved_precision: 1e-8,
        convergence: ConvergenceOutcome::ConvergedRelaxed {
            requested: 1e-10,
            achieved: 1e-8,
        },
    };
    if let ConvergenceOutcome::ConvergedRelaxed {
        requested,
        achieved,
    } = &diag.convergence
    {
        assert!(*achieved > *requested);
    } else {
        panic!("wrong variant");
    }
}

#[test]
fn solver_diagnostics_fallback() {
    let diag = SolverDiagnostics {
        iterations: 30,
        final_residual: 1e-9,
        achieved_precision: 1e-9,
        convergence: ConvergenceOutcome::FallbackMethod {
            method: "bisection",
        },
    };
    if let ConvergenceOutcome::FallbackMethod { method } = &diag.convergence {
        assert_eq!(*method, "bisection");
    } else {
        panic!("wrong variant");
    }
}

#[test]
fn solver_diagnostics_precision_escalated() {
    let diag = SolverDiagnostics {
        iterations: 40,
        final_residual: 1e-11,
        achieved_precision: 1e-11,
        convergence: ConvergenceOutcome::PrecisionEscalated,
    };
    assert!(matches!(
        diag.convergence,
        ConvergenceOutcome::PrecisionEscalated
    ));
}

#[test]
fn solver_diagnostics_failure_reasons() {
    let reasons = [
        FailureReason::MaxIterations,
        FailureReason::Divergence,
        FailureReason::Oscillation,
        FailureReason::IllConditioned {
            condition_estimate: 1e12,
        },
        FailureReason::Singular,
    ];
    for reason in reasons {
        let diag = SolverDiagnostics {
            iterations: 100,
            final_residual: 1.0,
            achieved_precision: 1.0,
            convergence: ConvergenceOutcome::Failed {
                reason: reason.clone(),
            },
        };
        assert!(matches!(
            diag.convergence,
            ConvergenceOutcome::Failed { .. }
        ));
    }
}
