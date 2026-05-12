use super::*;
use crate::euler::{DiagnosticTolerance, DiagnosticToleranceSource};

fn tolerance() -> DiagnosticTolerance {
    DiagnosticTolerance {
        source: DiagnosticToleranceSource::LocalGeometryPredicate,
        nanometers: Some(10),
    }
}

#[test]
fn integration_validation_rows_carry_actions_and_tolerances() {
    let integration = classify_integration_validation(
        IntegrationValidationKind::Continuity,
        Some(false),
        tolerance(),
    );
    assert_eq!(
        integration.state,
        IntegrationValidationState::DowngradeEligible("continuity-downgrade")
    );
    assert_eq!(
        integration.action,
        IntegrationValidationAction::Downgrade("continuity-downgrade")
    );
    assert_eq!(integration.tolerance, tolerance());

    let residual = classify_integration_validation(
        IntegrationValidationKind::EdgeOnSurface,
        Some(false),
        tolerance(),
    );
    assert_eq!(
        residual.state,
        IntegrationValidationState::RepairRequired("edge-on-surface-residual-repair")
    );
    assert_eq!(
        residual.action,
        IntegrationValidationAction::Repair("edge-on-surface-residual-repair")
    );
    assert_eq!(residual.tolerance.source, tolerance().source);
    assert_eq!(residual.tolerance.nanometers, Some(10));

    let self_intersection = classify_integration_validation(
        IntegrationValidationKind::LocalSelfIntersection,
        Some(false),
        tolerance(),
    );
    assert_eq!(
        self_intersection.action,
        IntegrationValidationAction::Invalid("integration-self-intersection")
    );

    let relaxed = classify_integration_validation(
        IntegrationValidationKind::ToleranceExceeded,
        Some(false),
        tolerance(),
    );
    assert_eq!(
        relaxed.action,
        IntegrationValidationAction::RelaxedAcceptance("relaxed-tolerance-acceptance")
    );

    let semantic = classify_integration_validation(
        IntegrationValidationKind::FaceTracking,
        Some(false),
        tolerance(),
    );
    assert_eq!(
        semantic.action,
        IntegrationValidationAction::SemanticCheck("semantic-ownership-repair")
    );
}
