//! Integration validation classifiers for corner/apex diagnostics.
#![allow(missing_docs)]

use super::types::{IntegrationValidationKind, IntegrationValidationState};
use super::validation_types::IntegrationValidationAction;

pub(super) fn integration_validation_state_action(
    kind: IntegrationValidationKind,
    passes: Option<bool>,
) -> (IntegrationValidationState, IntegrationValidationAction) {
    match passes {
        Some(true) => (
            IntegrationValidationState::Ready,
            IntegrationValidationAction::None,
        ),
        None => (
            IntegrationValidationState::UnknownInputFact,
            IntegrationValidationAction::UnknownInputFact,
        ),
        Some(false) => failed_validation_state_action(kind),
    }
}

fn failed_validation_state_action(
    kind: IntegrationValidationKind,
) -> (IntegrationValidationState, IntegrationValidationAction) {
    match kind {
        IntegrationValidationKind::Watertightness => repair("watertightness-repair"),
        IntegrationValidationKind::EdgeOnSurface => repair("edge-on-surface-residual-repair"),
        IntegrationValidationKind::VertexClosure => repair("vertex-closure-residual-repair"),
        IntegrationValidationKind::LocalSelfIntersection
        | IntegrationValidationKind::AdjacentSurfaceIntersection => (
            IntegrationValidationState::Invalid("integration-self-intersection"),
            IntegrationValidationAction::Invalid("integration-self-intersection"),
        ),
        IntegrationValidationKind::ToleranceExceeded => (
            IntegrationValidationState::DowngradeEligible("relaxed-tolerance-acceptance"),
            IntegrationValidationAction::RelaxedAcceptance("relaxed-tolerance-acceptance"),
        ),
        IntegrationValidationKind::Continuity => (
            IntegrationValidationState::DowngradeEligible("continuity-downgrade"),
            IntegrationValidationAction::Downgrade("continuity-downgrade"),
        ),
        IntegrationValidationKind::SemanticHistory
        | IntegrationValidationKind::FaceTracking
        | IntegrationValidationKind::MergeReadiness => (
            IntegrationValidationState::RepairRequired("semantic-ownership-repair"),
            IntegrationValidationAction::SemanticCheck("semantic-ownership-repair"),
        ),
    }
}

fn repair(reason: &'static str) -> (IntegrationValidationState, IntegrationValidationAction) {
    (
        IntegrationValidationState::RepairRequired(reason),
        IntegrationValidationAction::Repair(reason),
    )
}
