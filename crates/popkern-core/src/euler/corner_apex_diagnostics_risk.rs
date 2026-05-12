//! Construction-risk classifiers for corner/apex diagnostics.
#![allow(missing_docs)]

use super::risk_types::ConstructionRiskAction;
use super::types::{ConstructionRiskKind, ConstructionRiskState};

pub(super) fn construction_risk_state_action(
    kind: ConstructionRiskKind,
    primitive_available: bool,
) -> (ConstructionRiskState, ConstructionRiskAction) {
    if primitive_available {
        return (ConstructionRiskState::NoRisk, ConstructionRiskAction::None);
    }
    match kind {
        ConstructionRiskKind::IntersectionDegenerate => (
            ConstructionRiskState::ManagedRetry("intersection-refinement"),
            ConstructionRiskAction::Retry("intersection-refinement"),
        ),
        ConstructionRiskKind::TangencyUnsolved => (
            ConstructionRiskState::RequiresMissingPrimitive("three-face-tangency-solver"),
            ConstructionRiskAction::MissingPrimitive("three-face-tangency-solver"),
        ),
        ConstructionRiskKind::PatchBoundaryIncompatible => (
            ConstructionRiskState::RequiresMissingPrimitive("patch-boundary-repair"),
            ConstructionRiskAction::MissingPrimitive("patch-boundary-repair"),
        ),
        ConstructionRiskKind::CurvatureExplosion => (
            ConstructionRiskState::ManagedRetry("continuity-relaxation"),
            ConstructionRiskAction::Reclassify("continuity-relaxation"),
        ),
        ConstructionRiskKind::SelfIntersection => (
            ConstructionRiskState::Invalid("self-intersection"),
            ConstructionRiskAction::Invalid("self-intersection"),
        ),
        ConstructionRiskKind::VariableRadiusCollapse => (
            ConstructionRiskState::Invalid("variable-radius-collapse"),
            ConstructionRiskAction::Invalid("variable-radius-collapse"),
        ),
        ConstructionRiskKind::NumericalSolverFailure => (
            ConstructionRiskState::RequiresMissingPrimitive("numerical-solver-failure"),
            ConstructionRiskAction::RollbackRequired("transactional-construction-rollback"),
        ),
        _ => missing_primitive_action(kind),
    }
}

fn missing_primitive_action(
    kind: ConstructionRiskKind,
) -> (ConstructionRiskState, ConstructionRiskAction) {
    let primitive = match kind {
        ConstructionRiskKind::BoundaryCurveMissing => "boundary-curve-extraction",
        ConstructionRiskKind::SetbackSplitInfeasible => "setback-split-construction",
        ConstructionRiskKind::BoundaryConditionsIncomplete => "boundary-condition-assembly",
        ConstructionRiskKind::TwistCompatibility => "twist-compatibility-solver",
        ConstructionRiskKind::ChamferRoundMeshing => "mixed-chamfer-round-meshing",
        ConstructionRiskKind::FacetedSmoothMismatch => "geometry-mode-conversion",
        ConstructionRiskKind::ApproximationResidual => "approximation-refinement",
        _ => "construction-risk-refinement",
    };
    (
        ConstructionRiskState::RequiresMissingPrimitive(primitive),
        ConstructionRiskAction::MissingPrimitive(primitive),
    )
}
