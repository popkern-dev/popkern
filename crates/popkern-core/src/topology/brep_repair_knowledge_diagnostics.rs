//! Diagnostic rendering for shared BRep repair knowledge facts.

use crate::topology::{
    BrepRepairAuthorityReason, BrepRepairDiagnostic, BrepRepairDiagnosticCode,
    BrepRepairDiagnosticId, BrepRepairDiagnosticSeverity, BrepRepairFact, BrepRepairFactKind,
    BrepRepairFactState, BrepRepairSuggestedRepair,
};

pub(super) fn build_diagnostic_records(facts: &[BrepRepairFact]) -> Vec<BrepRepairDiagnostic> {
    facts
        .iter()
        .enumerate()
        .map(|(index, fact)| BrepRepairDiagnostic {
            id: BrepRepairDiagnosticId(index),
            code: BrepRepairDiagnosticCode(fact_code(fact).to_string()),
            severity: fact_severity(fact.header.state),
            facts: vec![fact.header.id],
            handles: fact.handles.clone(),
            suggested_repair: fact_suggested_repair(fact.header.kind),
            render_key: fact_code(fact).to_string(),
            migration_label: fact.migration_label.clone(),
        })
        .collect()
}

fn fact_code(fact: &BrepRepairFact) -> &'static str {
    if let Some(authority) = &fact.authority {
        return authority_code(authority.reason);
    }
    match fact.header.kind {
        BrepRepairFactKind::Fingerprint => "brep.scope.fingerprint",
        BrepRepairFactKind::Loop => "brep.loop.facts",
        BrepRepairFactKind::LoopAdjacency => "brep.loop.adjacency",
        BrepRepairFactKind::EndpointMembership => "brep.endpoint.membership",
        BrepRepairFactKind::TrimEndpointAgreement => "brep.trim.endpoint_agreement",
        BrepRepairFactKind::SurfaceBoundary => "brep.surface.boundary",
        BrepRepairFactKind::EdgeSurfaceCandidate => "brep.edge_surface.candidate",
        BrepRepairFactKind::Authority => "brep.authority",
        BrepRepairFactKind::OperatorPolicy => "brep.operator_policy",
    }
}

fn authority_code(reason: BrepRepairAuthorityReason) -> &'static str {
    match reason {
        BrepRepairAuthorityReason::MissingPcurve => "brep.trim.pcurve_missing",
        BrepRepairAuthorityReason::StalePcurve => "brep.trim.pcurve_stale",
        BrepRepairAuthorityReason::SideDisagreement
        | BrepRepairAuthorityReason::SurfaceBoundaryTrimDisagreement => {
            "brep.trim.endpoint_side_disagreement"
        }
        BrepRepairAuthorityReason::MissingSurfaceBoundaryEvidence => {
            "brep.surface.boundary_hint_stale"
        }
        BrepRepairAuthorityReason::SurfaceBoundaryStaleHintCleared => {
            "brep.surface.boundary_hint_stale"
        }
        BrepRepairAuthorityReason::SurfaceBoundarySharedTrimAccepted
        | BrepRepairAuthorityReason::SurfaceBoundaryClassifiedHint => {
            "brep.surface.boundary_accepted"
        }
        BrepRepairAuthorityReason::SurfaceBoundaryClassifiedDecline => {
            "brep.surface.boundary_declined"
        }
        BrepRepairAuthorityReason::SurfaceBoundarySharedTrimMissingTopology
        | BrepRepairAuthorityReason::SurfaceBoundarySharedTrimHalfEdgeMismatch
        | BrepRepairAuthorityReason::SurfaceBoundarySharedTrimNotTwins => {
            "brep.surface.boundary_shared_trim_topology"
        }
        BrepRepairAuthorityReason::SurfaceBoundarySharedTrimEndpointMismatch => {
            "brep.surface.boundary_shared_trim_endpoint"
        }
        BrepRepairAuthorityReason::SurfaceBoundarySharedTrimCurveSurfaceMismatch => {
            "brep.surface.boundary_shared_trim_curve_surface"
        }
        BrepRepairAuthorityReason::SurfaceBoundarySharedTrimSideAgreementMismatch => {
            "brep.surface.boundary_shared_trim_side_agreement"
        }
        BrepRepairAuthorityReason::SurfaceBoundarySharedTrimValidationFailed => {
            "brep.surface.boundary_shared_trim_validation"
        }
        BrepRepairAuthorityReason::MissingTopology => "brep.scope.seed_stale",
        BrepRepairAuthorityReason::DegenerateZeroLengthVertexEdge
        | BrepRepairAuthorityReason::DegenerateAmbiguousCollapsedTrim
        | BrepRepairAuthorityReason::DegenerateCollapsedSurfaceTrim => "brep.edge.degenerate",
        BrepRepairAuthorityReason::MissingSurface => "brep.surface.missing",
        BrepRepairAuthorityReason::EvaluationFailed => "brep.evaluation.failed",
        _ => "brep.authority",
    }
}

fn fact_severity(state: BrepRepairFactState) -> BrepRepairDiagnosticSeverity {
    match state {
        BrepRepairFactState::Observed
        | BrepRepairFactState::Derived
        | BrepRepairFactState::Candidate
        | BrepRepairFactState::Decision => BrepRepairDiagnosticSeverity::Info,
        BrepRepairFactState::Stale | BrepRepairFactState::Contradictory => {
            BrepRepairDiagnosticSeverity::Warning
        }
        BrepRepairFactState::Unavailable => BrepRepairDiagnosticSeverity::RepairableError,
    }
}

fn fact_suggested_repair(kind: BrepRepairFactKind) -> BrepRepairSuggestedRepair {
    match kind {
        BrepRepairFactKind::TrimEndpointAgreement
        | BrepRepairFactKind::SurfaceBoundary
        | BrepRepairFactKind::EdgeSurfaceCandidate
        | BrepRepairFactKind::Authority => BrepRepairSuggestedRepair::TrimReconstruction,
        BrepRepairFactKind::Loop
        | BrepRepairFactKind::LoopAdjacency
        | BrepRepairFactKind::EndpointMembership => BrepRepairSuggestedRepair::Classification,
        BrepRepairFactKind::Fingerprint | BrepRepairFactKind::OperatorPolicy => {
            BrepRepairSuggestedRepair::None
        }
    }
}
