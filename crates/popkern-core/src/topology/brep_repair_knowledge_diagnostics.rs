//! Diagnostic rendering for shared BRep repair knowledge facts.

use crate::topology::{
    BrepRepairDiagnostic, BrepRepairDiagnosticCode, BrepRepairDiagnosticId,
    BrepRepairDiagnosticSeverity, BrepRepairFact, BrepRepairFactKind, BrepRepairFactState,
    BrepRepairSuggestedRepair,
};

pub(super) fn build_diagnostic_records(facts: &[BrepRepairFact]) -> Vec<BrepRepairDiagnostic> {
    facts
        .iter()
        .enumerate()
        .map(|(index, fact)| BrepRepairDiagnostic {
            id: BrepRepairDiagnosticId(index),
            code: BrepRepairDiagnosticCode(fact_code(fact.header.kind).to_string()),
            severity: fact_severity(fact.header.state),
            facts: vec![fact.header.id],
            handles: fact.handles.clone(),
            suggested_repair: fact_suggested_repair(fact.header.kind),
            render_key: fact_code(fact.header.kind).to_string(),
            legacy_label: fact.legacy_label.clone(),
        })
        .collect()
}

fn fact_code(kind: BrepRepairFactKind) -> &'static str {
    match kind {
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
