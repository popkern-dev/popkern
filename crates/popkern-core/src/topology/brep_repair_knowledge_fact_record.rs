//! Fact-record constructors for shared BRep repair knowledge.

use crate::topology::{
    BrepRepairAuthorityFact, BrepRepairFact, BrepRepairFactDependencies, BrepRepairFactHeader,
    BrepRepairFactId, BrepRepairFactKind, BrepRepairFactState, BrepRepairTopologyHandle,
    FaceHandle, HalfEdgeHandle, LocalToleranceEvidence,
};

pub(super) fn surface_boundary_label(face: FaceHandle, half_edge: HalfEdgeHandle) -> String {
    format!("{face:?}@{half_edge:?} surface_boundary")
}

pub(super) fn fact(
    id: usize,
    kind: BrepRepairFactKind,
    state: BrepRepairFactState,
    handles: Vec<BrepRepairTopologyHandle>,
    migration_label: Option<String>,
    tolerance: Option<LocalToleranceEvidence>,
) -> BrepRepairFact {
    authority_fact(id, kind, state, handles, migration_label, tolerance, None)
}

pub(super) fn authority_fact(
    id: usize,
    kind: BrepRepairFactKind,
    state: BrepRepairFactState,
    handles: Vec<BrepRepairTopologyHandle>,
    migration_label: Option<String>,
    tolerance: Option<LocalToleranceEvidence>,
    authority: Option<BrepRepairAuthorityFact>,
) -> BrepRepairFact {
    BrepRepairFact {
        header: BrepRepairFactHeader {
            id: BrepRepairFactId(id),
            kind,
            state,
            dependencies: BrepRepairFactDependencies::default(),
            provenance: Vec::new(),
            tolerance,
        },
        handles,
        migration_label,
        authority,
    }
}

pub(super) fn trim_fact_state(
    agreement: &crate::topology::BrepTrimEndpointAgreementFacts,
) -> BrepRepairFactState {
    use crate::topology::BrepTrimEndpointAgreementStatus as Status;
    match agreement.status {
        Status::AcceptedDirect | Status::AcceptedReconciled => BrepRepairFactState::Derived,
        Status::SideDisagreement => BrepRepairFactState::Contradictory,
        Status::MissingSurface
        | Status::MissingPcurve
        | Status::MissingTopology
        | Status::EvaluationFailed => BrepRepairFactState::Unavailable,
    }
}
