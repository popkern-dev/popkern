//! Shared fact helpers for curve-backed trim reconstruction.

use crate::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepRepairFactId, BrepRepairFactKind,
    BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeSeed,
    BrepRepairTopologyHandle, HalfEdgeHandle, LocalToleranceSource,
};

pub(super) fn curve_trim_reconstruction_knowledge(
    brep: &BRep,
    options: BrepCurveTrimReconstructionOptions,
) -> BrepRepairKnowledgeGraph {
    BrepRepairKnowledgeGraph::build_with_local_tolerance(
        brep,
        BrepRepairKnowledgeScope::seeds_only(
            brep.half_edges
                .keys()
                .map(BrepRepairScopeSeed::HalfEdge)
                .collect::<Vec<_>>(),
        ),
        options.tolerance,
        LocalToleranceSource::RepairPolicy,
        options.samples,
    )
}

pub(super) fn half_edge_reconstruction_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Vec<BrepRepairFactId> {
    let mut ids = Vec::new();
    ids.extend(graph.fact_ids_for_handle(
        BrepRepairTopologyHandle::HalfEdge(half_edge),
        BrepRepairFactKind::TrimEndpointAgreement,
    ));
    ids.extend(graph.fact_ids_for_handle(
        BrepRepairTopologyHandle::HalfEdge(half_edge),
        BrepRepairFactKind::Authority,
    ));
    if let Some(face) = brep
        .half_edges
        .get(half_edge)
        .and_then(|he| brep.loops.get(he.loop_))
        .map(|loop_data| loop_data.face)
        && let Some(fact) = graph.surface_boundary_fact(face, half_edge)
    {
        ids.push(fact.header.id);
    }
    ids.sort();
    ids.dedup();
    ids
}

pub(super) fn shared_fact_ids_from_reconstruction_decisions(
    decisions: &[crate::topology::BrepCurveTrimReconstructionDecision],
) -> Vec<BrepRepairFactId> {
    let mut ids = decisions
        .iter()
        .flat_map(|decision| decision.shared_fact_ids.iter().copied())
        .collect::<Vec<_>>();
    ids.sort();
    ids.dedup();
    ids
}
