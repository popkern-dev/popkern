//! Shared fact-id helpers for edge/surface vertex-candidate diagnostics.

use crate::topology::{
    BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairTopologyHandle,
    HalfEdgeHandle,
};

pub(super) fn trim_endpoint_fact_ids(
    knowledge: &BrepRepairKnowledgeGraph,
    half_edge: HalfEdgeHandle,
    kind: BrepRepairFactKind,
) -> Vec<BrepRepairFactId> {
    let mut fact_ids =
        knowledge.fact_ids_for_handle(BrepRepairTopologyHandle::HalfEdge(half_edge), kind);
    fact_ids.extend(authority_fact_ids(knowledge, half_edge));
    fact_ids.sort();
    fact_ids.dedup();
    fact_ids
}

pub(super) fn surface_boundary_fact_ids(
    knowledge: &BrepRepairKnowledgeGraph,
    half_edge: HalfEdgeHandle,
    kind: BrepRepairFactKind,
) -> Vec<BrepRepairFactId> {
    let mut fact_ids =
        knowledge.fact_ids_for_handle(BrepRepairTopologyHandle::HalfEdge(half_edge), kind);
    fact_ids.extend(authority_fact_ids(knowledge, half_edge));
    fact_ids.sort();
    fact_ids.dedup();
    fact_ids
}

pub(super) fn edge_surface_candidate_fact_ids(
    knowledge: &BrepRepairKnowledgeGraph,
    half_edges: &[HalfEdgeHandle; 2],
) -> Vec<BrepRepairFactId> {
    let mut fact_ids = half_edges
        .iter()
        .flat_map(|half_edge| {
            let mut ids = knowledge.fact_ids_for_handle(
                BrepRepairTopologyHandle::HalfEdge(*half_edge),
                BrepRepairFactKind::EdgeSurfaceCandidate,
            );
            ids.extend(authority_fact_ids(knowledge, *half_edge));
            ids
        })
        .collect::<Vec<_>>();
    fact_ids.sort();
    fact_ids.dedup();
    fact_ids
}

fn authority_fact_ids(
    knowledge: &BrepRepairKnowledgeGraph,
    half_edge: HalfEdgeHandle,
) -> Vec<BrepRepairFactId> {
    knowledge.fact_ids_for_handle(
        BrepRepairTopologyHandle::HalfEdge(half_edge),
        BrepRepairFactKind::Authority,
    )
}
