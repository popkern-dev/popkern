//! Shared fact helpers for missing edge-curve completion.

use crate::topology::{
    BRep, BrepEdgeCurveCompletionDecision, BrepEdgeCurveCompletionOptions, BrepRepairFactId,
    BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeSeed,
    BrepRepairTopologyHandle, EdgeHandle, LocalToleranceSource,
};

pub(super) fn edge_curve_completion_knowledge(
    brep: &BRep,
    options: BrepEdgeCurveCompletionOptions,
) -> BrepRepairKnowledgeGraph {
    BrepRepairKnowledgeGraph::build_with_local_tolerance(
        brep,
        BrepRepairKnowledgeScope::seeds_only(
            brep.edges
                .keys()
                .map(BrepRepairScopeSeed::Edge)
                .collect::<Vec<_>>(),
        ),
        options.tolerance,
        LocalToleranceSource::RepairPolicy,
        options.samples,
    )
}

pub(super) fn edge_completion_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    edge: EdgeHandle,
) -> Vec<BrepRepairFactId> {
    let mut ids = Vec::new();
    ids.extend(graph.fact_ids_for_handle(
        BrepRepairTopologyHandle::Edge(edge),
        BrepRepairFactKind::Authority,
    ));
    for half_edge in edge_half_edges(brep, edge) {
        if let Some(data) = brep.half_edges.get(half_edge) {
            ids.extend(graph.fact_ids_for_handle(
                BrepRepairTopologyHandle::Loop(data.loop_),
                BrepRepairFactKind::Loop,
            ));
            ids.extend(graph.fact_ids_for_handle(
                BrepRepairTopologyHandle::Vertex(data.origin),
                BrepRepairFactKind::EndpointMembership,
            ));
            if let Some(next) = brep.half_edges.get(data.next) {
                ids.extend(graph.fact_ids_for_handle(
                    BrepRepairTopologyHandle::Vertex(next.origin),
                    BrepRepairFactKind::EndpointMembership,
                ));
            }
        }
        ids.extend(graph.fact_ids_for_handle(
            BrepRepairTopologyHandle::HalfEdge(half_edge),
            BrepRepairFactKind::TrimEndpointAgreement,
        ));
        ids.extend(graph.fact_ids_for_handle(
            BrepRepairTopologyHandle::HalfEdge(half_edge),
            BrepRepairFactKind::EdgeSurfaceCandidate,
        ));
    }
    ids.sort();
    ids.dedup();
    ids
}

pub(super) fn shared_fact_ids_from_edge_completion_decisions(
    decisions: &[BrepEdgeCurveCompletionDecision],
) -> Vec<BrepRepairFactId> {
    let mut ids = decisions
        .iter()
        .flat_map(|decision| decision.shared_fact_ids.iter().copied())
        .collect::<Vec<_>>();
    ids.sort();
    ids.dedup();
    ids
}

fn edge_half_edges(brep: &BRep, edge: EdgeHandle) -> Vec<crate::topology::HalfEdgeHandle> {
    brep.half_edges
        .iter()
        .filter_map(|(half_edge, data)| (data.edge == edge).then_some(half_edge))
        .collect()
}
