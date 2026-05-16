//! Shared repair-knowledge fact plumbing for edge/surface compatibility reports.

use std::collections::BTreeSet;

use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairReport,
    BrepEdgeSurfaceCompatibilityReport,
};
use crate::topology::brep_edge_surface_compatibility_helpers::unique_edge_half_edges;
use crate::topology::{
    BRep, BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairTopologyHandle, EdgeHandle,
};

pub(super) fn edge_surface_repair_knowledge(
    brep: &BRep,
    report: &BrepEdgeSurfaceCompatibilityReport,
) -> BrepRepairKnowledgeGraph {
    let mut edges = report
        .findings
        .iter()
        .map(|finding| finding.edge)
        .collect::<Vec<_>>();
    edges.sort_by_key(|edge| format!("{edge:?}"));
    edges.dedup();
    let mut scope = BrepRepairKnowledgeScope::seeds_only(
        edges.into_iter().map(BrepRepairScopeSeed::Edge).collect(),
    );
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    BrepRepairKnowledgeGraph::build(brep, scope)
}

pub(super) fn with_shared_fact_ids(
    mut decision: BrepEdgeSurfaceCompatibilityRepairDecision,
    shared_fact_ids: Vec<BrepRepairFactId>,
) -> BrepEdgeSurfaceCompatibilityRepairDecision {
    decision.shared_fact_ids = shared_fact_ids;
    decision
}

pub(super) fn edge_repair_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    edge: EdgeHandle,
) -> Vec<BrepRepairFactId> {
    let mut fact_ids = Vec::new();
    if let Some(edge_data) = brep.edges.get(edge) {
        for half_edge in unique_edge_half_edges(brep, edge_data.half_edge) {
            fact_ids.extend(graph.fact_ids_for_handle(
                BrepRepairTopologyHandle::HalfEdge(half_edge),
                BrepRepairFactKind::SurfaceBoundary,
            ));
            fact_ids.extend(graph.fact_ids_for_handle(
                BrepRepairTopologyHandle::HalfEdge(half_edge),
                BrepRepairFactKind::TrimEndpointAgreement,
            ));
            fact_ids.extend(graph.fact_ids_for_handle(
                BrepRepairTopologyHandle::HalfEdge(half_edge),
                BrepRepairFactKind::EdgeSurfaceCandidate,
            ));
            fact_ids.extend(graph.fact_ids_for_handle(
                BrepRepairTopologyHandle::HalfEdge(half_edge),
                BrepRepairFactKind::Authority,
            ));
        }
    }
    fact_ids.sort();
    fact_ids.dedup();
    fact_ids
}

pub(super) fn shared_fact_ids_from_report(
    report: &BrepEdgeSurfaceCompatibilityRepairReport,
) -> Vec<BrepRepairFactId> {
    let fact_ids = report
        .vertex_stars
        .decisions
        .iter()
        .flat_map(|decision| decision.shared_fact_ids.iter().copied())
        .chain(
            report
                .decisions
                .iter()
                .flat_map(|decision| decision.shared_fact_ids.iter().copied()),
        )
        .collect::<BTreeSet<_>>();
    fact_ids.into_iter().collect()
}
