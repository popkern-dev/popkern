//! Shared repair knowledge attribution for sewing operators.

use crate::topology::{
    BRep, BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairTopologyHandle, BrepSewingEdgePair,
    EdgeHandle, LocalToleranceSource,
};

pub(super) fn sewing_pair_knowledge(
    brep: &BRep,
    pairs: &[BrepSewingEdgePair],
    tolerance: f64,
    samples: usize,
) -> BrepRepairKnowledgeGraph {
    let mut seeds = pairs
        .iter()
        .flat_map(|pair| [pair.first, pair.second])
        .map(BrepRepairScopeSeed::Edge)
        .collect::<Vec<_>>();
    seeds.sort();
    seeds.dedup();
    let mut scope = BrepRepairKnowledgeScope::seeds_only(seeds);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    BrepRepairKnowledgeGraph::build_with_local_tolerance(
        brep,
        scope,
        tolerance,
        LocalToleranceSource::RepairPolicy,
        samples,
    )
}

pub(super) fn sewing_pair_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    pair: &BrepSewingEdgePair,
) -> Vec<BrepRepairFactId> {
    let mut ids = Vec::new();
    append_edge_fact_ids(graph, brep, pair.first, &mut ids);
    append_edge_fact_ids(graph, brep, pair.second, &mut ids);
    ids.sort();
    ids.dedup();
    ids
}

pub(super) fn shared_fact_ids_from_pairs(
    pairs: impl Iterator<Item = Vec<BrepRepairFactId>>,
) -> Vec<BrepRepairFactId> {
    let mut ids = pairs.flatten().collect::<Vec<_>>();
    ids.sort();
    ids.dedup();
    ids
}

fn append_edge_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    edge: EdgeHandle,
    ids: &mut Vec<BrepRepairFactId>,
) {
    append_handle_fact_ids(graph, BrepRepairTopologyHandle::Edge(edge), ids);
    let Some(edge_data) = brep.edges.get(edge) else {
        return;
    };
    append_half_edge_fact_ids(graph, brep, edge_data.half_edge, ids);
    if let Some(half_edge) = brep.half_edges.get(edge_data.half_edge) {
        append_half_edge_fact_ids(graph, brep, half_edge.twin, ids);
    }
}

fn append_half_edge_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    half_edge: crate::topology::HalfEdgeHandle,
    ids: &mut Vec<BrepRepairFactId>,
) {
    let Some(data) = brep.half_edges.get(half_edge) else {
        return;
    };
    append_handle_fact_ids(graph, BrepRepairTopologyHandle::HalfEdge(half_edge), ids);
    append_handle_fact_ids(graph, BrepRepairTopologyHandle::Loop(data.loop_), ids);
    append_handle_fact_ids(graph, BrepRepairTopologyHandle::Vertex(data.origin), ids);
}

fn append_handle_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    handle: BrepRepairTopologyHandle,
    ids: &mut Vec<BrepRepairFactId>,
) {
    for kind in [
        BrepRepairFactKind::Loop,
        BrepRepairFactKind::LoopAdjacency,
        BrepRepairFactKind::EndpointMembership,
        BrepRepairFactKind::TrimEndpointAgreement,
        BrepRepairFactKind::SurfaceBoundary,
        BrepRepairFactKind::EdgeSurfaceCandidate,
        BrepRepairFactKind::Authority,
    ] {
        ids.extend(graph.fact_ids_for_handle(handle, kind));
    }
}
