//! Shared BRep repair fact attribution for cylinder-mediated support planning.

use std::collections::BTreeSet;

use crate::euler::CylinderMediatedSupportContinuationRequest;
use crate::topology::{
    BRep, BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairTopologyHandle, HalfEdgeHandle,
};

pub(super) fn request_source_fact_ids(
    brep: &BRep,
    request: &CylinderMediatedSupportContinuationRequest,
) -> Vec<BrepRepairFactId> {
    let graph = request_knowledge(brep, request);
    let handles = request_handles(request);
    fact_ids_for_handles(&graph, &handles)
}

pub(super) fn batch_source_fact_ids(plan_ids: &[Vec<BrepRepairFactId>]) -> Vec<BrepRepairFactId> {
    let mut ids = plan_ids
        .iter()
        .flat_map(|ids| ids.iter().copied())
        .collect::<Vec<_>>();
    ids.sort();
    ids.dedup();
    ids
}

pub(super) fn segment_source_fact_ids(
    brep: &BRep,
    half_edges: &[HalfEdgeHandle],
) -> Vec<BrepRepairFactId> {
    let mut seeds = half_edges
        .iter()
        .copied()
        .map(BrepRepairScopeSeed::HalfEdge)
        .collect::<Vec<_>>();
    seeds.sort();
    seeds.dedup();
    let mut scope = BrepRepairKnowledgeScope::seeds_only(seeds);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    let graph = BrepRepairKnowledgeGraph::build(brep, scope);
    let handles = segment_handles(brep, half_edges);
    fact_ids_for_handles(&graph, &handles)
}

fn request_knowledge(
    brep: &BRep,
    request: &CylinderMediatedSupportContinuationRequest,
) -> BrepRepairKnowledgeGraph {
    let mut seeds = vec![
        BrepRepairScopeSeed::Face(request.support_face),
        BrepRepairScopeSeed::Face(request.outside_face),
        BrepRepairScopeSeed::Face(request.cylinder_face),
        BrepRepairScopeSeed::Loop(request.outside_loop),
        BrepRepairScopeSeed::Loop(request.support_loop),
        BrepRepairScopeSeed::HalfEdge(request.outside_start_half_edge),
        BrepRepairScopeSeed::HalfEdge(request.outside_end_half_edge),
        BrepRepairScopeSeed::HalfEdge(request.support_cylinder_half_edge),
        BrepRepairScopeSeed::HalfEdge(request.outside_cylinder_half_edge),
        BrepRepairScopeSeed::Vertex(request.start),
        BrepRepairScopeSeed::Vertex(request.end),
    ];
    seeds.extend(
        request
            .cylinder_interval
            .iter()
            .copied()
            .map(BrepRepairScopeSeed::HalfEdge),
    );
    seeds.extend(
        request
            .context_owner_faces
            .iter()
            .copied()
            .map(BrepRepairScopeSeed::Face),
    );
    seeds.sort();
    seeds.dedup();
    let mut scope = BrepRepairKnowledgeScope::seeds_only(seeds);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    BrepRepairKnowledgeGraph::build(brep, scope)
}

fn request_handles(
    request: &CylinderMediatedSupportContinuationRequest,
) -> Vec<BrepRepairTopologyHandle> {
    let mut handles = vec![
        BrepRepairTopologyHandle::Face(request.support_face),
        BrepRepairTopologyHandle::Face(request.outside_face),
        BrepRepairTopologyHandle::Face(request.cylinder_face),
        BrepRepairTopologyHandle::Loop(request.outside_loop),
        BrepRepairTopologyHandle::Loop(request.support_loop),
        BrepRepairTopologyHandle::HalfEdge(request.outside_start_half_edge),
        BrepRepairTopologyHandle::HalfEdge(request.outside_end_half_edge),
        BrepRepairTopologyHandle::HalfEdge(request.support_cylinder_half_edge),
        BrepRepairTopologyHandle::HalfEdge(request.outside_cylinder_half_edge),
        BrepRepairTopologyHandle::Vertex(request.start),
        BrepRepairTopologyHandle::Vertex(request.end),
    ];
    handles.extend(
        request
            .cylinder_interval
            .iter()
            .copied()
            .map(BrepRepairTopologyHandle::HalfEdge),
    );
    handles.extend(
        request
            .context_owner_faces
            .iter()
            .copied()
            .map(BrepRepairTopologyHandle::Face),
    );
    handles
}

fn segment_handles(brep: &BRep, half_edges: &[HalfEdgeHandle]) -> Vec<BrepRepairTopologyHandle> {
    let mut handles = Vec::new();
    for &half_edge in half_edges {
        handles.push(BrepRepairTopologyHandle::HalfEdge(half_edge));
        let Some(data) = brep.half_edges.get(half_edge) else {
            continue;
        };
        handles.push(BrepRepairTopologyHandle::Edge(data.edge));
        handles.push(BrepRepairTopologyHandle::Vertex(data.origin));
        handles.push(BrepRepairTopologyHandle::Loop(data.loop_));
        if let Some(loop_data) = brep.loops.get(data.loop_) {
            handles.push(BrepRepairTopologyHandle::Face(loop_data.face));
        }
    }
    handles
}

fn fact_ids_for_handles(
    graph: &BrepRepairKnowledgeGraph,
    handles: &[BrepRepairTopologyHandle],
) -> Vec<BrepRepairFactId> {
    let mut ids = BTreeSet::new();
    for &handle in handles {
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
    ids.into_iter().collect()
}
