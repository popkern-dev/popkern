//! Shared fact-id helpers for face-use trim authority diagnostics.

use crate::topology::{
    BRep, BrepFaceUseTrimAuthorityFinding, BrepFaceUseTrimAuthorityOptions, BrepRepairFactId,
    BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairTopologyHandle, LocalToleranceSource,
};

pub(super) fn face_use_trim_authority_knowledge(
    brep: &BRep,
    options: BrepFaceUseTrimAuthorityOptions,
) -> BrepRepairKnowledgeGraph {
    let mut scope = BrepRepairKnowledgeScope::seeds_only(
        brep.faces
            .keys()
            .map(BrepRepairScopeSeed::Face)
            .collect::<Vec<_>>(),
    );
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    BrepRepairKnowledgeGraph::build_with_local_tolerance(
        brep,
        scope,
        options.tolerance,
        LocalToleranceSource::RepairPolicy,
        options.samples,
    )
}

pub(super) fn attach_face_use_trim_authority_fact_ids(
    brep: &BRep,
    graph: &BrepRepairKnowledgeGraph,
    findings: &mut [BrepFaceUseTrimAuthorityFinding],
) -> Vec<BrepRepairFactId> {
    let mut aggregate = Vec::new();
    for finding in findings {
        finding.shared_fact_ids = finding_fact_ids(brep, graph, finding);
        aggregate.extend(finding.shared_fact_ids.iter().copied());
    }
    aggregate.sort();
    aggregate.dedup();
    aggregate
}

fn finding_fact_ids(
    brep: &BRep,
    graph: &BrepRepairKnowledgeGraph,
    finding: &BrepFaceUseTrimAuthorityFinding,
) -> Vec<BrepRepairFactId> {
    let mut ids = Vec::new();
    append_handle_fact_ids(
        graph,
        BrepRepairTopologyHandle::Face(finding.face),
        &mut ids,
    );
    append_handle_fact_ids(
        graph,
        BrepRepairTopologyHandle::Loop(finding.loop_),
        &mut ids,
    );
    append_handle_fact_ids(
        graph,
        BrepRepairTopologyHandle::HalfEdge(finding.half_edge),
        &mut ids,
    );
    append_handle_fact_ids(
        graph,
        BrepRepairTopologyHandle::Edge(finding.edge),
        &mut ids,
    );
    append_edge_endpoint_fact_ids(brep, graph, finding.half_edge, &mut ids);
    if let Some(surface_boundary) = graph.surface_boundary_fact(finding.face, finding.half_edge) {
        ids.push(surface_boundary.header.id);
    }
    ids.sort();
    ids.dedup();
    ids
}

fn append_edge_endpoint_fact_ids(
    brep: &BRep,
    graph: &BrepRepairKnowledgeGraph,
    half_edge: crate::topology::HalfEdgeHandle,
    ids: &mut Vec<BrepRepairFactId>,
) {
    let Some(data) = brep.half_edges.get(half_edge) else {
        return;
    };
    append_handle_fact_ids(graph, BrepRepairTopologyHandle::Vertex(data.origin), ids);
    if let Some(next) = brep.half_edges.get(data.next) {
        append_handle_fact_ids(graph, BrepRepairTopologyHandle::Vertex(next.origin), ids);
    }
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
