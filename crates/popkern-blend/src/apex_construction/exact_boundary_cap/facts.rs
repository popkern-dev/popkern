//! Shared BRep repair-knowledge attribution for exact-boundary cap diagnostics.

use std::collections::BTreeSet;

use popkern_core::topology::{
    BRep, BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairTopologyHandle,
};

use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

/// Attach shared repair fact ids to exact-boundary diagnostics that reference
/// live source topology.
pub(crate) fn attach_shared_repair_fact_ids(
    brep: &BRep,
    diagnostics: &mut [ApexConstructionDiagnostic],
) {
    let seeds = exact_boundary_scope_seeds(diagnostics);
    if seeds.is_empty() {
        return;
    }
    let mut scope = BrepRepairKnowledgeScope::seeds_only(seeds);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    let graph = BrepRepairKnowledgeGraph::build(brep, scope);
    for diagnostic in diagnostics {
        let ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) = diagnostic else {
            continue;
        };
        diagnostic.shared_fact_ids = diagnostic_fact_ids(&graph, diagnostic);
    }
}

fn exact_boundary_scope_seeds(
    diagnostics: &[ApexConstructionDiagnostic],
) -> Vec<BrepRepairScopeSeed> {
    let mut seeds = BTreeSet::new();
    for diagnostic in diagnostics {
        let ApexConstructionDiagnostic::ExactBoundaryCap(diagnostic) = diagnostic else {
            continue;
        };
        if let Some(edge) = diagnostic.edge {
            seeds.insert(BrepRepairScopeSeed::Edge(edge));
        }
        if let Some(half_edge) = diagnostic.support_half_edge {
            seeds.insert(BrepRepairScopeSeed::HalfEdge(half_edge));
        }
        if let Some(face) = diagnostic.support_face {
            seeds.insert(BrepRepairScopeSeed::Face(face));
        }
    }
    seeds.into_iter().collect()
}

fn diagnostic_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    diagnostic: &ExactBoundaryCapDiagnostic,
) -> Vec<BrepRepairFactId> {
    let mut ids = Vec::new();
    if let Some(edge) = diagnostic.edge {
        append_fact_ids(graph, BrepRepairTopologyHandle::Edge(edge), &mut ids);
    }
    if let Some(half_edge) = diagnostic.support_half_edge {
        append_fact_ids(
            graph,
            BrepRepairTopologyHandle::HalfEdge(half_edge),
            &mut ids,
        );
    }
    if let Some(face) = diagnostic.support_face {
        append_fact_ids(graph, BrepRepairTopologyHandle::Face(face), &mut ids);
    }
    ids.sort();
    ids.dedup();
    ids
}

fn append_fact_ids(
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
        BrepRepairFactKind::OperatorPolicy,
    ] {
        ids.extend(graph.fact_ids_for_handle(handle, kind));
    }
}
