//! Shared repair knowledge attribution for loop/shell diagnostics.

use std::collections::BTreeSet;

use crate::topology::{
    BRep, BrepLoopShellDiagnosticReport, BrepRepairFactId, BrepRepairFactKind,
    BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeExpansion,
    BrepRepairScopeSeed, BrepRepairTopologyHandle, BrepTopologyDiagnosticFinding,
};

pub(super) fn loop_shell_knowledge(brep: &BRep) -> BrepRepairKnowledgeGraph {
    let mut seeds = brep
        .faces
        .keys()
        .map(BrepRepairScopeSeed::Face)
        .chain(
            brep.boundary_edges
                .iter()
                .copied()
                .map(BrepRepairScopeSeed::Edge),
        )
        .collect::<Vec<_>>();
    seeds.sort();
    seeds.dedup();
    let mut scope = BrepRepairKnowledgeScope::seeds_only(seeds);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    BrepRepairKnowledgeGraph::build(brep, scope)
}

pub(super) fn attach_loop_shell_fact_ids(
    mut report: BrepLoopShellDiagnosticReport,
    brep: &BRep,
    knowledge: &BrepRepairKnowledgeGraph,
) -> BrepLoopShellDiagnosticReport {
    for finding in &mut report.findings {
        finding.shared_fact_ids = finding_fact_ids(knowledge, brep, finding);
    }
    report.shared_fact_ids = if report.findings.is_empty() {
        knowledge
            .facts
            .iter()
            .map(|fact| fact.header.id)
            .collect::<Vec<_>>()
    } else {
        report
            .findings
            .iter()
            .flat_map(|finding| finding.shared_fact_ids.iter().copied())
            .collect::<Vec<_>>()
    };
    report.shared_fact_ids.sort();
    report.shared_fact_ids.dedup();
    report
}

fn finding_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    brep: &BRep,
    finding: &BrepTopologyDiagnosticFinding,
) -> Vec<BrepRepairFactId> {
    let mut ids = BTreeSet::new();
    if let Some(shell) = finding.shell {
        for (face, data) in &brep.faces {
            if data.shell == shell {
                insert_fact_ids(graph, BrepRepairTopologyHandle::Face(face), &mut ids);
            }
        }
    }
    if let Some(face) = finding.face {
        insert_fact_ids(graph, BrepRepairTopologyHandle::Face(face), &mut ids);
    }
    if let Some(loop_) = finding.loop_ {
        insert_fact_ids(graph, BrepRepairTopologyHandle::Loop(loop_), &mut ids);
    }
    if let Some(edge) = finding.edge {
        insert_fact_ids(graph, BrepRepairTopologyHandle::Edge(edge), &mut ids);
    }
    if let Some(half_edge) = finding.half_edge {
        insert_fact_ids(
            graph,
            BrepRepairTopologyHandle::HalfEdge(half_edge),
            &mut ids,
        );
    }
    ids.into_iter().collect()
}

fn insert_fact_ids(
    graph: &BrepRepairKnowledgeGraph,
    handle: BrepRepairTopologyHandle,
    ids: &mut BTreeSet<BrepRepairFactId>,
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
