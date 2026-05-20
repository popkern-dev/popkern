use popkern_core::topology::{
    BRep, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeExpansion,
    BrepRepairScopeSeed, BrepRepairTopologyHandle, EdgeHandle, TopologyRole,
    classify_brep_fidelity_blocker_ownership,
};

#[derive(Debug, Clone, PartialEq)]
struct SurfaceBoundaryKnowledgeEvidence {
    edge: String,
    diagnostics: Vec<String>,
}

#[test]
fn gregory_surface_boundary_knowledge_covers_remaining_blockers() {
    let construction = super::mixed_radius_gregory_cube();
    let ownership = classify_brep_fidelity_blocker_ownership(&construction, 1e-6);
    let blocker_edges = ownership
        .edges
        .iter()
        .filter_map(|edge| {
            (edge.edge_role == TopologyRole::Production
                && edge
                    .half_edges
                    .iter()
                    .any(|side| side.face_role == TopologyRole::Production))
            .then_some(edge.edge)
        })
        .collect::<std::collections::BTreeSet<_>>();

    let evidence = blocker_edges
        .iter()
        .copied()
        .map(|edge| surface_boundary_knowledge_evidence(&construction, edge))
        .collect::<Vec<_>>();

    assert!(
        !evidence.is_empty(),
        "diagnostic must cover at least one production-owned pcurve blocker: {ownership:#?}"
    );
    assert!(
        evidence.iter().all(|entry| entry
            .diagnostics
            .iter()
            .any(|code| code == "brep.surface.boundary_hint_stale")),
        "surface-boundary repair knowledge should classify stale boundary hints for each remaining production-owned pcurve blocker: {evidence:#?}"
    );
    assert!(
        evidence.iter().all(|entry| entry
            .diagnostics
            .iter()
            .any(|code| code == "brep.trim.pcurve_stale")),
        "surface-boundary repair knowledge should classify stale pcurves for each remaining production-owned pcurve blocker: {evidence:#?}"
    );
}

fn surface_boundary_knowledge_evidence(
    brep: &BRep,
    edge: EdgeHandle,
) -> SurfaceBoundaryKnowledgeEvidence {
    let mut scope = BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Edge(edge)]);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    let graph = BrepRepairKnowledgeGraph::build(brep, scope);
    let edge_handle = BrepRepairTopologyHandle::Edge(edge);
    let diagnostics = graph
        .diagnostics_for_handle(edge_handle)
        .iter()
        .map(|diagnostic| diagnostic.code.0.clone())
        .collect::<Vec<_>>();
    SurfaceBoundaryKnowledgeEvidence {
        edge: format!("{edge:?}"),
        diagnostics,
    }
}
