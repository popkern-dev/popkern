use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepEdgeSurfaceKind, BrepRepairFactKind,
    BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeExpansion,
    BrepRepairScopeSeed, BrepRepairTopologyHandle, EdgeHandle,
    analyze_brep_edge_surface_compatibility,
};

#[derive(Debug, Clone, PartialEq)]
struct SurfaceBoundaryKnowledgeEvidence {
    edge: String,
    diagnostics: Vec<String>,
    audit_facts: Vec<String>,
}

#[test]
fn gregory_surface_boundary_knowledge_covers_remaining_blockers() {
    let construction = super::mixed_radius_gregory_cube();
    let edge_surface = analyze_brep_edge_surface_compatibility(
        &construction,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: 1e-6,
            max_anchor_motion: 0.0,
            samples: 9,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    let blocker_edges = edge_surface
        .findings
        .iter()
        .filter_map(|finding| {
            let kinds = adjacent_surface_kinds(&construction, finding.edge)?;
            matches!(
                kinds,
                (BrepEdgeSurfaceKind::Gregory, BrepEdgeSurfaceKind::Cylinder)
                    | (BrepEdgeSurfaceKind::Cylinder, BrepEdgeSurfaceKind::Gregory)
            )
            .then_some(finding.edge)
        })
        .collect::<std::collections::BTreeSet<_>>();

    let evidence = blocker_edges
        .iter()
        .copied()
        .map(|edge| surface_boundary_knowledge_evidence(&construction, edge))
        .collect::<Vec<_>>();
    if std::env::var_os("POPKERN_GREGORY_FIDELITY_TRACE").is_some() {
        eprintln!("surface-boundary knowledge evidence: {evidence:#?}");
    }

    assert!(
        !evidence.is_empty(),
        "diagnostic must cover at least one Gregory/cylinder edge/surface blocker"
    );
    assert!(
        evidence.iter().all(|entry| !entry.audit_facts.is_empty()),
        "surface-boundary repair knowledge should expose durable audit facts for each remaining Gregory/cylinder blocker: {evidence:#?}"
    );
}

fn adjacent_surface_kinds(
    brep: &BRep,
    edge: EdgeHandle,
) -> Option<(BrepEdgeSurfaceKind, BrepEdgeSurfaceKind)> {
    let edge_data = brep.edges.get(edge)?;
    let first = edge_data.half_edge;
    let second = brep.half_edges.get(first)?.twin;
    Some((
        surface_kind_for_half_edge(brep, first)?,
        surface_kind_for_half_edge(brep, second)?,
    ))
}

fn surface_kind_for_half_edge(
    brep: &BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
) -> Option<BrepEdgeSurfaceKind> {
    let face = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face)?;
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    Some(match &surface.geometry {
        SurfaceGeometry::Plane { .. } => BrepEdgeSurfaceKind::Plane,
        SurfaceGeometry::Cylinder { .. } => BrepEdgeSurfaceKind::Cylinder,
        SurfaceGeometry::Cone { .. } => BrepEdgeSurfaceKind::Cone,
        SurfaceGeometry::Sphere { .. } => BrepEdgeSurfaceKind::Sphere,
        SurfaceGeometry::Torus { .. } => BrepEdgeSurfaceKind::Torus,
        SurfaceGeometry::Nurbs(_) => BrepEdgeSurfaceKind::Nurbs,
        SurfaceGeometry::Gregory(_) => BrepEdgeSurfaceKind::Gregory,
    })
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
    let audit_facts = graph
        .facts
        .iter()
        .filter(|fact| {
            fact.header.kind == BrepRepairFactKind::SurfaceBoundary
                && fact.handles.contains(&edge_handle)
                && fact
                    .legacy_label
                    .as_deref()
                    .is_some_and(|label| label.contains("surface_boundary_audit"))
        })
        .filter_map(|fact| fact.legacy_label.clone())
        .collect::<Vec<_>>();
    SurfaceBoundaryKnowledgeEvidence {
        edge: format!("{edge:?}"),
        diagnostics,
        audit_facts,
    }
}
