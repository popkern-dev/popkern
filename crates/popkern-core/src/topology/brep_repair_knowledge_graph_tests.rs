use super::*;
use crate::tolerance::Tolerance;

#[test]
fn builds_queryable_graph_from_face_scope() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let face = brep.faces.keys().next().expect("cube has face");
    let scope = BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Face(face)]);

    let graph = BrepRepairKnowledgeGraph::build(&brep, scope);

    assert!(graph.validate_fresh(&brep));
    assert_eq!(
        graph.fingerprint,
        BrepKnowledgeFingerprint::from_brep(&brep)
    );
    assert!(graph.face_loops(face).next().is_some());
    assert!(!graph.facts.is_empty());
    assert!(
        graph
            .facts
            .iter()
            .any(|fact| fact.header.kind == BrepRepairFactKind::Loop)
    );
}

#[test]
fn adjacent_face_scope_exposes_endpoint_and_trim_queries() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let face = brep.faces.keys().next().expect("cube has face");
    let mut scope = BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Face(face)]);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;

    let graph = BrepRepairKnowledgeGraph::build(&brep, scope);
    let loop_facts = graph.face_loops(face).next().expect("face loop fact");
    let vertex = loop_facts
        .boundary_vertices()
        .first()
        .copied()
        .expect("loop has vertex");
    let adjacency = graph
        .loop_adjacencies
        .first()
        .expect("expanded cube face has loop adjacency");

    assert!(graph.endpoint_membership(vertex).is_some());
    assert!(graph.trim_endpoint_agreement(adjacency.half_edge).is_some());
    assert!(
        graph
            .trim_endpoint_agreement_for_face(adjacency.half_edge, adjacency.face)
            .is_some()
    );
    assert!(
        graph
            .surface_boundary_fact(adjacency.face, adjacency.half_edge)
            .is_some()
    );
    assert!(graph.edge_surface_candidate(adjacency.half_edge).is_some());
    let authority = graph.authority_for_edge_side(adjacency.half_edge);
    assert!(!authority.is_empty());
    assert!(authority.iter().all(|fact| fact.authority.is_some()));
    let edge = brep.half_edges[adjacency.half_edge].edge;
    let edge_face_authority = graph.authority_for_edge_face(edge, adjacency.face);
    assert!(!edge_face_authority.is_empty());
    assert!(edge_face_authority.iter().all(|fact| {
        fact.handles.contains(&BrepRepairTopologyHandle::Edge(edge))
            && fact
                .handles
                .contains(&BrepRepairTopologyHandle::HalfEdge(adjacency.half_edge))
    }));
    assert!(
        !graph
            .diagnostics_for_handle(BrepRepairTopologyHandle::HalfEdge(adjacency.half_edge))
            .is_empty()
    );
    assert!(
        !graph
            .fact_ids_for_handle(
                BrepRepairTopologyHandle::HalfEdge(adjacency.half_edge),
                BrepRepairFactKind::TrimEndpointAgreement,
            )
            .is_empty()
    );
    assert!(
        graph
            .facts
            .iter()
            .any(|fact| fact.header.kind == BrepRepairFactKind::TrimEndpointAgreement)
    );
    assert!(
        graph
            .facts
            .iter()
            .any(|fact| fact.header.kind == BrepRepairFactKind::SurfaceBoundary)
    );
    assert!(
        graph
            .facts
            .iter()
            .any(|fact| fact.header.kind == BrepRepairFactKind::EdgeSurfaceCandidate)
    );
}

#[test]
fn explicit_local_tolerance_controls_trim_agreement_state() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let face = brep.faces.keys().next().expect("cube has face");
    let mut scope = BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Face(face)]);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    let adjacency = BrepRepairKnowledgeGraph::build(&brep, scope.clone())
        .loop_adjacencies
        .first()
        .cloned()
        .expect("expanded cube face has loop adjacency");
    let half_edge = adjacency.half_edge;
    let pcurve = brep.half_edges[half_edge]
        .pcurve
        .as_ref()
        .expect("cube half-edge has pcurve");
    let (t_start, t_end) = pcurve.parameter_range();
    let start = pcurve.point_at(t_start);
    let end = pcurve.point_at(t_end);
    brep.half_edges[half_edge].pcurve = Some(crate::geometry::pcurve::Pcurve::linear(
        (start.x + 1e-7, start.y),
        (end.x + 1e-7, end.y),
        t_start,
        t_end,
    ));

    let model_graph = BrepRepairKnowledgeGraph::build(&brep, scope.clone());
    let strict_graph = BrepRepairKnowledgeGraph::build_with_local_tolerance(
        &brep,
        scope,
        1e-9,
        LocalToleranceSource::RepairPolicy,
        9,
    );

    assert_ne!(
        model_graph
            .trim_endpoint_agreement(half_edge)
            .expect("model agreement exists")
            .status,
        crate::topology::BrepTrimEndpointAgreementStatus::SideDisagreement
    );
    assert_eq!(
        strict_graph
            .trim_endpoint_agreement(half_edge)
            .expect("strict agreement exists")
            .status,
        crate::topology::BrepTrimEndpointAgreementStatus::SideDisagreement
    );
}

#[test]
fn surface_boundary_audit_events_are_promoted_to_queryable_facts() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let face = brep.faces.keys().next().expect("cube has face");
    let mut scope = BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Face(face)]);
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    let adjacency = BrepRepairKnowledgeGraph::build(&brep, scope.clone())
        .loop_adjacencies
        .first()
        .cloned()
        .expect("expanded cube face has loop adjacency");
    let edge = brep.half_edges[adjacency.half_edge].edge;
    brep.record_surface_boundary_audit_event(crate::topology::SurfaceBoundaryAuditEvent {
        half_edge: adjacency.half_edge,
        edge: Some(edge),
        face: Some(adjacency.face),
        twin_face: Some(adjacency.adjacent_face),
        outcome: crate::topology::SurfaceBoundaryAuditOutcome::SharedTrimDeclined,
        boundary: None,
        decline_reason: None,
        shared_trim_failure: Some(
            crate::topology::SurfaceBoundarySharedTrimFailureReason::CurveSurfaceMismatch,
        ),
        shared_trim_validation: Some(crate::topology::SurfaceBoundarySharedTrimValidation {
            max_curve_surface_deviation: 2.0e-4,
            max_side_agreement_deviation: 1.0e-5,
            max_endpoint_gap: 1.0e-6,
        }),
        tolerance: Some(1.0e-6),
        provenance: "unit_test".to_string(),
    });

    let graph = BrepRepairKnowledgeGraph::build(&brep, scope);

    let audit_fact = graph
        .facts
        .iter()
        .find(|fact| {
            fact.header.kind == BrepRepairFactKind::SurfaceBoundary
                && fact.authority.as_ref().is_some_and(|authority| {
                    authority.reason
                        == crate::topology::BrepRepairAuthorityReason::SurfaceBoundarySharedTrimCurveSurfaceMismatch
                })
        })
        .expect("audit event promoted");
    assert_eq!(
        audit_fact.header.state,
        crate::topology::BrepRepairFactState::Contradictory
    );
    assert_eq!(
        audit_fact
            .header
            .tolerance
            .expect("numeric evidence")
            .residual,
        2.0e-4
    );
    assert!(
        graph
            .diagnostics_for_handle(BrepRepairTopologyHandle::HalfEdge(adjacency.half_edge))
            .iter()
            .any(
                |diagnostic| diagnostic.code.0 == "brep.surface.boundary_shared_trim_curve_surface"
            )
    );
}
