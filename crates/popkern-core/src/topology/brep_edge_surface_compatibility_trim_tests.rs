use super::*;
use crate::topology::{
    BrepEdgeSurfaceVertexCandidateOutcome, BrepEdgeSurfaceVertexCandidateSource,
    BrepRepairAuthorityKind, BrepRepairAuthorityReason, BrepRepairFactKind, BrepRepairFactState,
    BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeSeed,
    SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason,
};

#[test]
fn vertex_star_accepts_noncanonical_trim_endpoint_candidates() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    let he_b = brep.half_edges[he_a].twin;
    let va = brep.half_edges[he_a].origin;
    let vb = brep.half_edges[he_b].origin;
    let start = Point3::new(0.25, 0.0, 0.0);
    let end = Point3::new(0.25, 1.0, 0.0);
    brep.vertices[va].position = start;
    brep.vertices[vb].position = end;
    brep.half_edges[he_a].pcurve = Some(Pcurve::linear((0.25, 0.0), (0.25, 1.0), 0.0, 1.0));
    brep.half_edges[he_b].pcurve = Some(Pcurve::linear((0.25, 1.0), (0.25, 0.0), 0.0, 1.0));
    brep.edges[edge].curve = CurveOnEdge::line_segment(
        start + Vec3::new(0.0, 0.0, 0.1),
        end + Vec3::new(0.0, 0.0, 0.1),
        1e-12,
    );

    let report = match repair_brep_edge_surface_compatibility_transactionally(
        &brep,
        no_validation_options(),
    ) {
        Ok(result) => result.report,
        Err(failure) => failure.report,
    };

    assert_eq!(report.vertex_stars.inspected_vertices, 2);
    assert_eq!(report.vertex_stars.coherent_vertices, 2);
    assert!(
        report
            .vertex_stars
            .decisions
            .iter()
            .all(|decision| decision.max_anchor_motion <= no_validation_options().tolerance)
    );
    assert!(report.vertex_stars.decisions.iter().all(|decision| {
        decision.candidate_diagnostics.iter().any(|diagnostic| {
            diagnostic.source == BrepEdgeSurfaceVertexCandidateSource::TrimEndpoint
                && diagnostic.outcome == BrepEdgeSurfaceVertexCandidateOutcome::Accepted
        })
    }));
}

#[test]
fn migration_matrix_edge_surface_candidate_contract() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    let he_b = brep.half_edges[he_a].twin;
    let va = brep.half_edges[he_a].origin;
    let vb = brep.half_edges[he_b].origin;
    let start = Point3::new(0.25, 0.0, 0.0);
    let end = Point3::new(0.25, 1.0, 0.0);
    brep.vertices[va].position = start;
    brep.vertices[vb].position = end;
    brep.half_edges[he_a].pcurve = Some(Pcurve::linear((0.25, 0.0), (0.25, 1.0), 0.0, 1.0));
    brep.half_edges[he_b].pcurve = Some(Pcurve::linear((0.25, 1.0), (0.25, 0.0), 0.0, 1.0));
    brep.edges[edge].curve = CurveOnEdge::line_segment(
        start + Vec3::new(0.0, 0.0, 0.1),
        end + Vec3::new(0.0, 0.0, 0.1),
        1e-12,
    );

    let report = match repair_brep_edge_surface_compatibility_transactionally(
        &brep,
        no_validation_options(),
    ) {
        Ok(result) => result.report,
        Err(failure) => failure.report,
    };

    assert_eq!(report.vertex_stars.inspected_vertices, 2);
    assert_eq!(report.vertex_stars.coherent_vertices, 2);
    assert_eq!(report.vertex_stars.decisions.len(), 2);
    assert!(
        !report.shared_fact_ids.is_empty(),
        "repair report must expose aggregate shared repair fact ids"
    );
    let knowledge =
        crate::topology::brep_edge_surface_compatibility_shared_facts::edge_surface_repair_knowledge(
            &brep,
            &report.before,
        );
    assert_shared_ids_are_from_graph(&report, &knowledge);
    assert!(knowledge.facts.iter().any(|fact| {
        fact.header.kind == BrepRepairFactKind::Authority
            && fact.authority.as_ref().is_some_and(|authority| {
                authority.kind == BrepRepairAuthorityKind::AuthoritativeEdgeCurve
                    || authority.kind == BrepRepairAuthorityKind::AuthoritativeTrimEndpoint
                    || authority.kind == BrepRepairAuthorityKind::AuthoritativeSurfaceBoundary
            })
    }));
    for decision in &report.vertex_stars.decisions {
        assert!(decision.max_anchor_motion <= no_validation_options().tolerance);
        assert!(
            !decision.shared_fact_ids.is_empty(),
            "vertex-star decision must expose shared repair fact ids"
        );
        assert!(decision.candidate_diagnostics.iter().any(|diagnostic| {
            diagnostic.source == BrepEdgeSurfaceVertexCandidateSource::TrimEndpoint
                && diagnostic.outcome == BrepEdgeSurfaceVertexCandidateOutcome::Accepted
                && !diagnostic.shared_fact_ids.is_empty()
        }));
        assert!(decision.candidate_diagnostics.iter().any(|diagnostic| {
            diagnostic.source == BrepEdgeSurfaceVertexCandidateSource::SurfaceBoundary
                && !diagnostic.shared_fact_ids.is_empty()
        }));
        assert!(decision.candidate_diagnostics.iter().any(|diagnostic| {
            diagnostic.source == BrepEdgeSurfaceVertexCandidateSource::SideAgreement
                && !diagnostic.shared_fact_ids.is_empty()
        }));
        assert!(decision.candidate_diagnostics.iter().any(|diagnostic| {
            diagnostic.shared_fact_ids.iter().any(|id| {
                knowledge
                    .facts
                    .get(id.0)
                    .is_some_and(|fact| fact.header.kind == BrepRepairFactKind::Authority)
            })
        }));
    }
}

#[test]
fn authority_keeps_missing_trim_and_edge_curve_sources_distinct() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    brep.half_edges[he_a].pcurve = None;

    let knowledge = BrepRepairKnowledgeGraph::build(
        &brep,
        BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Edge(edge)]),
    );
    let authority = knowledge.authority_for_edge_side(he_a);

    assert!(authority.iter().any(|fact| {
        fact.header.state == BrepRepairFactState::Unavailable
            && fact.authority.as_ref().is_some_and(|authority| {
                authority.kind == BrepRepairAuthorityKind::NoAuthoritativeSource
                    && authority.reason == BrepRepairAuthorityReason::MissingPcurve
            })
    }));
    assert!(authority.iter().any(|fact| {
        fact.authority.as_ref().is_some_and(|authority| {
            authority.kind == BrepRepairAuthorityKind::AuthoritativeEdgeCurve
                && authority.reason == BrepRepairAuthorityReason::EdgeCanonical3D
        })
    }));
}

#[test]
fn authority_marks_declined_surface_boundary_evidence_as_contradictory() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    assert!(brep.set_surface_boundary_decline(
        he_a,
        SurfaceBoundaryHintDecline {
            reason: SurfaceBoundaryHintDeclineReason::AuthoredHintBoundaryMismatch,
            best_boundary_deviation: 0.25,
            best_boundary: Some(crate::geometry::edge_triple::BoundarySide::UMin),
            endpoint_gap: 0.1,
        },
    ));

    let knowledge = BrepRepairKnowledgeGraph::build(
        &brep,
        BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Edge(edge)]),
    );
    let authority = knowledge.authority_for_edge_side(he_a);

    assert!(authority.iter().any(|fact| {
        fact.header.state == BrepRepairFactState::Contradictory
            && fact.authority.as_ref().is_some_and(|authority| {
                authority.kind == BrepRepairAuthorityKind::SideLocalSurfaceMismatch
                    && authority.reason
                        == BrepRepairAuthorityReason::SurfaceBoundaryTrimDisagreement
            })
    }));
}

#[test]
fn authority_marks_existing_pcurve_that_disagrees_with_edge_curve_as_stale() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    brep.half_edges[he_a].pcurve = Some(Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0));

    let knowledge = BrepRepairKnowledgeGraph::build(
        &brep,
        BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Edge(edge)]),
    );
    let authority = knowledge.authority_for_edge_side(he_a);

    assert!(authority.iter().any(|fact| {
        fact.header.state == BrepRepairFactState::Stale
            && fact.authority.as_ref().is_some_and(|authority| {
                authority.kind == BrepRepairAuthorityKind::StalePcurve
                    && authority.reason == BrepRepairAuthorityReason::StalePcurve
            })
    }));
}

#[test]
fn authority_records_generated_and_untrusted_trim_source_forms() {
    let mut brep = gregory_plane_boundary_fixture();
    let edge = brep.edges.keys().next().expect("fixture has one edge");
    let he_a = brep.edges[edge].half_edge;
    let he_b = brep.half_edges[he_a].twin;
    let nurbs_pcurve = brep.half_edges[he_b]
        .pcurve
        .as_ref()
        .expect("fixture has twin pcurve")
        .to_nurbs();
    brep.half_edges[he_b].pcurve = Some(Pcurve::from_nurbs(nurbs_pcurve));

    let knowledge = BrepRepairKnowledgeGraph::build(
        &brep,
        BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Edge(edge)]),
    );

    assert!(knowledge.authority_for_edge_side(he_a).iter().any(|fact| {
        fact.authority.as_ref().is_some_and(|authority| {
            authority.kind == BrepRepairAuthorityKind::GeneratedConstructionRecord
                && authority.reason == BrepRepairAuthorityReason::GeneratedAnalyticPcurve
        })
    }));
    assert!(knowledge.authority_for_edge_side(he_b).iter().any(|fact| {
        fact.authority.as_ref().is_some_and(|authority| {
            authority.kind == BrepRepairAuthorityKind::ImportedUntrustedTrim
                && authority.reason == BrepRepairAuthorityReason::NurbsTrimRequiresValidation
        })
    }));
}

fn assert_shared_ids_are_from_graph(
    report: &BrepEdgeSurfaceCompatibilityRepairReport,
    knowledge: &BrepRepairKnowledgeGraph,
) {
    for id in &report.shared_fact_ids {
        assert!(
            id.0 < knowledge.facts.len(),
            "report shared fact id {:?} must refer to the repair knowledge graph",
            id
        );
    }
    for decision in &report.decisions {
        for id in &decision.shared_fact_ids {
            assert!(
                id.0 < knowledge.facts.len(),
                "edge decision shared fact id {:?} must refer to the repair knowledge graph",
                id
            );
        }
    }
    for decision in &report.vertex_stars.decisions {
        for id in &decision.shared_fact_ids {
            assert!(
                id.0 < knowledge.facts.len(),
                "vertex-star shared fact id {:?} must refer to the repair knowledge graph",
                id
            );
        }
        for diagnostic in &decision.candidate_diagnostics {
            for id in &diagnostic.shared_fact_ids {
                assert!(
                    id.0 < knowledge.facts.len(),
                    "candidate diagnostic shared fact id {:?} must refer to the repair knowledge graph",
                    id
                );
            }
        }
    }
}
