//! Shared BRep repair knowledge fact construction.

use crate::topology::brep_degenerate_edges::analyze_brep_degenerate_edges_raw;
use crate::topology::brep_repair_knowledge_authority::{
    degenerate_edge_authority_fact, surface_boundary_authority_fact, trim_authority_facts,
};
use crate::topology::brep_repair_knowledge_fact_record::{
    fact, surface_boundary_label, trim_fact_state,
};
use crate::topology::{
    BRep, BrepDegenerateEdgeOptions, BrepRepairFact, BrepRepairFactKind, BrepRepairFactState,
    BrepRepairTopologyHandle, EdgeHandle, LocalToleranceEvidence, LocalToleranceSource,
};
use std::collections::BTreeSet;

pub(super) struct BrepRepairFactBuildInputs<'a> {
    pub loops: &'a [crate::topology::BrepLoopFacts],
    pub adjacencies: &'a [crate::topology::BrepLoopAdjacencyFacts],
    pub endpoints: &'a [crate::topology::BrepEndpointMembershipFacts],
    pub trim_agreements: &'a [crate::topology::BrepTrimEndpointAgreementFacts],
    pub tolerance: f64,
    pub tolerance_source: LocalToleranceSource,
    pub samples: usize,
}

pub(super) fn build_fact_records(
    brep: &BRep,
    inputs: BrepRepairFactBuildInputs<'_>,
) -> Vec<BrepRepairFact> {
    let BrepRepairFactBuildInputs {
        loops,
        adjacencies,
        endpoints,
        trim_agreements,
        tolerance,
        tolerance_source,
        samples,
    } = inputs;
    let mut facts = Vec::new();
    push_loop_facts(&mut facts, loops);
    push_adjacency_facts(&mut facts, brep, adjacencies);
    push_endpoint_facts(&mut facts, endpoints);
    push_trim_facts(
        &mut facts,
        brep,
        trim_agreements,
        tolerance,
        tolerance_source,
        samples,
    );
    push_degenerate_edge_authority_facts(
        &mut facts,
        brep,
        loops,
        adjacencies,
        tolerance,
        tolerance_source,
        samples,
    );
    facts
}

fn push_loop_facts(facts: &mut Vec<BrepRepairFact>, loops: &[crate::topology::BrepLoopFacts]) {
    for loop_facts in loops {
        facts.push(fact(
            facts.len(),
            BrepRepairFactKind::Loop,
            BrepRepairFactState::Observed,
            vec![
                BrepRepairTopologyHandle::Face(loop_facts.face),
                BrepRepairTopologyHandle::Loop(loop_facts.loop_handle),
            ],
            Some(loop_facts.diagnostic_label()),
            None,
        ));
    }
}

fn push_adjacency_facts(
    facts: &mut Vec<BrepRepairFact>,
    brep: &BRep,
    adjacencies: &[crate::topology::BrepLoopAdjacencyFacts],
) {
    for adjacency in adjacencies {
        facts.push(fact(
            facts.len(),
            BrepRepairFactKind::LoopAdjacency,
            BrepRepairFactState::Derived,
            vec![
                BrepRepairTopologyHandle::HalfEdge(adjacency.half_edge),
                BrepRepairTopologyHandle::HalfEdge(adjacency.adjacent_half_edge),
            ],
            Some(adjacency.diagnostic_label()),
            None,
        ));
        push_surface_boundary_side(facts, brep, adjacency.face, adjacency.half_edge);
        push_surface_boundary_side(
            facts,
            brep,
            adjacency.adjacent_face,
            adjacency.adjacent_half_edge,
        );
    }
}

fn push_surface_boundary_side(
    facts: &mut Vec<BrepRepairFact>,
    brep: &BRep,
    face: crate::topology::FaceHandle,
    half_edge: crate::topology::HalfEdgeHandle,
) {
    facts.push(fact(
        facts.len(),
        BrepRepairFactKind::SurfaceBoundary,
        BrepRepairFactState::Derived,
        vec![
            BrepRepairTopologyHandle::Face(face),
            BrepRepairTopologyHandle::HalfEdge(half_edge),
        ],
        Some(surface_boundary_label(face, half_edge)),
        None,
    ));
    facts.push(surface_boundary_authority_fact(
        facts.len(),
        brep,
        face,
        half_edge,
    ));
}

fn push_endpoint_facts(
    facts: &mut Vec<BrepRepairFact>,
    endpoints: &[crate::topology::BrepEndpointMembershipFacts],
) {
    for endpoint in endpoints {
        facts.push(fact(
            facts.len(),
            BrepRepairFactKind::EndpointMembership,
            BrepRepairFactState::Derived,
            vec![BrepRepairTopologyHandle::Vertex(endpoint.vertex)],
            Some(endpoint.diagnostic_label()),
            None,
        ));
    }
}

fn push_trim_facts(
    facts: &mut Vec<BrepRepairFact>,
    brep: &BRep,
    trim_agreements: &[crate::topology::BrepTrimEndpointAgreementFacts],
    tolerance: f64,
    tolerance_source: LocalToleranceSource,
    samples: usize,
) {
    for agreement in trim_agreements {
        let state = trim_fact_state(agreement);
        let evidence =
            tolerance_evidence(tolerance, tolerance_source, agreement.selected_disagreement);
        facts.push(fact(
            facts.len(),
            BrepRepairFactKind::TrimEndpointAgreement,
            state,
            vec![
                BrepRepairTopologyHandle::HalfEdge(agreement.half_edge),
                BrepRepairTopologyHandle::HalfEdge(agreement.adjacent_half_edge),
            ],
            Some(agreement.diagnostic_label()),
            evidence,
        ));
        facts.extend(trim_authority_facts(
            facts.len(),
            brep,
            agreement,
            state,
            tolerance,
            tolerance_source,
            samples,
        ));
        facts.push(fact(
            facts.len(),
            BrepRepairFactKind::EdgeSurfaceCandidate,
            state,
            vec![
                BrepRepairTopologyHandle::HalfEdge(agreement.half_edge),
                BrepRepairTopologyHandle::HalfEdge(agreement.adjacent_half_edge),
            ],
            Some(format!(
                "edge_surface_candidate {}",
                agreement.diagnostic_label()
            )),
            evidence,
        ));
    }
}

fn tolerance_evidence(
    tolerance: f64,
    source: LocalToleranceSource,
    residual: f64,
) -> Option<LocalToleranceEvidence> {
    residual.is_finite().then_some(LocalToleranceEvidence {
        tolerance,
        residual,
        source,
    })
}

fn push_degenerate_edge_authority_facts(
    facts: &mut Vec<BrepRepairFact>,
    brep: &BRep,
    loops: &[crate::topology::BrepLoopFacts],
    adjacencies: &[crate::topology::BrepLoopAdjacencyFacts],
    tolerance: f64,
    tolerance_source: LocalToleranceSource,
    samples: usize,
) {
    let scoped_edges = scoped_edges(brep, loops, adjacencies);
    if scoped_edges.is_empty() {
        return;
    }
    let report =
        analyze_brep_degenerate_edges_raw(brep, BrepDegenerateEdgeOptions { tolerance, samples });
    for finding in report
        .findings
        .iter()
        .filter(|finding| scoped_edges.contains(&finding.edge))
    {
        facts.push(degenerate_edge_authority_fact(
            facts.len(),
            finding,
            tolerance,
            tolerance_source,
        ));
    }
}

fn scoped_edges(
    brep: &BRep,
    loops: &[crate::topology::BrepLoopFacts],
    adjacencies: &[crate::topology::BrepLoopAdjacencyFacts],
) -> BTreeSet<EdgeHandle> {
    let mut edges = BTreeSet::new();
    for loop_facts in loops {
        for half_edge in loop_facts.boundary_half_edges() {
            if let Some(data) = brep.half_edges.get(*half_edge) {
                edges.insert(data.edge);
            }
        }
    }
    for adjacency in adjacencies {
        for half_edge in [adjacency.half_edge, adjacency.adjacent_half_edge] {
            if let Some(data) = brep.half_edges.get(half_edge) {
                edges.insert(data.edge);
            }
        }
    }
    edges
}
