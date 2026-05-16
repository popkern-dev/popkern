//! Shared BRep repair knowledge graph builder and queries.

use std::collections::{BTreeMap, BTreeSet};

use crate::topology::{
    BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION, BRep, BrepKnowledgeFingerprint, BrepRepairFact,
    BrepRepairFactDependencies, BrepRepairFactHeader, BrepRepairFactId, BrepRepairFactKind,
    BrepRepairFactState, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairTopologyHandle, FaceHandle,
    HalfEdgeHandle, LocalToleranceEvidence, LocalToleranceSource, VertexHandle,
    collect_brep_adjacent_faces, collect_brep_endpoint_memberships, collect_brep_face_loops,
    collect_brep_loop_adjacencies, collect_brep_trim_endpoint_agreements,
};

impl BrepRepairKnowledgeGraph {
    /// Build a scoped read-only graph from shared BRep repair collectors.
    pub fn build(brep: &BRep, mut scope: BrepRepairKnowledgeScope) -> Self {
        let faces = scoped_faces(brep, &scope);
        let mut loops_by_face = BTreeMap::new();
        for face in faces {
            loops_by_face.insert(face, collect_brep_face_loops(brep, face));
        }
        if matches!(scope.expansion, BrepRepairScopeExpansion::AdjacentFaces) {
            for face in collect_brep_adjacent_faces(brep, &loops_by_face) {
                loops_by_face
                    .entry(face)
                    .or_insert_with(|| collect_brep_face_loops(brep, face));
            }
        }

        let loops = loops_by_face
            .values()
            .flat_map(|face_loops| face_loops.iter().cloned())
            .collect::<Vec<_>>();
        let loop_adjacencies = collect_brep_loop_adjacencies(brep, &loops_by_face);
        let endpoint_vertices = scoped_endpoint_vertices(&loops);
        let endpoint_memberships =
            collect_brep_endpoint_memberships(&loops_by_face, endpoint_vertices);
        let trim_endpoint_agreements =
            collect_brep_trim_endpoint_agreements(brep, &loop_adjacencies);
        scope.expanded = expanded_handles(&loops, &loop_adjacencies);

        let facts = build_fact_records(
            &loops,
            &loop_adjacencies,
            &endpoint_memberships,
            &trim_endpoint_agreements,
            brep.tolerance.linear(),
        );

        Self {
            schema: BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION,
            fingerprint: BrepKnowledgeFingerprint::from_brep(brep),
            scope,
            loops,
            loop_adjacencies,
            endpoint_memberships,
            trim_endpoint_agreements,
            facts,
            diagnostics: Vec::new(),
        }
    }

    /// Whether this graph still matches the current coarse BRep topology.
    pub fn validate_fresh(&self, brep: &BRep) -> bool {
        self.fingerprint == BrepKnowledgeFingerprint::from_brep(brep)
    }

    /// Loop facts for a face.
    pub fn face_loops(
        &self,
        face: FaceHandle,
    ) -> impl Iterator<Item = &crate::topology::BrepLoopFacts> {
        self.loops.iter().filter(move |facts| facts.face == face)
    }

    /// Endpoint membership facts for a vertex.
    pub fn endpoint_membership(
        &self,
        vertex: VertexHandle,
    ) -> Option<&crate::topology::BrepEndpointMembershipFacts> {
        self.endpoint_memberships
            .iter()
            .find(|membership| membership.vertex == vertex)
    }

    /// Trim endpoint agreement for a half-edge, if collected.
    pub fn trim_endpoint_agreement(
        &self,
        half_edge: HalfEdgeHandle,
    ) -> Option<&crate::topology::BrepTrimEndpointAgreementFacts> {
        self.trim_endpoint_agreements.iter().find(|agreement| {
            agreement.half_edge == half_edge || agreement.adjacent_half_edge == half_edge
        })
    }
}

fn scoped_faces(brep: &BRep, scope: &BrepRepairKnowledgeScope) -> Vec<FaceHandle> {
    let mut faces = BTreeSet::new();
    for seed in &scope.seeds {
        match *seed {
            BrepRepairScopeSeed::Face(face) => {
                if brep.faces.contains_key(face) {
                    faces.insert(face);
                }
            }
            BrepRepairScopeSeed::Loop(loop_handle) => {
                if let Some(loop_data) = brep.loops.get(loop_handle) {
                    faces.insert(loop_data.face);
                }
            }
            BrepRepairScopeSeed::HalfEdge(half_edge) => {
                add_half_edge_face(brep, half_edge, &mut faces);
            }
            BrepRepairScopeSeed::Edge(edge) => {
                if let Some(edge_data) = brep.edges.get(edge) {
                    add_half_edge_face(brep, edge_data.half_edge, &mut faces);
                    if let Some(half_edge) = brep.half_edges.get(edge_data.half_edge) {
                        add_half_edge_face(brep, half_edge.twin, &mut faces);
                    }
                }
            }
            BrepRepairScopeSeed::Vertex(vertex) => {
                for (half_edge, half_edge_data) in &brep.half_edges {
                    if half_edge_data.origin == vertex {
                        add_half_edge_face(brep, half_edge, &mut faces);
                    }
                }
            }
        }
    }
    faces.into_iter().collect()
}

fn add_half_edge_face(brep: &BRep, half_edge: HalfEdgeHandle, faces: &mut BTreeSet<FaceHandle>) {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return;
    };
    let Some(loop_data) = brep.loops.get(half_edge_data.loop_) else {
        return;
    };
    faces.insert(loop_data.face);
}

fn scoped_endpoint_vertices(loops: &[crate::topology::BrepLoopFacts]) -> Vec<VertexHandle> {
    let mut vertices = BTreeSet::new();
    for loop_facts in loops {
        vertices.extend(loop_facts.boundary_vertices().iter().copied());
    }
    vertices.into_iter().collect()
}

fn expanded_handles(
    loops: &[crate::topology::BrepLoopFacts],
    adjacencies: &[crate::topology::BrepLoopAdjacencyFacts],
) -> crate::topology::BrepRepairScopeHandles {
    let mut handles = crate::topology::BrepRepairScopeHandles::default();
    let mut vertices = BTreeSet::new();
    let mut half_edges = BTreeSet::new();
    let mut loops_seen = BTreeSet::new();
    let mut faces = BTreeSet::new();
    for loop_facts in loops {
        faces.insert(loop_facts.face);
        loops_seen.insert(loop_facts.loop_handle);
        vertices.extend(loop_facts.boundary_vertices().iter().copied());
        half_edges.extend(loop_facts.boundary_half_edges().iter().copied());
    }
    for adjacency in adjacencies {
        faces.insert(adjacency.adjacent_face);
        loops_seen.insert(adjacency.adjacent_loop);
        half_edges.insert(adjacency.adjacent_half_edge);
    }
    handles.vertices = vertices.into_iter().collect();
    handles.half_edges = half_edges.into_iter().collect();
    handles.loops = loops_seen.into_iter().collect();
    handles.faces = faces.into_iter().collect();
    handles
}

fn build_fact_records(
    loops: &[crate::topology::BrepLoopFacts],
    adjacencies: &[crate::topology::BrepLoopAdjacencyFacts],
    endpoints: &[crate::topology::BrepEndpointMembershipFacts],
    trim_agreements: &[crate::topology::BrepTrimEndpointAgreementFacts],
    tolerance: f64,
) -> Vec<BrepRepairFact> {
    let mut facts = Vec::new();
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
    }
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
    for agreement in trim_agreements {
        facts.push(fact(
            facts.len(),
            BrepRepairFactKind::TrimEndpointAgreement,
            trim_fact_state(agreement),
            vec![
                BrepRepairTopologyHandle::HalfEdge(agreement.half_edge),
                BrepRepairTopologyHandle::HalfEdge(agreement.adjacent_half_edge),
            ],
            Some(agreement.diagnostic_label()),
            agreement
                .selected_disagreement
                .is_finite()
                .then_some(LocalToleranceEvidence {
                    tolerance,
                    residual: agreement.selected_disagreement,
                    source: LocalToleranceSource::Model,
                }),
        ));
    }
    facts
}

fn fact(
    id: usize,
    kind: BrepRepairFactKind,
    state: BrepRepairFactState,
    handles: Vec<BrepRepairTopologyHandle>,
    legacy_label: Option<String>,
    tolerance: Option<LocalToleranceEvidence>,
) -> BrepRepairFact {
    BrepRepairFact {
        header: BrepRepairFactHeader {
            id: BrepRepairFactId(id),
            kind,
            state,
            dependencies: BrepRepairFactDependencies::default(),
            provenance: Vec::new(),
            tolerance,
        },
        handles,
        legacy_label,
    }
}

fn trim_fact_state(
    agreement: &crate::topology::BrepTrimEndpointAgreementFacts,
) -> BrepRepairFactState {
    use crate::topology::BrepTrimEndpointAgreementStatus as Status;
    match agreement.status {
        Status::AcceptedDirect | Status::AcceptedReconciled => BrepRepairFactState::Derived,
        Status::SideDisagreement => BrepRepairFactState::Contradictory,
        Status::MissingSurface
        | Status::MissingPcurve
        | Status::MissingTopology
        | Status::EvaluationFailed => BrepRepairFactState::Unavailable,
    }
}

#[cfg(test)]
mod tests {
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
                .facts
                .iter()
                .any(|fact| fact.header.kind == BrepRepairFactKind::TrimEndpointAgreement)
        );
    }
}
