//! Shared BRep repair knowledge graph builder and queries.

use std::collections::{BTreeMap, BTreeSet};

use crate::topology::brep_repair_knowledge_diagnostics::build_diagnostic_records;
use crate::topology::brep_repair_knowledge_fact_build::build_fact_records;
use crate::topology::{
    BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION, BRep, BrepKnowledgeFingerprint, BrepRepairDiagnostic,
    BrepRepairFact, BrepRepairFactId, BrepRepairFactKind, BrepRepairKnowledgeGraph,
    BrepRepairKnowledgeScope, BrepRepairScopeExpansion, BrepRepairScopeSeed,
    BrepRepairTopologyHandle, FaceHandle, HalfEdgeHandle, LocalToleranceSource, VertexHandle,
    collect_brep_adjacent_faces, collect_brep_endpoint_memberships, collect_brep_face_loops,
    collect_brep_loop_adjacencies, collect_brep_trim_endpoint_agreements_with_tolerance,
};

impl BrepRepairKnowledgeGraph {
    /// Build a scoped read-only graph from shared BRep repair collectors.
    pub fn build(brep: &BRep, scope: BrepRepairKnowledgeScope) -> Self {
        Self::build_with_local_tolerance(
            brep,
            scope,
            brep.tolerance.linear(),
            LocalToleranceSource::Model,
            9,
        )
    }

    pub(crate) fn build_with_local_tolerance(
        brep: &BRep,
        mut scope: BrepRepairKnowledgeScope,
        tolerance: f64,
        tolerance_source: LocalToleranceSource,
        samples: usize,
    ) -> Self {
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
        let trim_endpoint_agreements = collect_brep_trim_endpoint_agreements_with_tolerance(
            brep,
            &loop_adjacencies,
            tolerance,
        );
        scope.expanded = expanded_handles(&loops, &loop_adjacencies);

        let facts = build_fact_records(
            brep,
            &loops,
            &loop_adjacencies,
            &endpoint_memberships,
            &trim_endpoint_agreements,
            tolerance,
            tolerance_source,
            samples,
        );
        let diagnostics = build_diagnostic_records(&facts);

        Self {
            schema: BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION,
            fingerprint: BrepKnowledgeFingerprint::from_brep(brep),
            scope,
            loops,
            loop_adjacencies,
            endpoint_memberships,
            trim_endpoint_agreements,
            facts,
            diagnostics,
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

    /// Trim endpoint agreement for a half-edge on a specific owning face, if collected.
    pub fn trim_endpoint_agreement_for_face(
        &self,
        half_edge: HalfEdgeHandle,
        face: FaceHandle,
    ) -> Option<&crate::topology::BrepTrimEndpointAgreementFacts> {
        self.loop_adjacencies
            .iter()
            .find(|adjacency| {
                (adjacency.half_edge == half_edge && adjacency.face == face)
                    || (adjacency.adjacent_half_edge == half_edge
                        && adjacency.adjacent_face == face)
            })
            .and_then(|_| self.trim_endpoint_agreement(half_edge))
    }

    /// Surface-boundary fact for a scoped face/half-edge use.
    pub fn surface_boundary_fact(
        &self,
        face: FaceHandle,
        half_edge: HalfEdgeHandle,
    ) -> Option<&BrepRepairFact> {
        self.facts.iter().find(|fact| {
            fact.header.kind == BrepRepairFactKind::SurfaceBoundary
                && fact.handles.contains(&BrepRepairTopologyHandle::Face(face))
                && fact
                    .handles
                    .contains(&BrepRepairTopologyHandle::HalfEdge(half_edge))
        })
    }

    /// Edge/surface candidate fact for a half-edge use.
    pub fn edge_surface_candidate(&self, half_edge: HalfEdgeHandle) -> Option<&BrepRepairFact> {
        self.facts.iter().find(|fact| {
            fact.header.kind == BrepRepairFactKind::EdgeSurfaceCandidate
                && fact
                    .handles
                    .contains(&BrepRepairTopologyHandle::HalfEdge(half_edge))
        })
    }

    /// Authority facts for a half-edge side.
    pub fn authority_for_edge_side(&self, half_edge: HalfEdgeHandle) -> Vec<&BrepRepairFact> {
        self.facts
            .iter()
            .filter(|fact| {
                fact.header.kind == BrepRepairFactKind::Authority
                    && fact
                        .handles
                        .contains(&BrepRepairTopologyHandle::HalfEdge(half_edge))
            })
            .collect()
    }

    /// Authority facts for an edge use on a specific face.
    pub fn authority_for_edge_face(
        &self,
        edge: crate::topology::EdgeHandle,
        face: FaceHandle,
    ) -> Vec<&BrepRepairFact> {
        let side_half_edges = self
            .loop_adjacencies
            .iter()
            .flat_map(|adjacency| {
                [
                    (adjacency.half_edge, adjacency.face),
                    (adjacency.adjacent_half_edge, adjacency.adjacent_face),
                ]
            })
            .filter_map(|(half_edge, side_face)| (side_face == face).then_some(half_edge))
            .collect::<Vec<_>>();
        self.facts
            .iter()
            .filter(|fact| {
                fact.header.kind == BrepRepairFactKind::Authority
                    && fact.handles.contains(&BrepRepairTopologyHandle::Edge(edge))
                    && side_half_edges.iter().any(|half_edge| {
                        fact.handles
                            .contains(&BrepRepairTopologyHandle::HalfEdge(*half_edge))
                    })
            })
            .collect()
    }

    /// Diagnostics that reference a topology handle.
    pub fn diagnostics_for_handle(
        &self,
        handle: BrepRepairTopologyHandle,
    ) -> Vec<&BrepRepairDiagnostic> {
        self.diagnostics
            .iter()
            .filter(|diagnostic| diagnostic.handles.contains(&handle))
            .collect()
    }

    /// Shared fact ids that reference a topology handle and belong to a fact family.
    pub fn fact_ids_for_handle(
        &self,
        handle: BrepRepairTopologyHandle,
        kind: BrepRepairFactKind,
    ) -> Vec<BrepRepairFactId> {
        self.facts
            .iter()
            .filter(|fact| fact.header.kind == kind && fact.handles.contains(&handle))
            .map(|fact| fact.header.id)
            .collect()
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
}
