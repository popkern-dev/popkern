//! Read-only topology/geometry knowledge for repeated-apex side-support planning.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{
    BRep, BrepKnowledgeFingerprint, FaceHandle, LoopHandle, VertexHandle,
};
pub(super) use popkern_core::topology::{
    BrepEndpointLoopMembership as EndpointLoopMembership,
    BrepEndpointMembershipFacts as EndpointMembershipFacts,
    BrepLoopAdjacencyFacts as LoopAdjacencyFacts, BrepLoopFacts as LoopFacts,
    BrepLoopKind as LoopKind,
};

use super::repeated_apex_topology_knowledge_collect::{
    collect_adjacent_faces, collect_face_loops, collect_loop_adjacencies,
};
use super::repeated_apex_topology_knowledge_trim::{
    TrimEndpointAgreementFacts, collect_trim_endpoint_agreements,
};
use super::repeated_apex_trim_subcomplex::{TrimRewriteIntentFacts, TrimRewriteSubcomplexPlan};

#[derive(Clone, Debug)]
pub(super) struct TopologyKnowledgeGraph {
    fingerprint: BrepKnowledgeFingerprint,
    context_faces: Vec<FaceHandle>,
    support_families: Vec<SupportFamilyFacts>,
    loop_adjacencies: Vec<LoopAdjacencyFacts>,
    trim_endpoint_agreements: Vec<TrimEndpointAgreementFacts>,
    endpoint_memberships: Vec<EndpointMembershipFacts>,
    trim_intents: Vec<TrimRewriteIntentFacts>,
    counters: KnowledgeGraphCounters,
}

#[derive(Clone, Debug)]
pub(super) struct KnowledgeGraphCounters {
    seeded_faces: usize,
    loops_visited: usize,
    half_edges_visited: usize,
    vertices_visited: usize,
    trim_records_read: usize,
}

#[derive(Clone, Debug)]
pub(super) struct SupportFamilyFacts {
    support_face: FaceHandle,
    generated_trim_count: usize,
    cylinder_faces: Vec<FaceHandle>,
    loops: Vec<LoopFacts>,
}

impl TopologyKnowledgeGraph {
    pub(super) fn build_local(brep: &BRep, plan: &TrimRewriteSubcomplexPlan) -> Self {
        let seed_faces = plan.seed_faces();

        let mut loops_by_face = BTreeMap::new();
        let mut half_edges_visited = 0usize;
        let mut vertices_seen = BTreeSet::new();
        for face in seed_faces.iter().copied() {
            let loops = collect_face_loops(brep, face, &mut half_edges_visited, &mut vertices_seen);
            loops_by_face.insert(face, loops);
        }
        let mut context_faces = seed_faces.clone();
        let adjacent_faces = collect_adjacent_faces(brep, &loops_by_face);
        context_faces.extend(adjacent_faces);
        context_faces.sort();
        context_faces.dedup();

        let support_families = plan
            .support_groups()
            .iter()
            .map(|group| SupportFamilyFacts {
                support_face: group.support_face(),
                generated_trim_count: group.generated_trim_count(),
                cylinder_faces: group.cylinder_faces().to_vec(),
                loops: loops_by_face
                    .get(&group.support_face())
                    .cloned()
                    .unwrap_or_default(),
            })
            .collect();
        let loop_adjacencies = collect_loop_adjacencies(brep, &loops_by_face);
        let trim_endpoint_agreements = collect_trim_endpoint_agreements(brep, &loop_adjacencies);

        let trim_intents = plan.trim_intents();
        let trim_records_read = trim_intents.len();

        let mut endpoint_vertices = BTreeSet::new();
        for trim in &trim_intents {
            endpoint_vertices.insert(trim.start);
            endpoint_vertices.insert(trim.end);
        }
        let endpoint_memberships = endpoint_vertices
            .into_iter()
            .map(|vertex| EndpointMembershipFacts {
                vertex,
                containing_loops: loops_by_face
                    .values()
                    .flat_map(|loops| loops.iter())
                    .filter_map(|loop_facts| {
                        let position = loop_facts
                            .boundary_vertices
                            .iter()
                            .position(|loop_vertex| *loop_vertex == vertex)?;
                        Some(EndpointLoopMembership {
                            face: loop_facts.face,
                            loop_handle: loop_facts.loop_handle,
                            half_edge: loop_facts.boundary_half_edges[position],
                            kind: loop_facts.kind,
                        })
                    })
                    .collect(),
            })
            .collect();

        Self {
            fingerprint: BrepKnowledgeFingerprint::from_brep(brep),
            context_faces,
            support_families,
            loop_adjacencies,
            trim_endpoint_agreements,
            endpoint_memberships,
            trim_intents,
            counters: KnowledgeGraphCounters {
                seeded_faces: seed_faces.len(),
                loops_visited: loops_by_face.values().map(Vec::len).sum(),
                half_edges_visited,
                vertices_visited: vertices_seen.len(),
                trim_records_read,
            },
        }
    }

    pub(super) fn validate_fresh(&self, brep: &BRep) -> bool {
        self.fingerprint == BrepKnowledgeFingerprint::from_brep(brep)
    }

    pub(super) fn support_families(&self) -> impl Iterator<Item = &SupportFamilyFacts> {
        self.support_families.iter()
    }

    pub(super) fn context_faces(&self) -> &[FaceHandle] {
        &self.context_faces
    }

    pub(super) fn trim_intents(&self) -> impl Iterator<Item = &TrimRewriteIntentFacts> {
        self.trim_intents.iter()
    }

    pub(super) fn loop_adjacencies(&self) -> impl Iterator<Item = &LoopAdjacencyFacts> {
        self.loop_adjacencies.iter()
    }

    pub(super) fn endpoint_membership(
        &self,
        vertex: VertexHandle,
    ) -> Option<&EndpointMembershipFacts> {
        self.endpoint_memberships
            .iter()
            .find(|membership| membership.vertex == vertex)
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "knowledge fingerprint {}; counters seeded_faces={} loops={} half_edges={} vertices={} trim_records={}; support_families={}; loop_adjacencies={}; trim_endpoint_agreements={}; endpoints={}",
            self.fingerprint.diagnostic_label(),
            self.counters.seeded_faces,
            self.counters.loops_visited,
            self.counters.half_edges_visited,
            self.counters.vertices_visited,
            self.counters.trim_records_read,
            self.support_families
                .iter()
                .map(SupportFamilyFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.loop_adjacencies
                .iter()
                .map(LoopAdjacencyFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.trim_endpoint_agreements
                .iter()
                .map(TrimEndpointAgreementFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.endpoint_memberships
                .iter()
                .map(EndpointMembershipFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl SupportFamilyFacts {
    pub(super) fn support_face(&self) -> FaceHandle {
        self.support_face
    }

    pub(super) fn loops(&self) -> &[LoopFacts] {
        &self.loops
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} generated_trims={} cylinders={:?} loops=[{}]",
            self.support_face,
            self.generated_trim_count,
            self.cylinder_faces,
            self.loops
                .iter()
                .map(LoopFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

pub(super) fn collect_loop_facts(
    brep: &BRep,
    face: FaceHandle,
    loop_handle: LoopHandle,
    kind: LoopKind,
    half_edges_visited: &mut usize,
    vertices_seen: &mut BTreeSet<VertexHandle>,
) -> LoopFacts {
    let mut boundary_vertices = Vec::new();
    let mut boundary_half_edges = Vec::new();
    let mut seen_half_edges = BTreeSet::new();
    let mut repeated_half_edge = false;
    let mut repeated_vertex = false;
    let Some(loop_data) = brep.loops.get(loop_handle) else {
        return LoopFacts {
            face,
            loop_handle,
            kind: LoopKind::Unknown,
            boundary_half_edges,
            boundary_vertices,
            repeated_vertex,
            repeated_half_edge,
        };
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    loop {
        if !seen_half_edges.insert(cursor) {
            repeated_half_edge = true;
            break;
        }
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            break;
        };
        *half_edges_visited += 1;
        if boundary_vertices.contains(&half_edge.origin) {
            repeated_vertex = true;
        }
        boundary_half_edges.push(cursor);
        boundary_vertices.push(half_edge.origin);
        vertices_seen.insert(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    LoopFacts {
        face,
        loop_handle,
        kind,
        boundary_half_edges,
        boundary_vertices,
        repeated_vertex,
        repeated_half_edge,
    }
}
