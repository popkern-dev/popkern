//! Read-only topology/geometry knowledge for repeated-apex side-support planning.

use std::collections::BTreeSet;

use super::repeated_apex_trim_subcomplex::{TrimRewriteIntentFacts, TrimRewriteSubcomplexPlan};
use popkern_core::topology::{
    BRep, BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeSeed,
    BrepTrimEndpointAgreementFacts as TrimEndpointAgreementFacts, FaceHandle, VertexHandle,
    collect_brep_adjacent_faces,
};
pub(super) use popkern_core::topology::{
    BrepEndpointMembershipFacts as EndpointMembershipFacts,
    BrepLoopAdjacencyFacts as LoopAdjacencyFacts, BrepLoopFacts as LoopFacts,
    BrepLoopKind as LoopKind,
};

#[derive(Clone, Debug)]
pub(super) struct TopologyKnowledgeGraph {
    core: BrepRepairKnowledgeGraph,
    overlay: ApexRepairKnowledgeOverlay,
}

#[derive(Clone, Debug)]
pub(super) struct ApexRepairKnowledgeOverlay {
    context_faces: Vec<FaceHandle>,
    support_families: Vec<SupportFamilyFacts>,
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
        let core = BrepRepairKnowledgeGraph::build(
            brep,
            BrepRepairKnowledgeScope::seeds_only(
                seed_faces
                    .iter()
                    .copied()
                    .map(BrepRepairScopeSeed::Face)
                    .collect(),
            ),
        );
        let overlay = ApexRepairKnowledgeOverlay::from_core(brep, &core, plan, seed_faces);

        Self { core, overlay }
    }

    pub(super) fn validate_fresh(&self, brep: &BRep) -> bool {
        self.core.validate_fresh(brep)
    }

    pub(super) fn support_families(&self) -> impl Iterator<Item = &SupportFamilyFacts> {
        self.overlay.support_families.iter()
    }

    pub(super) fn context_faces(&self) -> &[FaceHandle] {
        &self.overlay.context_faces
    }

    pub(super) fn trim_intents(&self) -> impl Iterator<Item = &TrimRewriteIntentFacts> {
        self.overlay.trim_intents.iter()
    }

    pub(super) fn loop_adjacencies(&self) -> impl Iterator<Item = &LoopAdjacencyFacts> {
        self.core.loop_adjacencies.iter()
    }

    pub(super) fn endpoint_membership(
        &self,
        vertex: VertexHandle,
    ) -> Option<&EndpointMembershipFacts> {
        self.core.endpoint_membership(vertex)
    }

    pub(super) fn diagnostic_label(&self) -> String {
        self.overlay.diagnostic_label(&self.core)
    }
}

impl ApexRepairKnowledgeOverlay {
    fn from_core(
        brep: &BRep,
        core: &BrepRepairKnowledgeGraph,
        plan: &TrimRewriteSubcomplexPlan,
        seed_faces: Vec<FaceHandle>,
    ) -> Self {
        let half_edges_visited = core
            .loops
            .iter()
            .map(|loop_facts| loop_facts.boundary_half_edges().len())
            .sum();
        let vertices_seen = core
            .loops
            .iter()
            .flat_map(|loop_facts| loop_facts.boundary_vertices().iter().copied())
            .collect::<BTreeSet<_>>();
        let mut context_faces = seed_faces.clone();
        let adjacent_faces = collect_brep_adjacent_faces(brep, &loops_by_face(core));
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
                loops: core.face_loops(group.support_face()).cloned().collect(),
            })
            .collect();

        let trim_intents = plan.trim_intents();
        let trim_records_read = trim_intents.len();

        Self {
            context_faces,
            support_families,
            trim_intents,
            counters: KnowledgeGraphCounters {
                seeded_faces: seed_faces.len(),
                loops_visited: core.loops.len(),
                half_edges_visited,
                vertices_visited: vertices_seen.len(),
                trim_records_read,
            },
        }
    }

    fn diagnostic_label(&self, core: &BrepRepairKnowledgeGraph) -> String {
        let trim_endpoint_vertices = self
            .trim_intents
            .iter()
            .flat_map(|trim| [trim.start, trim.end])
            .collect::<BTreeSet<_>>();
        format!(
            "knowledge fingerprint {}; counters seeded_faces={} loops={} half_edges={} vertices={} trim_records={}; support_families={}; loop_adjacencies={}; trim_endpoint_agreements={}; endpoints={}",
            core.fingerprint.diagnostic_label(),
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
            core.loop_adjacencies
                .iter()
                .map(LoopAdjacencyFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            core.trim_endpoint_agreements
                .iter()
                .map(TrimEndpointAgreementFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            core.endpoint_memberships
                .iter()
                .filter(|membership| trim_endpoint_vertices.contains(&membership.vertex))
                .map(EndpointMembershipFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

fn loops_by_face(
    core: &BrepRepairKnowledgeGraph,
) -> std::collections::BTreeMap<FaceHandle, Vec<LoopFacts>> {
    let mut loops_by_face = std::collections::BTreeMap::new();
    for loop_facts in &core.loops {
        loops_by_face
            .entry(loop_facts.face)
            .or_insert_with(Vec::new)
            .push(loop_facts.clone());
    }
    loops_by_face
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
