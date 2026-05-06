//! Read-only topology/geometry knowledge for repeated-apex side-support planning.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_trim_subcomplex::{TrimRewriteIntentFacts, TrimRewriteSubcomplexPlan};

#[derive(Clone, Debug)]
pub(super) struct TopologyKnowledgeGraph {
    fingerprint: TopologyFingerprint,
    support_families: Vec<SupportFamilyFacts>,
    endpoint_memberships: Vec<EndpointMembershipFacts>,
    trim_intents: Vec<TrimRewriteIntentFacts>,
    counters: KnowledgeGraphCounters,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct TopologyFingerprint {
    vertices: usize,
    edges: usize,
    half_edges: usize,
    loops: usize,
    faces: usize,
    shells: usize,
}

impl TopologyFingerprint {
    pub(super) fn from_brep(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len(),
            edges: brep.edges.len(),
            half_edges: brep.half_edges.len(),
            loops: brep.loops.len(),
            faces: brep.faces.len(),
            shells: brep.shells.len(),
        }
    }

    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "V={} E={} HE={} L={} F={} S={}",
            self.vertices, self.edges, self.half_edges, self.loops, self.faces, self.shells
        )
    }
}

#[derive(Clone, Debug)]
pub(super) struct KnowledgeGraphCounters {
    pub(super) seeded_faces: usize,
    pub(super) loops_visited: usize,
    pub(super) half_edges_visited: usize,
    pub(super) vertices_visited: usize,
    pub(super) trim_records_read: usize,
}

#[derive(Clone, Debug)]
pub(super) struct SupportFamilyFacts {
    pub(super) support_face: FaceHandle,
    pub(super) generated_trim_count: usize,
    pub(super) cylinder_faces: Vec<FaceHandle>,
    pub(super) loops: Vec<LoopFacts>,
}

#[derive(Clone, Debug)]
pub(super) struct LoopFacts {
    pub(super) face: FaceHandle,
    pub(super) loop_handle: LoopHandle,
    pub(super) kind: LoopKind,
    pub(super) boundary_half_edges: Vec<HalfEdgeHandle>,
    pub(super) boundary_vertices: Vec<VertexHandle>,
    pub(super) repeated_vertex: bool,
    pub(super) repeated_half_edge: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum LoopKind {
    Outer,
    Inner,
    Unknown,
}

#[derive(Clone, Debug)]
pub(super) struct EndpointMembershipFacts {
    pub(super) vertex: VertexHandle,
    pub(super) containing_loops: Vec<EndpointLoopMembership>,
}

#[derive(Clone, Debug)]
pub(super) struct EndpointLoopMembership {
    pub(super) face: FaceHandle,
    pub(super) loop_handle: LoopHandle,
    pub(super) half_edge: HalfEdgeHandle,
    pub(super) kind: LoopKind,
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
            fingerprint: TopologyFingerprint::from_brep(brep),
            support_families,
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
        self.fingerprint == TopologyFingerprint::from_brep(brep)
    }

    pub(super) fn support_families(&self) -> impl Iterator<Item = &SupportFamilyFacts> {
        self.support_families.iter()
    }

    pub(super) fn trim_intents(&self) -> impl Iterator<Item = &TrimRewriteIntentFacts> {
        self.trim_intents.iter()
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
            "knowledge fingerprint {}; counters seeded_faces={} loops={} half_edges={} vertices={} trim_records={}; support_families={}; endpoints={}",
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
            self.endpoint_memberships
                .iter()
                .map(EndpointMembershipFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl SupportFamilyFacts {
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

impl LoopFacts {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}/{:?} vertices={:?} repeated_vertex={} repeated_half_edge={}",
            self.face,
            self.loop_handle,
            self.kind,
            self.boundary_vertices,
            self.repeated_vertex,
            self.repeated_half_edge
        )
    }
}

impl EndpointMembershipFacts {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?} in [{}]",
            self.vertex,
            self.containing_loops
                .iter()
                .map(EndpointLoopMembership::diagnostic_label)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

impl EndpointLoopMembership {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}/{:?}@{:?}",
            self.face, self.loop_handle, self.kind, self.half_edge
        )
    }
}

fn collect_face_loops(
    brep: &BRep,
    face_handle: FaceHandle,
    half_edges_visited: &mut usize,
    vertices_seen: &mut BTreeSet<VertexHandle>,
) -> Vec<LoopFacts> {
    let Some(face) = brep.faces.get(face_handle) else {
        return Vec::new();
    };
    let mut loop_handles = Vec::with_capacity(1 + face.inner_loops.len());
    loop_handles.push((face.outer_loop, LoopKind::Outer));
    loop_handles.extend(
        face.inner_loops
            .iter()
            .copied()
            .map(|loop_handle| (loop_handle, LoopKind::Inner)),
    );
    loop_handles
        .into_iter()
        .map(|(loop_handle, kind)| {
            collect_loop_facts(
                brep,
                face_handle,
                loop_handle,
                kind,
                half_edges_visited,
                vertices_seen,
            )
        })
        .collect()
}

fn collect_loop_facts(
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
