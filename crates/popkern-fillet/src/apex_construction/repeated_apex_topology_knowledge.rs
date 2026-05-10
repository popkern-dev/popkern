//! Read-only topology/geometry knowledge for repeated-apex side-support planning.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_trim_subcomplex::{TrimRewriteIntentFacts, TrimRewriteSubcomplexPlan};

#[derive(Clone, Debug)]
pub(super) struct TopologyKnowledgeGraph {
    fingerprint: TopologyFingerprint,
    support_families: Vec<SupportFamilyFacts>,
    loop_adjacencies: Vec<LoopAdjacencyFacts>,
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

#[derive(Clone, Debug)]
pub(super) struct LoopFacts {
    face: FaceHandle,
    loop_handle: LoopHandle,
    kind: LoopKind,
    boundary_half_edges: Vec<HalfEdgeHandle>,
    boundary_vertices: Vec<VertexHandle>,
    repeated_vertex: bool,
    repeated_half_edge: bool,
}

#[derive(Clone, Debug)]
pub(super) struct LoopAdjacencyFacts {
    face: FaceHandle,
    loop_handle: LoopHandle,
    half_edge: HalfEdgeHandle,
    adjacent_face: FaceHandle,
    adjacent_loop: LoopHandle,
    adjacent_half_edge: HalfEdgeHandle,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum LoopKind {
    Outer,
    Inner,
    Unknown,
}

#[derive(Clone, Debug)]
pub(super) struct EndpointMembershipFacts {
    vertex: VertexHandle,
    containing_loops: Vec<EndpointLoopMembership>,
}

#[derive(Clone, Debug)]
pub(super) struct EndpointLoopMembership {
    face: FaceHandle,
    loop_handle: LoopHandle,
    half_edge: HalfEdgeHandle,
    kind: LoopKind,
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
        let loop_adjacencies = collect_loop_adjacencies(brep, &loops_by_face);

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
            loop_adjacencies,
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
            "knowledge fingerprint {}; counters seeded_faces={} loops={} half_edges={} vertices={} trim_records={}; support_families={}; loop_adjacencies={}; endpoints={}",
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
            self.endpoint_memberships
                .iter()
                .map(EndpointMembershipFacts::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl LoopAdjacencyFacts {
    pub(super) fn face(&self) -> FaceHandle {
        self.face
    }

    pub(super) fn loop_handle(&self) -> LoopHandle {
        self.loop_handle
    }

    pub(super) fn half_edge(&self) -> HalfEdgeHandle {
        self.half_edge
    }

    pub(super) fn adjacent_face(&self) -> FaceHandle {
        self.adjacent_face
    }

    pub(super) fn adjacent_loop(&self) -> LoopHandle {
        self.adjacent_loop
    }

    pub(super) fn adjacent_half_edge(&self) -> HalfEdgeHandle {
        self.adjacent_half_edge
    }

    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}/{:?}@{:?}<->{:?}/{:?}@{:?}",
            self.face,
            self.loop_handle,
            self.half_edge,
            self.adjacent_face,
            self.adjacent_loop,
            self.adjacent_half_edge
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

impl LoopFacts {
    pub(super) fn face(&self) -> FaceHandle {
        self.face
    }

    pub(super) fn loop_handle(&self) -> LoopHandle {
        self.loop_handle
    }

    pub(super) fn kind(&self) -> LoopKind {
        self.kind
    }

    pub(super) fn boundary_half_edges(&self) -> &[HalfEdgeHandle] {
        &self.boundary_half_edges
    }

    pub(super) fn boundary_vertices(&self) -> &[VertexHandle] {
        &self.boundary_vertices
    }

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
    pub(super) fn containing_loops(&self) -> &[EndpointLoopMembership] {
        &self.containing_loops
    }

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
    pub(super) fn face(&self) -> FaceHandle {
        self.face
    }

    pub(super) fn loop_handle(&self) -> LoopHandle {
        self.loop_handle
    }

    pub(super) fn half_edge(&self) -> HalfEdgeHandle {
        self.half_edge
    }

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

fn collect_loop_adjacencies(
    brep: &BRep,
    loops_by_face: &BTreeMap<FaceHandle, Vec<LoopFacts>>,
) -> Vec<LoopAdjacencyFacts> {
    let seeded_loops = loops_by_face
        .values()
        .flat_map(|loops| loops.iter())
        .map(|loop_facts| (loop_facts.loop_handle, loop_facts.face))
        .collect::<BTreeMap<_, _>>();
    let mut adjacencies = Vec::new();
    for loops in loops_by_face.values() {
        for loop_facts in loops {
            for &half_edge in loop_facts.boundary_half_edges() {
                let Some(data) = brep.half_edges.get(half_edge) else {
                    continue;
                };
                let Some(twin) = brep.half_edges.get(data.twin) else {
                    continue;
                };
                let Some(adjacent_face) = seeded_loops.get(&twin.loop_).copied() else {
                    continue;
                };
                adjacencies.push(LoopAdjacencyFacts {
                    face: loop_facts.face,
                    loop_handle: loop_facts.loop_handle,
                    half_edge,
                    adjacent_face,
                    adjacent_loop: twin.loop_,
                    adjacent_half_edge: data.twin,
                });
            }
        }
    }
    adjacencies.sort_by_key(|adjacency| {
        (
            adjacency.face,
            adjacency.loop_handle,
            adjacency.half_edge,
            adjacency.adjacent_face,
            adjacency.adjacent_loop,
            adjacency.adjacent_half_edge,
        )
    });
    adjacencies
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
