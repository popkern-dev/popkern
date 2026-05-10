//! Collector helpers for repeated-apex topology knowledge.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::repeated_apex_topology_knowledge::{LoopAdjacencyFacts, LoopFacts, LoopKind};

pub(super) fn collect_face_loops(
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
            super::repeated_apex_topology_knowledge::collect_loop_facts(
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

pub(super) fn collect_loop_adjacencies(
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

pub(super) fn collect_adjacent_faces(
    brep: &BRep,
    loops_by_face: &BTreeMap<FaceHandle, Vec<LoopFacts>>,
) -> Vec<FaceHandle> {
    let mut adjacent_faces = BTreeSet::new();
    for loops in loops_by_face.values() {
        for loop_facts in loops {
            for &half_edge in loop_facts.boundary_half_edges() {
                let Some(data) = brep.half_edges.get(half_edge) else {
                    continue;
                };
                let Some(twin) = brep.half_edges.get(data.twin) else {
                    continue;
                };
                let Some(twin_loop) = brep.loops.get(twin.loop_) else {
                    continue;
                };
                adjacent_faces.insert(twin_loop.face);
            }
        }
    }
    adjacent_faces.into_iter().collect()
}
