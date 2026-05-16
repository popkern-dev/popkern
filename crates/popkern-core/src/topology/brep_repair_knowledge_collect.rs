//! Shared read-only collectors for BRep repair knowledge.

use std::collections::{BTreeMap, BTreeSet};

use crate::topology::{
    BRep, BrepEndpointLoopMembership, BrepEndpointMembershipFacts, BrepLoopAdjacencyFacts,
    BrepLoopFacts, BrepLoopKind, FaceHandle, LoopHandle, VertexHandle,
};

/// Collect loop facts for every loop owned by a face.
pub fn collect_brep_face_loops(brep: &BRep, face_handle: FaceHandle) -> Vec<BrepLoopFacts> {
    let Some(face) = brep.faces.get(face_handle) else {
        return Vec::new();
    };
    let mut loop_handles = Vec::with_capacity(1 + face.inner_loops.len());
    loop_handles.push((face.outer_loop, BrepLoopKind::Outer));
    loop_handles.extend(
        face.inner_loops
            .iter()
            .copied()
            .map(|loop_handle| (loop_handle, BrepLoopKind::Inner)),
    );
    loop_handles
        .into_iter()
        .map(|(loop_handle, kind)| collect_brep_loop_facts(brep, face_handle, loop_handle, kind))
        .collect()
}

/// Collect facts for one loop walk.
pub fn collect_brep_loop_facts(
    brep: &BRep,
    face: FaceHandle,
    loop_handle: LoopHandle,
    kind: BrepLoopKind,
) -> BrepLoopFacts {
    let mut boundary_vertices = Vec::new();
    let mut boundary_half_edges = Vec::new();
    let mut seen_half_edges = BTreeSet::new();
    let mut repeated_half_edge = false;
    let mut repeated_vertex = false;
    let Some(loop_data) = brep.loops.get(loop_handle) else {
        return BrepLoopFacts {
            face,
            loop_handle,
            kind: BrepLoopKind::Unknown,
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
        if boundary_vertices.contains(&half_edge.origin) {
            repeated_vertex = true;
        }
        boundary_half_edges.push(cursor);
        boundary_vertices.push(half_edge.origin);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
    BrepLoopFacts {
        face,
        loop_handle,
        kind,
        boundary_half_edges,
        boundary_vertices,
        repeated_vertex,
        repeated_half_edge,
    }
}

/// Collect loop adjacencies whose twin loops are also present in the scoped map.
pub fn collect_brep_loop_adjacencies(
    brep: &BRep,
    loops_by_face: &BTreeMap<FaceHandle, Vec<BrepLoopFacts>>,
) -> Vec<BrepLoopAdjacencyFacts> {
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
                adjacencies.push(BrepLoopAdjacencyFacts {
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

/// Collect faces adjacent to the supplied scoped loop facts.
pub fn collect_brep_adjacent_faces(
    brep: &BRep,
    loops_by_face: &BTreeMap<FaceHandle, Vec<BrepLoopFacts>>,
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

/// Collect endpoint loop memberships for a supplied endpoint vertex set.
pub fn collect_brep_endpoint_memberships(
    loops_by_face: &BTreeMap<FaceHandle, Vec<BrepLoopFacts>>,
    endpoint_vertices: impl IntoIterator<Item = VertexHandle>,
) -> Vec<BrepEndpointMembershipFacts> {
    endpoint_vertices
        .into_iter()
        .map(|vertex| BrepEndpointMembershipFacts {
            vertex,
            containing_loops: loops_by_face
                .values()
                .flat_map(|loops| loops.iter())
                .filter_map(|loop_facts| {
                    let position = loop_facts
                        .boundary_vertices
                        .iter()
                        .position(|loop_vertex| *loop_vertex == vertex)?;
                    Some(BrepEndpointLoopMembership {
                        face: loop_facts.face,
                        loop_handle: loop_facts.loop_handle,
                        half_edge: loop_facts.boundary_half_edges[position],
                        kind: loop_facts.kind,
                    })
                })
                .collect(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;

    #[test]
    fn unit_cube_face_loop_collector_returns_outer_loop() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has face");

        let loops = collect_brep_face_loops(&brep, face);

        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].face, face);
        assert_eq!(loops[0].kind, BrepLoopKind::Outer);
        assert!(!loops[0].boundary_half_edges.is_empty());
    }
}
