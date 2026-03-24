//! Loop walking for face fragment extraction.
//!
//! Walks the intersection graph to extract closed loops of half-edges that
//! form face fragments. Each fragment consists of an outer loop (CCW) and
//! zero or more inner loops (CW holes).

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use std::collections::HashSet;

use popkern_core::geometry::Vec3;

use super::graph::{FaceRef, IntersectionGraph};

/// A face fragment — loops forming a subset of an original face.
#[derive(Debug)]
pub(crate) struct Fragment {
    /// The face this fragment belongs to.
    pub face: FaceRef,
    /// Half-edge indices forming the outer loop (CCW).
    pub outer_loop: Vec<usize>,
    /// Half-edge indices forming inner loops (CW holes).
    pub inner_loops: Vec<Vec<usize>>,
}

/// Walk all face fragments in the graph.
///
/// For each face in the graph's `face_edges`, walks closed loops of
/// half-edges using the radial ordering at each node. Returns one
/// [`Fragment`] per distinct outer loop found.
pub(crate) fn walk_all_fragments(graph: &IntersectionGraph) -> Vec<Fragment> {
    let mut fragments = Vec::new();

    for (face_ref, he_indices) in &graph.face_edges {
        let face_hes: HashSet<usize> = he_indices.iter().copied().collect();
        let mut visited: HashSet<usize> = HashSet::new();
        let mut loops: Vec<Vec<usize>> = Vec::new();

        for &start_he in he_indices {
            if visited.contains(&start_he) {
                continue;
            }

            let mut current_loop = Vec::new();
            let mut he = start_he;

            loop {
                if visited.contains(&he) {
                    break;
                }
                visited.insert(he);
                current_loop.push(he);

                // Find next half-edge around this face.
                match next_around_face(graph, he, face_ref, &face_hes) {
                    Some(next_he) => he = next_he,
                    None => break,
                }

                if he == start_he {
                    break;
                }
            }

            if current_loop.len() >= 2 {
                loops.push(current_loop);
            }
        }

        if loops.is_empty() {
            continue;
        }

        // Classify loops by signed area.
        let mut outer_idx = 0;
        let mut best_area = f64::NEG_INFINITY;

        for (i, lp) in loops.iter().enumerate() {
            let area = signed_area(graph, lp);
            // The loop with the largest absolute area is the outer loop.
            if area.abs() > best_area.abs() {
                best_area = area;
                outer_idx = i;
            }
        }

        let outer_loop = loops.remove(outer_idx);

        fragments.push(Fragment {
            face: *face_ref,
            outer_loop,
            inner_loops: loops,
        });
    }

    fragments
}

/// Find the next half-edge around a face after `he_idx`.
///
/// At the end node of `he_idx`, looks at the radial ordering of outgoing
/// edges and finds the next one that borders the same face.
fn next_around_face(
    graph: &IntersectionGraph,
    he_idx: usize,
    face: &FaceRef,
    face_hes: &HashSet<usize>,
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let outgoing = &graph.nodes[end_node].outgoing;

    // Find twin_idx's position in the radial ring.
    let twin_pos = outgoing.iter().position(|&idx| idx == twin_idx)?;

    // Walk the ring starting from after twin, looking for next edge bordering face.
    let n = outgoing.len();
    for offset in 1..n {
        let candidate = outgoing[(twin_pos + offset) % n];
        if face_hes.contains(&candidate) && graph.half_edges[candidate].face.as_ref() == Some(face)
        {
            return Some(candidate);
        }
    }

    None
}

/// Compute a signed area proxy for a loop using Newell's method.
///
/// Positive = CCW (outer), negative = CW (inner/hole).
fn signed_area(graph: &IntersectionGraph, loop_hes: &[usize]) -> f64 {
    // Sum cross products of consecutive vertex positions projected
    // onto an approximate normal.
    let positions: Vec<_> = loop_hes
        .iter()
        .map(|&he_idx| graph.nodes[graph.half_edges[he_idx].origin].position)
        .collect();

    if positions.len() < 3 {
        return 0.0;
    }

    // Newell's method for 3D polygon normal (proportional to area).
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;

    let n = positions.len();
    for i in 0..n {
        let cur = &positions[i];
        let next = &positions[(i + 1) % n];
        nx += (cur.y - next.y) * (cur.z + next.z);
        ny += (cur.z - next.z) * (cur.x + next.x);
        nz += (cur.x - next.x) * (cur.y + next.y);
    }

    let normal = Vec3::new(nx, ny, nz);
    // The length of the normal vector is proportional to twice the area.
    // The sign indicates orientation.
    normal.length() * normal.z.signum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph, NodeSource, Operand};
    use popkern_core::geometry::Point3;

    /// Build a triangle face in the graph with 3 nodes and 3 edge pairs.
    fn triangle_graph() -> (IntersectionGraph, FaceRef) {
        let mut graph = IntersectionGraph::new();

        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let n2 = graph.add_node(Point3::new(0.5, 1.0, 0.0), NodeSource::Intersection);

        let eh = popkern_core::topology::EdgeHandle::default();
        let face_ref = FaceRef {
            operand: Operand::A,
            face: popkern_core::topology::FaceHandle::default(),
        };

        let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
        let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
        let (fwd20, _) = graph.add_half_edge_pair(n2, n0, EdgeSource::OriginalA(eh), None);

        graph.set_half_edge_face(fwd01, face_ref);
        graph.set_half_edge_face(fwd12, face_ref);
        graph.set_half_edge_face(fwd20, face_ref);

        (graph, face_ref)
    }

    #[test]
    fn triangle_produces_one_outer_loop() {
        let (graph, _face_ref) = triangle_graph();

        let fragments = walk_all_fragments(&graph);

        assert_eq!(fragments.len(), 1, "Expected 1 fragment for triangle");
        assert_eq!(
            fragments[0].outer_loop.len(),
            3,
            "Expected 3 half-edges in outer loop"
        );
        assert!(fragments[0].inner_loops.is_empty());
    }

    #[test]
    fn square_face_produces_one_outer_loop() {
        let mut graph = IntersectionGraph::new();

        let n0 = graph.add_node(Point3::new(0.0, 0.0, 0.0), NodeSource::Intersection);
        let n1 = graph.add_node(Point3::new(1.0, 0.0, 0.0), NodeSource::Intersection);
        let n2 = graph.add_node(Point3::new(1.0, 1.0, 0.0), NodeSource::Intersection);
        let n3 = graph.add_node(Point3::new(0.0, 1.0, 0.0), NodeSource::Intersection);

        let eh = popkern_core::topology::EdgeHandle::default();
        let face_ref = FaceRef {
            operand: Operand::A,
            face: popkern_core::topology::FaceHandle::default(),
        };

        let (fwd01, _) = graph.add_half_edge_pair(n0, n1, EdgeSource::OriginalA(eh), None);
        let (fwd12, _) = graph.add_half_edge_pair(n1, n2, EdgeSource::OriginalA(eh), None);
        let (fwd23, _) = graph.add_half_edge_pair(n2, n3, EdgeSource::OriginalA(eh), None);
        let (fwd30, _) = graph.add_half_edge_pair(n3, n0, EdgeSource::OriginalA(eh), None);

        graph.set_half_edge_face(fwd01, face_ref);
        graph.set_half_edge_face(fwd12, face_ref);
        graph.set_half_edge_face(fwd23, face_ref);
        graph.set_half_edge_face(fwd30, face_ref);

        let fragments = walk_all_fragments(&graph);

        assert_eq!(fragments.len(), 1, "Expected 1 fragment for square");
        assert_eq!(
            fragments[0].outer_loop.len(),
            4,
            "Expected 4 half-edges in outer loop"
        );
        assert!(fragments[0].inner_loops.is_empty());
    }
}
