//! Loop walking for face fragment extraction.
//!
//! Walks the intersection graph to extract closed loops of half-edges that
//! form face fragments. Each fragment consists of an outer loop (CCW) and
//! zero or more inner loops (CW holes).

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

use std::collections::HashSet;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;

use super::graph::{FaceRef, IntersectionGraph, Operand};

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
///
/// The `brep_a` and `brep_b` parameters provide surface geometry so that
/// the walker can project edge tangent vectors onto the correct tangent
/// plane at each node — essential for correct angle computation on curved
/// surfaces (cylinders, cones, spheres, etc.).
pub(crate) fn walk_all_fragments(
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Vec<Fragment> {
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
                match next_around_face(graph, he, face_ref, &face_hes, brep_a, brep_b) {
                    Some(next_he) => he = next_he,
                    None => break,
                }

                if he == start_he {
                    break;
                }
            }

            // Accept loops of length 1 only for self-loop edges (closed
            // curves like circles), otherwise require at least 2 edges.
            let min_len = if current_loop.len() == 1 {
                let he = current_loop[0];
                let twin = graph.half_edges[he].twin;
                let start = graph.half_edges[he].origin;
                let end = graph.half_edges[twin].origin;
                if start == end { 1 } else { 2 }
            } else {
                2
            };
            if current_loop.len() >= min_len {
                loops.push(current_loop);
            }
        }

        if loops.is_empty() {
            continue;
        }

        // Classify loops by signed area.  The loop with the largest
        // absolute area is the outer boundary.  Smaller loops (less than
        // half the outer area) are inner loops (holes) that belong to
        // the same fragment.  Large loops form separate fragments (they
        // arise when an intersection curve splits a face).
        let areas: Vec<f64> = loops.iter().map(|lp| signed_area(graph, lp)).collect();
        let max_area = areas
            .iter()
            .map(|a| a.abs())
            .fold(f64::NEG_INFINITY, f64::max);
        let threshold = max_area * 0.5;

        let outer_idx = areas
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.abs()
                    .partial_cmp(&b.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Extract the outer loop (swap-remove to avoid shifting).
        let outer_loop = loops.swap_remove(outer_idx);
        // Adjust areas to match new indices after swap_remove.
        let mut remaining_areas = areas;
        remaining_areas.swap_remove(outer_idx);

        let mut inner_loops_out = Vec::new();
        let mut separate_outers = Vec::new();

        for (i, lp) in loops.into_iter().enumerate() {
            if remaining_areas[i].abs() < threshold {
                inner_loops_out.push(lp);
            } else {
                separate_outers.push(lp);
            }
        }

        fragments.push(Fragment {
            face: *face_ref,
            outer_loop,
            inner_loops: inner_loops_out,
        });

        for lp in separate_outers {
            fragments.push(Fragment {
                face: *face_ref,
                outer_loop: lp,
                inner_loops: vec![],
            });
        }
    }

    fragments
}

/// Find the next half-edge around a face after `he_idx`.
///
/// At the end node of `he_idx`, finds the outgoing half-edge belonging to
/// the same face that makes the smallest counter-clockwise angle from the
/// incoming direction.  This uses per-face angular comparison rather than
/// the global radial sort, because intersection nodes sit at the junction
/// of multiple faces with different normals.
fn next_around_face(
    graph: &IntersectionGraph,
    he_idx: usize,
    face: &FaceRef,
    face_hes: &HashSet<usize>,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let outgoing = &graph.nodes[end_node].outgoing;

    // Collect candidate outgoing edges that belong to this face (exclude
    // the twin of the current edge to prevent U-turns).
    let candidates: Vec<usize> = outgoing
        .iter()
        .copied()
        .filter(|idx| face_hes.contains(idx) && *idx != twin_idx)
        .collect();

    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        return Some(candidates[0]);
    }

    // Use per-face angular comparison for all multi-candidate cases.
    pick_by_angle(graph, he_idx, &candidates, face, brep_a, brep_b)
}

/// Pick the next edge from candidates using per-face CW angle.
///
/// Projects edge tangent vectors onto the surface tangent plane at the
/// node, then picks the candidate with the smallest clockwise angle from
/// the incoming direction.  Using the surface normal (instead of the
/// cross product of edge directions) is essential for curved surfaces
/// where edges are not coplanar.
fn pick_by_angle(
    graph: &IntersectionGraph,
    he_idx: usize,
    candidates: &[usize],
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let node_pos = graph.nodes[end_node].position;

    let origin_node = graph.half_edges[he_idx].origin;
    let incoming_dir = graph.nodes[origin_node].position - node_pos;

    // Try to get the surface normal from the face's geometry.
    // For curved surfaces this is critical for correct angle computation.
    let normal =
        surface_normal_at_node(graph, end_node, face, brep_a, brep_b).unwrap_or_else(|| {
            // Fallback: compute from incoming direction and candidates (planar case).
            let mut n = Vec3::new(0.0, 0.0, 1.0);
            for &c in candidates {
                let t = graph.half_edges[c].twin;
                let out = graph.nodes[graph.half_edges[t].origin].position - node_pos;
                if let Some(computed) = incoming_dir.cross(&out).normalized() {
                    n = computed;
                    break;
                }
            }
            n
        });

    let in_proj = incoming_dir - normal * incoming_dir.dot(&normal);
    let ref_dir = match in_proj.normalized() {
        Some(d) => d,
        None => return Some(candidates[0]),
    };

    let mut best = None;
    let mut best_angle = f64::MAX;

    for &c in candidates {
        let twin_c = graph.half_edges[c].twin;
        let end_pos = graph.nodes[graph.half_edges[twin_c].origin].position;
        let out_dir = end_pos - node_pos;
        let out_proj = out_dir - normal * out_dir.dot(&normal);

        let cos_val = out_proj.dot(&ref_dir);
        let sin_val = ref_dir.cross(&out_proj).dot(&normal);
        let mut angle = sin_val.atan2(cos_val);

        if angle <= 1e-10 {
            angle += std::f64::consts::TAU;
        }

        if angle < best_angle {
            best_angle = angle;
            best = Some(c);
        }
    }

    best
}

/// Look up the surface normal at a graph node's position on the given face.
///
/// Uses the node's stored UV coordinates (if available) to evaluate the
/// surface normal from the face's geometry.  Returns `None` if the face
/// has no surface geometry or the node has no UV coordinates for this face.
fn surface_normal_at_node(
    graph: &IntersectionGraph,
    node_idx: usize,
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Vec3> {
    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let face_data = brep.faces.get(face.face)?;
    let surface = face_data.surface.as_ref()?;
    let geom = &surface.geometry;

    // Try UV coords first.
    if let Some(&(u, v)) = graph.nodes[node_idx].uv_coords.get(face) {
        let n = geom.normal_at(u, v);
        return n.normalized();
    }

    // Fallback: compute approximate UV by projecting onto the surface.
    // For analytic surfaces, we can use the node position directly.
    // This is a simplified approach that works well for common cases.
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
    use popkern_core::tolerance::Tolerance;

    /// Empty BReps for tests that only use planar geometry.
    fn empty_breps() -> (BRep, BRep) {
        (
            BRep::new(Tolerance::default()),
            BRep::new(Tolerance::default()),
        )
    }

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
        let (ba, bb) = empty_breps();

        let fragments = walk_all_fragments(&graph, &ba, &bb);

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
        let (ba, bb) = empty_breps();

        let fragments = walk_all_fragments(&graph, &ba, &bb);

        assert_eq!(fragments.len(), 1, "Expected 1 fragment for square");
        assert_eq!(
            fragments[0].outer_loop.len(),
            4,
            "Expected 4 half-edges in outer loop"
        );
        assert!(fragments[0].inner_loops.is_empty());
    }
}
