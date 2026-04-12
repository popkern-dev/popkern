//! Loop walking for face fragment extraction.
//!
//! Walks the intersection graph to extract closed loops of half-edges that
//! form face fragments. Each fragment consists of an outer loop (CCW) and
//! zero or more inner loops (CW holes).

// Used by later pipeline stages; not yet called from the main merge path.
#![allow(dead_code)]

mod walk_geometry;

use std::collections::HashSet;

use popkern_core::topology::BRep;

use super::graph::{FaceRef, IntersectionGraph};
use walk_geometry::{pick_by_angle, signed_area};

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

        // Sort starting edges: seam-like first, then other non-shared, then shared.
        // Seam-like edges (where the twin is also on this face) form the outer
        // boundary "spine" that connects self-loop circles. Starting from them
        // ensures the outer boundary walk includes expanded cap circles rather
        // than consuming them as standalone loops.
        let mut sorted_hes: Vec<usize> = he_indices.clone();
        sorted_hes.sort_by_key(|&idx| {
            let is_shared = matches!(
                graph.half_edges[idx].source,
                super::graph::EdgeSource::Shared { .. }
            );
            // Primary key: seam-like (0) → non-shared (1) → shared (2).
            // Secondary key: half-edge index (lower first) for deterministic
            // fragment ordering — ensures the walk starts from the same edge
            // regardless of face_edges insertion order.
            let priority = if is_shared {
                2u8
            } else {
                let twin = graph.half_edges[idx].twin;
                let twin_on_same_face = graph.half_edges[twin].face.as_ref() == Some(face_ref);
                if twin_on_same_face { 0 } else { 1 }
            };
            (priority, idx)
        });

        for &start_he in &sorted_hes {
            if visited.contains(&start_he) {
                continue;
            }
            // Skip cross-registered self-loop half-edges as starting
            // points. They appear in face_edges due to cross-registration
            // but should not form standalone loops on this face.
            let twin_he = graph.half_edges[start_he].twin;
            let is_self_loop =
                graph.half_edges[start_he].origin == graph.half_edges[twin_he].origin;
            if is_self_loop && graph.half_edges[start_he].face.as_ref() != Some(face_ref) {
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
                match next_around_face(
                    graph,
                    he,
                    face_ref,
                    &face_hes,
                    brep_a,
                    brep_b,
                    Some(start_he),
                    &visited,
                    &current_loop,
                ) {
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

        // Check if any loop contains self-loop edges (closed curves).
        // When self-loops are present, the Newell area for the boundary
        // loop can be zero (all vertices on the seam line). In this case,
        // use edge count as the primary metric: the longest loop is the
        // outer boundary.
        let has_self_loop_edges = loops.iter().any(|lp| {
            lp.iter().any(|&he_idx| {
                let twin = graph.half_edges[he_idx].twin;
                graph.half_edges[he_idx].origin == graph.half_edges[twin].origin
            })
        });

        let outer_idx = if has_self_loop_edges {
            loops
                .iter()
                .enumerate()
                .max_by_key(|(_, lp)| lp.len())
                .map(|(i, _)| i)
                .unwrap_or(0)
        } else {
            areas
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| {
                    a.abs()
                        .partial_cmp(&b.abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0)
        };

        // Extract the outer loop (swap-remove to avoid shifting).
        let outer_loop = loops.swap_remove(outer_idx);
        // Adjust areas to match new indices after swap_remove.
        let mut remaining_areas = areas;
        remaining_areas.swap_remove(outer_idx);

        let mut inner_loops_out: Vec<Vec<usize>> = Vec::new();
        let mut separate_outers = Vec::new();

        // Collect the half-edge indices in the outer loop for
        // twin-overlap detection.
        let outer_he_set: HashSet<usize> = outer_loop.iter().copied().collect();

        for (i, lp) in loops.into_iter().enumerate() {
            // A loop whose half-edges share an underlying graph edge
            // with the outer loop is a face-split boundary, not a hole.
            // Detect this by checking if any half-edge in the loop has
            // its twin in the outer loop.
            let shares_edge_with_outer = lp.iter().any(|&he_idx| {
                let twin = graph.half_edges[he_idx].twin;
                outer_he_set.contains(&twin)
            });
            if shares_edge_with_outer || remaining_areas[i].abs() >= threshold {
                separate_outers.push(lp);
            } else {
                // Skip duplicate inner loops: if this loop's edges are all
                // twins of an already-collected inner loop, it's a duplicate
                // (both directions of the same hole boundary were walked).
                let is_twin_duplicate = inner_loops_out.iter().any(|existing| {
                    existing.len() == lp.len()
                        && lp.iter().all(|&he_idx| {
                            let twin = graph.half_edges[he_idx].twin;
                            existing.contains(&twin)
                        })
                });
                if !is_twin_duplicate {
                    inner_loops_out.push(lp);
                }
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
///
/// **Cut behavior:** when `he_idx` is a closed self-loop on a `Shared`
/// edge — i.e. an intersection circle that the boolean inserted as a
/// cut through the face — and the previous edge in `current_loop`
/// has its twin among the candidates, the walker bounces back along
/// that twin. This treats the intersection circle as a topological
/// cut: walk in via some boundary edge, traverse the circle, walk
/// back out the same way. Without this rule, the walker keeps going
/// past the cut and merges multiple regions into a single mega-loop
/// (the cube-minus-cylinder Euler-violation bug).
#[allow(clippy::too_many_arguments)]
fn next_around_face(
    graph: &IntersectionGraph,
    he_idx: usize,
    face: &FaceRef,
    face_hes: &HashSet<usize>,
    brep_a: &BRep,
    brep_b: &BRep,
    start_he: Option<usize>,
    visited: &HashSet<usize>,
    current_loop: &[usize],
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let outgoing = &graph.nodes[end_node].outgoing;

    // For self-loop edges (origin == twin.origin), exclude both the current
    // edge AND its twin to avoid re-traversing the same self-loop or U-turning.
    // For normal edges, only exclude the twin (standard U-turn prevention).
    let is_self_loop = graph.half_edges[he_idx].origin == end_node;

    let candidates: Vec<usize> = outgoing
        .iter()
        .copied()
        .filter(|&idx| {
            if !face_hes.contains(&idx) {
                return false;
            }
            // Allow the start_he (loop closure) even if it would normally
            // be excluded by U-turn prevention or visited check.
            if start_he == Some(idx) {
                return true;
            }
            // Exclude already-visited half-edges to prevent the walk
            // from crossing into edges consumed by a previous iteration.
            if visited.contains(&idx) {
                return false;
            }
            // Exclude the twin (U-turn prevention)
            if idx == twin_idx {
                return false;
            }
            // If we're on a self-loop, also exclude ourselves
            if is_self_loop && idx == he_idx {
                return false;
            }
            // When walking a NON-self-loop path, exclude self-loop
            // candidates from other faces. Self-loops whose face matches
            // the current face are boundary edges (e.g., intersection
            // circles on the barrel face); cross-registered self-loops
            // from other faces should still be excluded.
            if !is_self_loop {
                let twin_c = graph.half_edges[idx].twin;
                let is_candidate_self_loop =
                    graph.half_edges[idx].origin == graph.half_edges[twin_c].origin;
                if is_candidate_self_loop && graph.half_edges[idx].face.as_ref() != Some(face) {
                    return false;
                }
            }
            true
        })
        .collect();

    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        return Some(candidates[0]);
    }

    // §7.4 cut behavior: when the previous edge is a closed self-loop
    // on a Shared edge (an intersection circle that the boolean
    // inserted as a cut through this face), bounce back along the
    // twin of the edge that brought us TO this circle. The
    // intersection circle acts as a topological barrier — without
    // this rule, the walker continues past the cut and merges
    // multiple disjoint regions into a single mega-loop.
    //
    // The detection: `he_idx` is a closed self-loop (origin == end),
    // its source is `Shared`, and the SECOND-most-recent edge in
    // `current_loop` (the one before he_idx) has its twin in the
    // candidate set.
    if is_self_loop_shared(graph, he_idx)
        && current_loop.len() >= 2
        && let Some(prev_prev) = current_loop.get(current_loop.len().wrapping_sub(2))
    {
        let prev_prev_twin = graph.half_edges[*prev_prev].twin;
        if candidates.contains(&prev_prev_twin) {
            return Some(prev_prev_twin);
        }
    }

    // Prefer loop closure when start_he is a candidate AND all other
    // candidates are already visited. This avoids closing the loop
    // prematurely (skipping unvisited boundary edges like cap circles)
    // while still preventing extension into cross-registered duplicates.
    if let Some(sh) = start_he
        && candidates.contains(&sh)
    {
        // Only defer closure if there are unvisited candidates whose face
        // matches (actual boundary edges, not cross-registered duplicates).
        // Only defer closure for unvisited NON-SHARED face-matching
        // candidates (actual boundary edges like cap circles). Shared
        // edges (even face-matching ones from cross-registration) should
        // not prevent loop closure.
        let unvisited_face_matching = candidates.iter().any(|&c| {
            c != sh
                && !visited.contains(&c)
                && graph.half_edges[c].face.as_ref() == Some(face)
                && !matches!(
                    graph.half_edges[c].source,
                    super::graph::EdgeSource::Shared { .. }
                )
        });
        if !unvisited_face_matching {
            return Some(sh);
        }
    }

    // At seam T-junctions on curved surfaces (where seam continuation
    // competes with a self-loop circle), prefer the unvisited self-loop.
    // Self-loops with matching face are actual boundary edges; the seam
    // is just the internal spine of a periodic surface.
    let self_loop_candidates: Vec<usize> = candidates
        .iter()
        .copied()
        .filter(|&c| {
            !visited.contains(&c) && {
                let twin_c = graph.half_edges[c].twin;
                graph.half_edges[c].origin == graph.half_edges[twin_c].origin
            }
        })
        .collect();
    if !self_loop_candidates.is_empty() && self_loop_candidates.len() < candidates.len() {
        if self_loop_candidates.len() == 1 {
            return Some(self_loop_candidates[0]);
        }
        return pick_by_angle(graph, he_idx, &self_loop_candidates, face, brep_a, brep_b);
    }

    // Use per-face angular comparison for all multi-candidate cases.
    pick_by_angle(graph, he_idx, &candidates, face, brep_a, brep_b)
}

/// True if `he_idx` is a closed self-loop half-edge whose source is
/// `Shared` (an intersection circle).
fn is_self_loop_shared(graph: &IntersectionGraph, he_idx: usize) -> bool {
    let he = &graph.half_edges[he_idx];
    let twin = &graph.half_edges[he.twin];
    he.origin == twin.origin && matches!(he.source, super::graph::EdgeSource::Shared { .. })
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
