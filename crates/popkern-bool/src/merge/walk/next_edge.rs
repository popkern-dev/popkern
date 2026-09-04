//! Per-node next-edge selection for the 3D-Newell fragment walk.
//!
//! Split from `walk/mod.rs` under the structural line ratchet.

use std::collections::HashSet;

use popkern_core::topology::BRep;

use super::boundary_predicates::{self, is_self_loop_shared, node_on_face_boundary};
use super::walk_geometry::{WalkAxis, pick_by_angle, pick_by_angle_ext};
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph};

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
pub(super) fn next_around_face(
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
                // ZERO-EXTENT self-loops are point-contact records (a
                // tangency between two faces). They bound nothing, and
                // entering one hijacks the host loop: the closure
                // deferral held a ceiling loop open at its closure
                // node, the walk swallowed the tangency plus the twin
                // of the flank it had just traversed, and the spike's
                // flipped claim collided with the neighbouring flank
                // fragment's (measured on the uniform-width flush
                // arm). Their seed walks become one-half-edge
                // zero-extent fragments, which the degenerate cull
                // already removes.
                if is_candidate_self_loop {
                    use popkern_core::geometry::curve::CurveEval as _;
                    let extent = graph.half_edges[idx]
                        .curve
                        .as_ref()
                        .or_else(|| graph.half_edges[twin_c].curve.as_ref())
                        .map(|c| {
                            let origin = graph.nodes[graph.half_edges[idx].origin].position;
                            (1..8)
                                .map(|i| {
                                    let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                                    c.geometry.point_at(t).distance_to(&origin)
                                })
                                .fold(0.0, f64::max)
                        })
                        .unwrap_or(0.0);
                    if extent < 1e-9 {
                        return false;
                    }
                }
                // A self-loop candidate belongs to this walk when
                // EITHER half is assigned to this face: intersection
                // circles are cross-registered exactly so both of a
                // face's fragments (each side of the cut) can use one
                // direction each. Excluding by primary face alone
                // starved the second fragment of its circle boundary.
                let this_face = graph.half_edges[idx].face.as_ref() == Some(face)
                    || graph.half_edges[twin_c].face.as_ref() == Some(face);
                if is_candidate_self_loop && !this_face {
                    return false;
                }
                // INTERIOR-GRAZE FILTER: a candidate that dead-ends
                // (its far node has no onward half-edge on this face)
                // comes in two shapes that need opposite handling.
                // A HALF-SPLIT trim's far endpoint lies ON one of the
                // face's boundary edges — the boundary just was not
                // split there in this phase. Entering it is required:
                // the walk records the cut and downstream completes it
                // (measured: a cube bisected by a plane sheet stopped
                // splitting when these were filtered). An INTERIOR
                // graze record's far endpoint dead-ends strictly
                // inside the face — an embedded body's face tracing a
                // segment whose both sides bound the SAME kept
                // material (a flush gusset top's trace on the arm
                // wall). No result edge should exist there and no rule
                // may enter it: the junction rule already declined
                // these spurs, but the angular fallback picked them
                // anyway, weaving a 10-he unclosable outer that minted
                // a once-used edge.
                if !is_candidate_self_loop {
                    let far = graph.half_edges[twin_c].origin;
                    let has_onward = graph.nodes[far]
                        .outgoing
                        .iter()
                        .any(|&o| o != twin_c && o != idx && face_hes.contains(&o));
                    if !has_onward && !node_on_face_boundary(graph, far, face_hes, idx, twin_c) {
                        return false;
                    }
                }
            }
            true
        })
        .collect();

    if popkern_core::trace::trace_enabled("fragment-walk") {
        let p = graph.nodes[end_node].position;
        eprintln!(
            "    [next] at ({:.3},{:.3},{:.3}) from he {he_idx}: candidates {candidates:?}",
            p.x, p.y, p.z
        );
    }
    if candidates.is_empty() {
        // APEX TURNAROUND (PK-343) — see `apex_turnaround_candidate`.
        return boundary_predicates::apex_turnaround_candidate(
            graph, end_node, he_idx, twin_idx, face_hes, visited,
        );
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
                && !matches!(graph.half_edges[c].source, EdgeSource::Shared { .. })
        });
        // An unvisited SELF-LOOP intersection curve at this node is a
        // topological barrier that bounds this fragment — closing the
        // loop early would leave the fragment without its cut
        // boundary (the wall of a boss/plate union lost its seam
        // circle and the circle became a stray one-edge fragment).
        let unvisited_cut_here = candidates.iter().any(|&c| {
            let twin_c = graph.half_edges[c].twin;
            c != sh
                && !visited.contains(&c)
                && graph.half_edges[c].origin == graph.half_edges[twin_c].origin
                && (graph.half_edges[c].face.as_ref() == Some(face)
                    || graph.half_edges[twin_c].face.as_ref() == Some(face))
        });
        if !unvisited_face_matching && !unvisited_cut_here {
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
        // Choosing between the TWO DIRECTIONS of one circle is a
        // CLAIM, not a turn: either direction closes the same region
        // (assembly's flip fixes winding), but the claim must
        // anti-align with what the neighbouring face's walker claims
        // for the twin. The chart walker and the assembly's reversal
        // machinery are tuned around the historical outward-axis
        // claim, so this pick keeps that hand (measured both ways:
        // the material-left axis claimed the boss-wall rim direction
        // the plate chart had taken — duplicate minted edge; a
        // face-registration rule inverted the torus-equator claims
        // and starved the annulus cap of both rims). Transversal
        // junction TURNS keep the material-left axis — there the
        // candidates bound different regions and the hand decides
        // which region the loop closes.
        return pick_by_angle_ext(
            graph,
            he_idx,
            &self_loop_candidates,
            face,
            brep_a,
            brep_b,
            WalkAxis::OutwardClaim,
        );
    }

    // TOPOLOGICAL JUNCTION RULE (convention-free): a transversal
    // boundary-to-boundary trim separates two fragments, so every
    // loop through one of its junction nodes must ALTERNATE — a walk
    // arriving along the BOUNDARY turns onto the trim (it bounds the
    // fragment on this side), and a walk arriving along the TRIM
    // continues onto the boundary. Continuing straight past a
    // crossing trim encloses both fragments in one loop (measured: a
    // rest-plane trim on a reflex block/plate junction wall was left
    // unused, both wall regions merged into one 10-half-edge
    // fragment carrying removed-region geometry, and the trim's own
    // walk dead-ended). This is exact whenever the junction has
    // exactly one non-self-loop Shared candidate against non-shared
    // boundary continuations — no surface normal or angle involved.
    // Branch points (several trims) and circle junctions keep the
    // angular pick below.
    let incoming_shared =
        matches!(graph.half_edges[he_idx].source, EdgeSource::Shared { .. }) && !is_self_loop;
    if !incoming_shared {
        let shared_straight: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|&c| {
                let twin_c = graph.half_edges[c].twin;
                let self_loop_c = graph.half_edges[c].origin == graph.half_edges[twin_c].origin;
                !self_loop_c && matches!(graph.half_edges[c].source, EdgeSource::Shared { .. })
            })
            .collect();
        if shared_straight.len() == 1 && shared_straight.len() < candidates.len() {
            // SPUR GUARD: a closed loop can never traverse a trim
            // whose far node has no onward half-edge on this face —
            // entering it forces a forbidden twin u-turn. A grazing
            // record (a flush gusset's top-edge line lying ON the
            // wall face) dead-ends mid-face; the rule walked the
            // outer boundary straight into it and the walk stuck.
            // Only take the trim when it connects onward.
            let c = shared_straight[0];
            let twin_c = graph.half_edges[c].twin;
            let far = graph.half_edges[twin_c].origin;
            let has_onward = graph.nodes[far]
                .outgoing
                .iter()
                .any(|&o| o != twin_c && o != c && face_hes.contains(&o));
            if has_onward {
                return Some(c);
            }
        }
    }

    // Use per-face angular comparison for all multi-candidate cases.
    pick_by_angle(graph, he_idx, &candidates, face, brep_a, brep_b)
}
