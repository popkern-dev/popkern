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
    walk_all_fragments_ext(graph, brep_a, brep_b, false)
}

/// `coincident_interface`: the boolean carries coplanar face
/// overlaps — enables the disjoint-loop election override, whose
/// sampled-polygon containment test is currently only trusted on
/// that configuration (it mis-voted a through-slot ring on a long
/// tapered face — see the flush-gusset pin).
pub(crate) fn walk_all_fragments_ext(
    graph: &IntersectionGraph,
    brep_a: &BRep,
    brep_b: &BRep,
    coincident_interface: bool,
) -> Vec<Fragment> {
    let mut fragments = Vec::new();

    for (face_ref, he_indices) in &graph.face_edges {
        if std::env::var_os("PK_WALK_TRACE").is_some() {
            let n_int = he_indices
                .iter()
                .filter(|&&h| {
                    matches!(
                        graph.half_edges[h].source,
                        crate::merge::graph::EdgeSource::Shared { .. }
                    )
                })
                .count();
            eprintln!(
                "  [walk-pre] {face_ref:?}: {} hes ({n_int} intersection)",
                he_indices.len()
            );
            if std::env::var_os("PK_WALK_TRACE").is_some() {
                for &h in he_indices.iter() {
                    let he = &graph.half_edges[h];
                    let o = he.origin;
                    let e = graph.half_edges[he.twin].origin;
                    let po = graph.nodes[o].position;
                    let pe = graph.nodes[e].position;
                    eprintln!(
                        "      [he] {h} n{o}({:.9},{:.9},{:.9}) -> n{e}({:.9},{:.9},{:.9}) twin {} src {:?} curve {}",
                        po.x,
                        po.y,
                        po.z,
                        pe.x,
                        pe.y,
                        pe.z,
                        he.twin,
                        he.source,
                        he.curve
                            .as_ref()
                            .or(graph.half_edges[he.twin].curve.as_ref())
                            .map(|c| format!("[{:.2},{:.2}]", c.t_start, c.t_end))
                            .unwrap_or_else(|| "none".into())
                    );
                }
            }
        }
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
                    Some(next_he) => {
                        if std::env::var_os("PK_WALK_TRACE").is_some() {
                            eprintln!("      [step] {he} -> {next_he}");
                        }
                        he = next_he;
                    }
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

        // Outer = largest |area|. signed_area samples curve geometry
        // for closed/self-loop edges, so a 1-half-edge closed rim
        // (cylinder cap) reports its true disk area. The historical
        // fallback — "longest loop wins when self-loops are present"
        // — chose a 4-edge interior pocket loop over the 1-edge rim
        // and every pocket cut into a closed-rim cap lost its hole
        // (measured: cap kept with inner 0, pocket walls orphaned,
        // 307 open mesh edges). Length only breaks area ties now.
        let outer_idx = {
            let mut best = 0usize;
            for i in 1..loops.len() {
                let ai = areas[i].abs();
                let ab = areas[best].abs();
                if ai > ab * (1.0 + 1e-9)
                    || ((ai - ab).abs() <= ab * 1e-9 && loops[i].len() > loops[best].len())
                {
                    best = i;
                }
            }
            best
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

        let outer_newell = walk_geometry::loop_newell(graph, &outer_loop);
        // The winding override below applies to WORKPIECE (operand A)
        // faces only. On the TOOL side a large interior ring's both
        // directions serve as candidate caps for the result pieces (a
        // plane-split's cutting box face contributes the cap to each
        // half), and pre-assigning one as the tool-frame's hole
        // starved a result piece of its cap (measured on the
        // solid-by-plane splitter).
        let face_on_workpiece = matches!(face_ref.operand, super::graph::Operand::A);
        for (i, lp) in loops.into_iter().enumerate() {
            // A loop whose half-edges share an underlying graph edge
            // with the outer loop is a face-split boundary, not a hole.
            // Detect this by checking if any half-edge in the loop has
            // its twin in the outer loop.
            let shares_edge_with_outer = lp.iter().any(|&he_idx| {
                let twin = graph.half_edges[he_idx].twin;
                outer_he_set.contains(&twin)
            });
            // Winding RESCUE for large PAIRED holes: the 0.5·area
            // threshold alone called a 76 mm seat ring on a small
            // block face a "separate fragment" purely because the
            // hole was big relative to the face. When BOTH directions
            // of an interior ring were walked, the one winding
            // OPPOSITE the outer (Newell dot < 0) is the hole and its
            // twin is the interior piece — regardless of size. An
            // UNPAIRED big loop keeps the size rule (a sheet-split's
            // interior piece walks once and must stay a piece), and
            // promoting winding to the primary test for small loops
            // was measured to break six tuned suites.
            // Winding override for LARGE holes: the 0.5·area
            // threshold alone called a 76 mm seat ring on a small
            // block face a "separate fragment" purely because the
            // hole was big relative to the face. A big loop winding
            // clearly OPPOSITE the outer (Newell dot < 0) is its
            // hole. Small loops keep the size rule — demoting small
            // same-winding loops was measured to break six suites.
            let separate_by_shape = if remaining_areas[i].abs() >= threshold {
                if face_on_workpiece {
                    let winding_dot = walk_geometry::loop_newell(graph, &lp).dot(&outer_newell);
                    winding_dot >= -1e-9
                } else {
                    true
                }
            } else {
                // CONTAINMENT override: a hole must lie INSIDE the
                // outer loop's region. The size proxy alone called a
                // disjoint root-collar loop a "hole" of a swoop
                // fragment on one mirror side and not the other
                // (area 263 vs threshold 352 — pure size accident).
                // A loop whose nodes are strictly outside the outer
                // polygon can never be its hole; loops inside keep
                // the existing rules.
                coincident_interface
                    && walk_geometry::loop_outside_outer(graph, &lp, &outer_loop, &outer_newell)
            };
            if std::env::var_os("PK_WALK_TRACE").is_some() {
                eprintln!(
                    "    [elect] loop {} hes area {:.1} thr {:.1} shares={} sep={}",
                    lp.len(),
                    remaining_areas[i],
                    threshold,
                    shares_edge_with_outer,
                    separate_by_shape
                );
            }
            if shares_edge_with_outer || separate_by_shape {
                separate_outers.push(lp);
            } else {
                // Both directions of an interior closed ring reach here.
                // They are NOT duplicates: one (winding OPPOSITE the
                // outer loop) is the hole boundary of this fragment;
                // the other (winding WITH the outer) is the OUTER
                // boundary of the interior piece — a genuine separate
                // fragment. Discarding it as a "duplicate" deleted the
                // interior piece wholesale (measured: a rotated slab's
                // cutting face lost the quad inside the prism, so the
                // subtracted wedge had no rest plane — shell gap = the
                // full prism width).
                let twin_of = |target: &[usize]| {
                    inner_loops_out.iter().position(|existing| {
                        existing.len() == target.len()
                            && target.iter().all(|&he_idx| {
                                let twin = graph.half_edges[he_idx].twin;
                                existing.contains(&twin)
                            })
                    })
                };
                match twin_of(&lp) {
                    None => inner_loops_out.push(lp),
                    // The first-walked direction stays the fragment's
                    // hole boundary — identical to the long-standing
                    // behavior when the twin was discarded, so every
                    // configuration that worked before is untouched.
                    // Role-reassignment heuristics (area sign, face
                    // ownership) were each measured to break a
                    // configuration the other fixed.
                    Some(_) => separate_outers.push(lp),
                }
            }
        }

        if std::env::var_os("PK_LOOP_DUMP").is_some() && outer_loop.len() >= 20 {
            eprintln!("  [loopdump] {face_ref:?} outer {} hes:", outer_loop.len());
            for &h in &outer_loop {
                let n = graph.half_edges[h].origin;
                let p = graph.nodes[n].position;
                let shared = matches!(
                    graph.half_edges[h].source,
                    super::graph::EdgeSource::Shared { .. }
                );
                eprintln!(
                    "    he {h} ({:.1},{:.1},{:.1}) shared={shared}",
                    p.x, p.y, p.z
                );
            }
        }
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [walk] {face_ref:?}: outer {} hes, {} inner {:?}, {} separate",
                outer_loop.len(),
                inner_loops_out.len(),
                inner_loops_out.iter().map(|l| l.len()).collect::<Vec<_>>(),
                separate_outers.len()
            );
        }
        if std::env::var_os("PK_BOOL_TRACE").is_some() && !inner_loops_out.is_empty() {
            eprintln!(
                "  [walk] fragment for {face_ref:?}: outer len {} + {} inner loops",
                outer_loop.len(),
                inner_loops_out.len()
            );
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
            }
            true
        })
        .collect();

    if std::env::var_os("PK_WALK_TRACE").is_some() {
        let p = graph.nodes[end_node].position;
        eprintln!(
            "    [next] at ({:.3},{:.3},{:.3}) from he {he_idx}: candidates {candidates:?}",
            p.x, p.y, p.z
        );
    }
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
        return pick_by_angle(graph, he_idx, &self_loop_candidates, face, brep_a, brep_b);
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
    let incoming_shared = matches!(
        graph.half_edges[he_idx].source,
        super::graph::EdgeSource::Shared { .. }
    ) && !is_self_loop;
    if !incoming_shared {
        let shared_straight: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|&c| {
                let twin_c = graph.half_edges[c].twin;
                let self_loop_c = graph.half_edges[c].origin == graph.half_edges[twin_c].origin;
                !self_loop_c
                    && matches!(
                        graph.half_edges[c].source,
                        super::graph::EdgeSource::Shared { .. }
                    )
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

/// True if `he_idx` is a closed self-loop half-edge whose source is
/// `Shared` (an intersection circle).
fn is_self_loop_shared(graph: &IntersectionGraph, he_idx: usize) -> bool {
    let he = &graph.half_edges[he_idx];
    let twin = &graph.half_edges[he.twin];
    he.origin == twin.origin && matches!(he.source, super::graph::EdgeSource::Shared { .. })
}

#[cfg(test)]
mod tests;
