//! Boundary-proximity predicates for the fragment walk:
//! stored-boundary containment for graze nodes and shared self-loop
//! detection.
//!
//! Split from `walk/mod.rs` to honour the structural line ratchet.

use std::collections::HashSet;

use crate::merge::graph::IntersectionGraph;

/// True when the graph node lies ON one of the face's registered
/// half-edge spans (excluding the candidate being tested and its
/// twin). Used to tell a HALF-SPLIT trim (far endpoint on an unsplit
/// boundary edge — walkable, downstream completes the cut) from an
/// INTERIOR graze record (far endpoint strictly inside the face —
/// never walkable). Straight spans use exact point-segment distance;
/// curved spans are sampled. The epsilon is relative to the span
/// length: the discrimination this feeds is coarse (a half-split
/// endpoint sits exactly on its host edge; a graze endpoint sits
/// millimetres inside the face).
pub(super) fn node_on_face_boundary(
    graph: &IntersectionGraph,
    node: usize,
    face_hes: &HashSet<usize>,
    exclude_a: usize,
    exclude_b: usize,
) -> bool {
    let p = graph.nodes[node].position;
    for &h in face_hes {
        if h == exclude_a || h == exclude_b {
            continue;
        }
        let he = &graph.half_edges[h];
        let a = graph.nodes[he.origin].position;
        let b = graph.nodes[graph.half_edges[he.twin].origin].position;
        let curved = he.curve.as_ref().map(|c| {
            !matches!(
                c.geometry,
                popkern_core::geometry::curve::CurveGeometry::Line { .. }
            )
        });
        if curved == Some(true)
            && let Some(c) = he.curve.as_ref()
        {
            use popkern_core::geometry::curve::CurveEval;
            let span = a.distance_to(&b).max(1e-12);
            for i in 0..=16 {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 16.0;
                if c.geometry.point_at(t).distance_to(&p) < 1e-6 * (1.0 + span) {
                    return true;
                }
            }
        } else {
            let ab = b - a;
            let l2 = ab.dot(&ab);
            if l2 < 1e-24 {
                continue;
            }
            let t = ((p - a).dot(&ab) / l2).clamp(0.0, 1.0);
            let foot = a + ab * t;
            if foot.distance_to(&p) < 1e-6 * (1.0 + l2.sqrt()) {
                return true;
            }
        }
    }
    false
}

/// True if `he_idx` is a closed self-loop half-edge whose source is
/// `Shared` (an intersection circle).
pub(super) fn is_self_loop_shared(graph: &IntersectionGraph, he_idx: usize) -> bool {
    let he = &graph.half_edges[he_idx];
    let twin = &graph.half_edges[he.twin];
    he.origin == twin.origin && matches!(he.source, crate::merge::graph::EdgeSource::Shared { .. })
}

/// APEX TURNAROUND (PK-343): at an ORIGINAL vertex whose only
/// registered departure is the twin — a cone apex, where the seam's
/// two halves are twins inside the SAME loop — the legal continuation
/// IS the twin: the loop doubles back by construction. U-turn
/// prevention exists to stop premature backtracking where
/// alternatives exist; refusing the turn here tore the loop
/// (measured: cone-cube union dropped the seam-to-apex piece, loop
/// gap 0.5099, served silently until the PK-342 closure gate).
/// Half-split trim dead-ends keep their accepted-chain semantics:
/// their far endpoints are intersection-created nodes, never original
/// vertices.
pub(super) fn apex_turnaround_candidate(
    graph: &super::super::graph::IntersectionGraph,
    end_node: usize,
    he_idx: usize,
    twin_idx: usize,
    face_hes: &std::collections::HashSet<usize>,
    visited: &std::collections::HashSet<usize>,
) -> Option<usize> {
    let node_is_original = matches!(
        graph.nodes[end_node].source,
        super::super::graph::NodeSource::OriginalA(_)
            | super::super::graph::NodeSource::OriginalB(_)
    );
    let only_twin_departs = graph.nodes[end_node]
        .outgoing
        .iter()
        .all(|&o| o == twin_idx || o == he_idx || !face_hes.contains(&o));
    if node_is_original
        && only_twin_departs
        && face_hes.contains(&twin_idx)
        && !visited.contains(&twin_idx)
    {
        return Some(twin_idx);
    }
    None
}

/// Boundary half-edges of the face (both halves of the edge) that a
/// Shared record registered on the face COVERS: the same two nodes
/// and the same 3-D points at the quarter, middle and three-quarter
/// stations, in either direction. A cut that runs along an edge the
/// operand already has (a box whose bottom lies in a torus's equator
/// plane) gives the face two coincident closed pairs at one node; the
/// 3D-Newell walk chained all four halves into a zero-area loop and
/// lost the equator from BOTH bands (measured on a revolved circle's
/// torus, PK-420). Here the RECORD outranks the boundary edge — the
/// opposite of the chart walker's dedup — because this walker treats
/// a Shared record as a cut (in, around, back out) while a boundary
/// edge whose twin is on the same face is a seam it passes through:
/// keeping the original folded both bands into one slit fragment.
/// The assembly unifies the record with the original's edge, so the
/// served body keeps the original's exact curve.
fn boundary_covered_by_records(
    graph: &IntersectionGraph,
    he_indices: &[usize],
    gate: f64,
) -> Vec<usize> {
    use popkern_core::geometry::curve::CurveEval as _;
    let is_record = |h: usize| {
        matches!(
            graph.half_edges[h].source,
            crate::merge::graph::EdgeSource::Shared { .. }
        )
    };
    let ends = |h: usize| -> (usize, usize) {
        (
            graph.half_edges[h].origin,
            graph.half_edges[graph.half_edges[h].twin].origin,
        )
    };
    // The 3-D point a fraction `f` along the half-edge's own traversal.
    let probe = |h: usize, f: f64| -> Option<popkern_core::geometry::Point3> {
        let he = &graph.half_edges[h];
        if let Some(c) = he.curve.as_ref() {
            return Some(c.geometry.point_at(c.t_start + f * (c.t_end - c.t_start)));
        }
        let c = graph.half_edges[he.twin].curve.as_ref()?;
        Some(c.geometry.point_at(c.t_end - f * (c.t_end - c.t_start)))
    };
    let mut covered = Vec::new();
    for &b in he_indices.iter().filter(|&&h| !is_record(h)) {
        let (b0, b1) = ends(b);
        let by_record = he_indices.iter().any(|&r| {
            if !is_record(r) {
                return false;
            }
            let (r0, r1) = ends(r);
            let same = r0 == b0 && r1 == b1;
            let reversed = r0 == b1 && r1 == b0;
            if !(same || reversed) {
                return false;
            }
            [0.25_f64, 0.5, 0.75].iter().all(|&f| {
                let g = if same { f } else { 1.0 - f };
                match (probe(r, f), probe(b, g)) {
                    (Some(p), Some(q)) => p.distance_to(&q) <= gate,
                    // Curveless straight spans: the nodes decide, but a
                    // closed span (one node) cannot be judged.
                    _ => b0 != b1,
                }
            })
        });
        if by_record {
            covered.push(b);
            covered.push(graph.half_edges[b].twin);
        }
    }
    covered
}

/// Mark every boundary half-edge a record covers as visited, so the
/// 3D-Newell walk neither starts from nor steps onto it.
pub(super) fn skip_covered(
    graph: &IntersectionGraph,
    he_indices: &[usize],
    tolerance: &popkern_core::tolerance::Tolerance,
    visited: &mut HashSet<usize>,
) {
    visited.extend(boundary_covered_by_records(
        graph,
        he_indices,
        tolerance.linear() * 10.0,
    ));
}
