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
        if curved == Some(true) {
            use popkern_core::geometry::curve::CurveEval;
            let c = he.curve.as_ref().expect("curved implies curve");
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
