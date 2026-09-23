//! Sampling helpers shared by the cellular fragment classifier
//! (PK-496 split from `fragment_classify_one.rs` under the 500-line
//! budget): edge-proximity tests, boundary chord midpoints, and the
//! tuned edge-midpoint fragment sample shortcut.

use popkern_core::topology::BRep;

use crate::merge::graph::IntersectionGraph;

pub(super) fn point_near_any_edge(
    brep: &BRep,
    p: &popkern_core::geometry::Point3,
    band: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    for (_, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let a = brep.vertices[he.origin].position;
        let b = brep
            .half_edges
            .get(he.twin)
            .map(|t| brep.vertices[t.origin].position)
            .unwrap_or(a);
        if let Some(c) = edge.curve.as_ref() {
            let mut prev = c.geometry.point_at(c.t_start);
            for i in 1..=16 {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 16.0;
                let q = c.geometry.point_at(t);
                if seg_dist(&prev, &q, p) < band {
                    return true;
                }
                prev = q;
            }
        } else if seg_dist(&a, &b, p) < band {
            return true;
        }
    }
    false
}

pub(super) fn seg_dist(
    a: &popkern_core::geometry::Point3,
    b: &popkern_core::geometry::Point3,
    p: &popkern_core::geometry::Point3,
) -> f64 {
    let ab = *b - *a;
    let l2 = ab.dot(&ab);
    if l2 < 1e-24 {
        return p.distance_to(a);
    }
    let t = ((*p - *a).dot(&ab) / l2).clamp(0.0, 1.0);
    (*a + ab * t).distance_to(p)
}

/// Chord midpoints of the loop's non-Shared half-edges, in loop
/// order, capped at `cap`. The first entry equals the tuned
/// edge-midpoint shortcut, so a unanimous vote reproduces the uncalibrated shortcut
/// sampling exactly.
/// Midpoint of a graph half-edge, ON THE EDGE'S CURVE when it has
/// one (chord midpoint otherwise). A curved edge's chord midpoint is
/// not a point of the fragment at all — a shallow rim-arc sliver cut
/// off by a breached wall has its arc CHORD lying exactly IN the
/// cutting plane, so the chord midpoint sampled the parity coin-flip
/// point on the other body's surface while the arc's belly sits
/// decisively outside it (measured: the sliver classified Inside and
/// a subtract kept a flap beyond the body). Straight edges are
/// unchanged (chord == curve).
pub(super) fn he_curve_midpoint(
    graph: &IntersectionGraph,
    he_idx: usize,
) -> popkern_core::geometry::Point3 {
    use popkern_core::geometry::curve::CurveEval as _;
    let he = &graph.half_edges[he_idx];
    if let Some(c) = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
    {
        return c.geometry.point_at(0.5 * (c.t_start + c.t_end));
    }
    let a = graph.nodes[he.origin].position;
    let b = graph.nodes[graph.half_edges[he.twin].origin].position;
    popkern_core::geometry::Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}

pub(super) fn non_shared_chord_midpoints(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
    cap: usize,
) -> Vec<popkern_core::geometry::Point3> {
    use crate::merge::graph::EdgeSource;
    let mut out = Vec::new();
    for &he_idx in loop_hes {
        if out.len() >= cap {
            break;
        }
        let he = &graph.half_edges[he_idx];
        if matches!(he.source, EdgeSource::Shared { .. }) {
            continue;
        }
        out.push(he_curve_midpoint(graph, he_idx));
    }
    out
}

pub(super) fn sample_fragment_point(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
) -> popkern_core::geometry::Point3 {
    use crate::merge::graph::EdgeSource;

    // Prefer midpoint of a non-shared edge (away from intersection curve)
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        match &he.source {
            EdgeSource::OriginalA(_)
            | EdgeSource::OriginalB(_)
            | EdgeSource::SplitA { .. }
            | EdgeSource::SplitB { .. } => {
                return he_curve_midpoint(graph, he_idx);
            }
            EdgeSource::Shared { .. } => continue,
        }
    }

    // Fallback: centroid of all loop vertices
    if loop_hes.is_empty() {
        return popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for &he_idx in loop_hes {
        let p = graph.nodes[graph.half_edges[he_idx].origin].position;
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    let n = loop_hes.len() as f64;
    popkern_core::geometry::Point3::new(sx / n, sy / n, sz / n)
}
