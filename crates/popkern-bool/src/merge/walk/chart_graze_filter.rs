//! Interior-graze and thin-hairpin filters for chart subdivisions
//! (split from `chart_subdivision.rs` under the file-size ratchet).

use popkern_core::tolerance::Tolerance;

use super::chart_subdivision::{FaceChartSubdivision, FilteredRecord};
use crate::merge::graph::IntersectionGraph;

/// Run the thin-hairpin cull and the iterated interior-graze filter
/// over `sub.embeddings` (see the block comments inside).
pub(super) fn filter_hairpins_and_grazes(
    graph: &IntersectionGraph,
    sub: &mut FaceChartSubdivision,
    tolerance: &Tolerance,
    ends: &dyn Fn(usize) -> (usize, usize),
) {
    // THIN HAIRPINS: a same-node pair whose uv run RETRACES itself
    // (the reversed run coincides with the forward run) is a
    // touch-point record folded onto one node — its two directions
    // sit at identical wedge angles and collide the rotational
    // matching (measured: the angled tap's 0.0015-long floor×flank
    // stub at the hole-bottom corner sent the whole flank to the
    // 3D-Newell fallback, PK-465; the m6 nut's mouth-chamfer folds
    // regressed the moment the cull was removed). SCOPED to SMALL
    // pairs — 3D extent under 1e-3 of the node cloud's diagonal:
    // an unscoped retrace cull removed a LOAD-BEARING full-size
    // pair on the keyed-bore fixture and its merge stranded. The
    // RETRACE test (not area) is the fold discriminator: a thin but
    // genuine ring's reversed run walks the far side.
    {
        let diag = {
            let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
            for node in &graph.nodes {
                let q = node.position;
                for (i, c) in [q.x, q.y, q.z].into_iter().enumerate() {
                    lo[i] = lo[i].min(c);
                    hi[i] = hi[i].max(c);
                }
            }
            ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt()
        };
        let small_gate = (diag * 1.0e-3).max(tolerance.linear() * 20.0);
        let hairpins: Vec<usize> = sub
            .embeddings
            .iter()
            .filter(|&(&h, emb)| {
                let (a, b) = ends(h);
                if a != b || emb.uv.len() < 2 {
                    return false;
                }
                let extent3 = {
                    use popkern_core::geometry::curve::CurveEval as _;
                    let he = &graph.half_edges[h];
                    let c = he
                        .curve
                        .as_ref()
                        .or_else(|| graph.half_edges[he.twin].curve.as_ref());
                    match c {
                        Some(c) => {
                            let anchor3 = graph.nodes[a].position;
                            (0..=6)
                                .map(|k| {
                                    let t = c.t_start + (c.t_end - c.t_start) * k as f64 / 6.0;
                                    c.geometry.point_at(t).distance_to(&anchor3)
                                })
                                .fold(0.0_f64, f64::max)
                        }
                        None => 0.0,
                    }
                };
                if extent3 > small_gate {
                    return false;
                }
                let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) =
                    (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
                for &(u, v) in &emb.uv {
                    lo_u = lo_u.min(u);
                    hi_u = hi_u.max(u);
                    lo_v = lo_v.min(v);
                    hi_v = hi_v.max(v);
                }
                let extent = (hi_u - lo_u) + (hi_v - lo_v);
                if extent <= 0.0 {
                    return false;
                }
                let n = emb.uv.len();
                let retrace = (0..n)
                    .map(|i| {
                        let (u1, v1) = emb.uv[i];
                        let (u2, v2) = emb.uv[n - 1 - i];
                        ((u1 - u2).powi(2) + (v1 - v2).powi(2)).sqrt()
                    })
                    .fold(0.0_f64, f64::max);
                retrace < extent * 0.25
            })
            .map(|(&h, _)| h)
            .collect();
        for h in hairpins {
            sub.embeddings.remove(&h);
            sub.report.filtered.push(FilteredRecord::ZeroExtent(h));
        }
    }

    // INTERIOR-GRAZE filter (design doc §5): a record pair whose
    // endpoint node has NO other registered half-edge dead-ends
    // strictly inside the face — both sides of the trace bound the
    // same region there, and walking it weaves an out-and-back spur
    // into a fragment loop (measured on the flush-gusset wall: the
    // gusset's sloped-face trace ends mid-wall at a degree-1 node,
    // and the chart outer loop carried the spur as its two extra
    // half-edges vs the 3D-Newell walk, which filters grazes at
    // candidate time). Iterated: removing one graze can expose a
    // chain.
    for _ in 0..8 {
        let mut incident: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for &h in sub.embeddings.keys() {
            let (a, b) = ends(h);
            incident.entry(a).or_default().push(h);
            incident.entry(b).or_default().push(h);
        }
        let mut grazes: Vec<(usize, usize)> = Vec::new();
        for &h in sub.embeddings.keys() {
            let tw = graph.half_edges[h].twin;
            if h > tw || !sub.embeddings.contains_key(&tw) {
                continue;
            }
            if !matches!(
                graph.half_edges[h].source,
                crate::merge::graph::EdgeSource::Shared { .. }
            ) {
                continue;
            }
            let (a, b) = ends(h);
            if a == b {
                continue; // closed rings are legitimate cuts
            }
            let dead_end = |n: usize| -> bool {
                incident
                    .get(&n)
                    .map(|v| v.iter().all(|&x| x == h || x == tw))
                    .unwrap_or(true)
            };
            if dead_end(a) || dead_end(b) {
                if popkern_core::trace::trace_enabled("chart-dump") {
                    eprintln!(
                        "    [chart] graze-root he {h}: dead a={} (n{a}) dead b={} (n{b})",
                        dead_end(a),
                        dead_end(b)
                    );
                }
                grazes.push((h, tw));
            }
        }
        if grazes.is_empty() {
            break;
        }
        for (h, tw) in grazes {
            if popkern_core::trace::trace_enabled("chart-dump") {
                let (a, b) = ends(h);
                let pa = graph.nodes[a].position;
                let pb = graph.nodes[b].position;
                eprintln!(
                    "    [chart] graze-drop he {h} tw {tw} n{a}({:.4},{:.4},{:.4})->n{b}({:.4},{:.4},{:.4})",
                    pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
                );
            }
            sub.embeddings.remove(&h);
            sub.embeddings.remove(&tw);
            sub.report.filtered.push(FilteredRecord::InteriorGraze(h));
        }
    }
}
