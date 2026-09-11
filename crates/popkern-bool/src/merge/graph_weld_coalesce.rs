//! Split-junction coalescing sweep (split from `graph_weld.rs`
//! under the 500-line structural budget). See the rationale comment
//! on [`split_junction_coalesce`].

use super::super::graph::{EdgeSource, IntersectionGraph};

/// SPLIT-JUNCTION COALESCING (see `graph_weld.rs` for the tier
/// layout): collect welds of record-lighter junction siblings onto
/// their directly-adjacent boundary junction.
pub(super) fn split_junction_coalesce(
    graph: &IntersectionGraph,
    incident: &[Vec<usize>],
    band: f64,
    welds: &mut Vec<(usize, usize)>,
) {
    let n_nodes = graph.nodes.len();
    let is_shared = |h: usize| matches!(graph.half_edges[h].source, EdgeSource::Shared { .. });
    // SPLIT-JUNCTION COALESCING: a record endpoint landing fit-slop
    // off its true boundary junction can also SPLIT the boundary edge
    // at its own sloppy position — the graph then carries TWO
    // junction nodes a sliver apart on one boundary edge, and on the
    // OTHER face's chart (where the boundary edge does not exist)
    // the ring gaps and the graze cascade unwalks it (measured,
    // PK-464: the chamfered nut's crest arc minted a junction 7.1e-3
    // from the flank's exact one). Weld the record-lighter node onto
    // its DIRECTLY-ADJACENT junction sibling when the record's
    // tangent ray points at it; the boundary sliver between them
    // collapses to zero length for the degenerate cleanup.
    for (n, inc_n) in incident.iter().enumerate().take(n_nodes) {
        if inc_n.is_empty() {
            continue;
        }
        let records_n: Vec<usize> = inc_n.iter().copied().filter(|&h| is_shared(h)).collect();
        let boundary_n: Vec<usize> = inc_n.iter().copied().filter(|&h| !is_shared(h)).collect();
        if records_n.len() != 1 || boundary_n.is_empty() {
            continue;
        }
        let p = graph.nodes[n].position;
        // Directly-adjacent nodes along boundary half-edges.
        for &bh in &boundary_n {
            let m = graph.half_edges[graph.half_edges[bh].twin].origin;
            if m == n {
                continue;
            }
            let d = graph.nodes[m].position.distance_to(&p);
            // Never weld across an existing record connection: two
            // junctions JOINED by a record piece are distinct by
            // construction (the tuned exact-connector class — the
            // m6 nut's mouth-chamfer chain carries real 1.5e-3
            // connectors between junctions 1.0e-3 apart; welding
            // their ends collapsed the connector into a self-loop,
            // the stub cleanup retired it, and the chain gapped into
            // a graze cascade). A weld target must be an UNCONNECTED
            // duplicate of this node.
            let connected = incident[n]
                .iter()
                .any(|&h| graph.half_edges[graph.half_edges[h].twin].origin == m);
            if connected {
                continue;
            }
            let m_has_record = graph
                .half_edges
                .iter()
                .any(|hd| hd.origin == m && matches!(hd.source, EdgeSource::Shared { .. }));
            // Verified by the RECORD's own tangent ray — the boundary
            // hes trivially point along the boundary at any sibling.
            let record_points_at_m = {
                use popkern_core::geometry::curve::CurveEval as _;
                let rh = records_n[0];
                let he = &graph.half_edges[rh];
                let q = graph.nodes[graph.half_edges[he.twin].origin].position;
                let tangent_from = |outward: popkern_core::geometry::Vec3| {
                    let to_m = graph.nodes[m].position - p;
                    let along = to_m.dot(&outward);
                    let perp = (to_m - outward * along).length();
                    along >= -band && perp < band
                };
                match he
                    .curve
                    .as_ref()
                    .or_else(|| graph.half_edges[he.twin].curve.as_ref())
                {
                    Some(c) => {
                        let d_start = c.geometry.point_at(c.t_start).distance_to(&p);
                        let d_end = c.geometry.point_at(c.t_end).distance_to(&p);
                        let span = c.t_end - c.t_start;
                        let (t_at, t_in) = if d_start <= d_end {
                            (c.t_start, c.t_start + span * 1e-3)
                        } else {
                            (c.t_end, c.t_end - span * 1e-3)
                        };
                        let dvec = c.geometry.point_at(t_at) - c.geometry.point_at(t_in);
                        dvec.normalized().map(tangent_from).unwrap_or(false)
                    }
                    None => (p - q).normalized().map(tangent_from).unwrap_or(false),
                }
            };
            if d > band * 100.0 || !m_has_record || !record_points_at_m {
                continue;
            }
            if popkern_core::trace::trace_enabled("boolean") {
                let q = graph.nodes[m].position;
                eprintln!(
                    "  [weld-junction] node {n} ({:.6},{:.6},{:.6}) -> {m} ({:.6},{:.6},{:.6}) d={d:.2e}",
                    p.x, p.y, p.z, q.x, q.y, q.z
                );
            }
            welds.push((n, m));
            break;
        }
    }
}
