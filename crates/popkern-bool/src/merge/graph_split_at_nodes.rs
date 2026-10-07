//! Splitting shared (intersection-curve) graph edges at nodes that lie
//! on them, with a conservative bounds prefilter (PK-708).

use popkern_core::geometry::Point3;

use super::graph::{EdgeSource, IntersectionGraph};
use super::graph_split::point_on_segment;

/// A box that certainly contains shared half-edge `i`: the edge's
/// chord for an edge without a curve, the span's ends for a line, the
/// whole circle or ellipse for conics, the control hull for a NURBS
/// curve; `None` (test everything) for curves without a cheap bound.
fn shared_edge_bounds(graph: &IntersectionGraph, i: usize) -> Option<(Point3, Point3)> {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
    let he = &graph.half_edges[i];
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut take = |p: Point3| {
        lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    };
    let sphere = |c: Point3, r: f64| {
        (
            Point3::new(c.x - r, c.y - r, c.z - r),
            Point3::new(c.x + r, c.y + r, c.z + r),
        )
    };
    match &he.curve {
        None => {
            take(graph.nodes[he.origin].position);
            take(graph.nodes[graph.half_edges[he.twin].origin].position);
        }
        Some(curve) => match &curve.geometry {
            CurveGeometry::Line { .. } => {
                take(curve.geometry.point_at(curve.t_start));
                take(curve.geometry.point_at(curve.t_end));
            }
            CurveGeometry::Circle { center, radius, .. }
            | CurveGeometry::Arc { center, radius, .. } => return Some(sphere(*center, *radius)),
            CurveGeometry::Ellipse {
                center,
                u_basis,
                v_basis,
                ..
            } => return Some(sphere(*center, u_basis.length() + v_basis.length())),
            CurveGeometry::Nurbs(n) => n.control_points.iter().copied().for_each(&mut take),
            _ => return None,
        },
    }
    lo.x.is_finite().then_some((lo, hi))
}

/// Split a graph edge at an existing node (no new node creation).
/// Split SHARED (intersection-curve) edges at the given operand
/// nodes: an edge-face hit puts a node on an operand edge, and any
/// intersection curve passing THROUGH that node (a boss seam
/// crossing a boss/plate intersection circle) must be split there or
/// the loop walk cannot route through the junction. This is the
/// inverse of [`split_edges_at_shared_nodes`].
pub(crate) fn split_shared_edges_at_nodes(
    graph: &mut IntersectionGraph,
    nodes: &[usize],
    tolerance: f64,
) {
    // Conservative bounds per shared half-edge (PK-708): a node outside
    // an edge's bounds grown by the tolerance can never test on its
    // curve, so the expensive curve test runs only where it can match —
    // the same decisions in the same order, so the same splits. Bounds
    // cover the whole curve (or, for a line, the edge's span), so a
    // cached box stays conservative after a split shrinks its edge.
    // Without them every node met every shared edge on every restart
    // (a 60-knuckle piano hinge spent 1.0 s of 1.1 s here).
    let mut bounds: std::collections::HashMap<usize, Option<(Point3, Point3)>> =
        std::collections::HashMap::new();
    let mut did_split = true;
    while did_split {
        did_split = false;
        for &node_idx in nodes {
            let pos = graph.nodes[node_idx].position;
            let mut i = 0;
            while i < graph.half_edges.len() {
                let he = &graph.half_edges[i];
                if i >= he.twin || !matches!(he.source, EdgeSource::Shared { .. }) {
                    i += 2;
                    continue;
                }
                let twin_origin = graph.half_edges[he.twin].origin;
                if he.origin == node_idx || twin_origin == node_idx {
                    i += 2;
                    continue;
                }
                let bound = *bounds
                    .entry(i)
                    .or_insert_with(|| shared_edge_bounds(graph, i));
                if let Some((lo, hi)) = bound
                    && (pos.x < lo.x - tolerance
                        || pos.y < lo.y - tolerance
                        || pos.z < lo.z - tolerance
                        || pos.x > hi.x + tolerance
                        || pos.y > hi.y + tolerance
                        || pos.z > hi.z + tolerance)
                {
                    i += 2;
                    continue;
                }
                let he = &graph.half_edges[i];
                // Distance from the node to the edge's CURVE (shared
                // edges carry true intersection geometry — circles
                // included; the chord test would miss them).
                let on_curve = if let Some(curve) = &he.curve {
                    let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
                    let (t0, t1) = (curve.t_start, curve.t_end);
                    let mut best = f64::INFINITY;
                    let mut best_t = t0;
                    for k in 0..=64 {
                        let t = t0 + (t1 - t0) * k as f64 / 64.0;
                        let d = eval.point_at(t).distance_to(&pos);
                        if d < best {
                            best = d;
                            best_t = t;
                        }
                    }
                    let (mut lo, mut hi) = (
                        (best_t - (t1 - t0) / 64.0).max(t0.min(t1)),
                        (best_t + (t1 - t0) / 64.0).min(t0.max(t1)),
                    );
                    for _ in 0..40 {
                        let m1 = lo + (hi - lo) / 3.0;
                        let m2 = hi - (hi - lo) / 3.0;
                        if eval.point_at(m1).distance_to(&pos) < eval.point_at(m2).distance_to(&pos)
                        {
                            hi = m2;
                        } else {
                            lo = m1;
                        }
                    }
                    let t = 0.5 * (lo + hi);
                    // Interior only — endpoint touches are already
                    // nodes. POSITIONAL, not parametric: a record's
                    // genuine crossing can sit arbitrarily close to
                    // its end in parameter while being a distinct 3D
                    // junction (measured, PK-491: the m5 tip-cone
                    // groove wall crosses the cone seam 1.2e-3 rad
                    // before its endpoint — relative parameter
                    // 8.5e-4 — and the old 1e-3 parametric exclusion
                    // vetoed the split; the seam-hit junction never
                    // formed, the unsplit seam edge self-overlapped
                    // the chart polygon, and the tessellator's
                    // grid fallback tore the served tip band).
                    let p_split = eval.point_at(t);
                    let start_pos = graph.nodes[he.origin].position;
                    let end_pos = graph.nodes[twin_origin].position;
                    let interior = p_split.distance_to(&start_pos) > tolerance
                        && p_split.distance_to(&end_pos) > tolerance;
                    interior && p_split.distance_to(&pos) < tolerance
                } else {
                    let p0 = graph.nodes[he.origin].position;
                    let p1 = graph.nodes[twin_origin].position;
                    point_on_segment(&pos, &p0, &p1, tolerance).is_some()
                };
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [shared-split] node {node_idx} vs shared he {i}: on_curve={on_curve}"
                    );
                }
                if on_curve {
                    // A REFUSED split (window-boundary minimizer)
                    // must not count as progress: at coarse
                    // tolerance the same (node, piece) re-matches
                    // every restart and the loop never terminates
                    // (measured, PK-342 envelope at linear 1e-2:
                    // two hours spinning here).
                    did_split = super::graph_split_node::split_at_existing_node(graph, i, node_idx);
                    if !did_split {
                        i += 2;
                        continue;
                    }
                    break;
                }
                i += 2;
            }
            if did_split {
                break;
            }
        }
    }
}
