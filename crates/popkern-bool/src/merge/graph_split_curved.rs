//! A point on a curved graph edge's interior, split from
//! `graph_split.rs` under the 500-line file cap.

use popkern_core::geometry::Point3;

use super::super::graph::IntersectionGraph;

/// Whether any knot span of `n` overlapping `[ta, tb]` has `p`
/// inside its LOCAL control hull's AABB inflated by `pad`. Each
/// span's curve piece is a convex combination of its `degree + 1`
/// control points, so a point beyond every span box provably cannot
/// lie within `pad` of the queried curve segment — an EXACT reject
/// for on-curve decisions, at three float compares per span instead
/// of a knot-seeded global inversion. Multi-wrap tool boundary coils
/// carry hundreds of spans, and the splitters' node × edge scans
/// paid the full inversion per pair — O(nodes × spans) at 80+
/// thread turns, the cellular merge stage's residual quadratic term
/// after PK-496 verdict propagation (measured by stack sampling:
/// `nearest_param_on_curve` under `split_edges_at_shared_nodes`).
fn span_hull_within(
    n: &popkern_core::geometry::NurbsCurve,
    ta: f64,
    tb: f64,
    p: &Point3,
    pad: f64,
) -> bool {
    let d = n.degree;
    let (lo_t, hi_t) = if ta <= tb { (ta, tb) } else { (tb, ta) };
    if n.knots.len() < 2 * d + 2 || n.control_points.len() + d + 1 != n.knots.len() {
        return true; // malformed shape: no opinion, run the exact test
    }
    for j in d..(n.knots.len() - 1 - d) {
        let (k0, k1) = (n.knots[j], n.knots[j + 1]);
        if k1 <= k0 || k1 < lo_t || k0 > hi_t {
            continue;
        }
        let mut inside = true;
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for cp in &n.control_points[j - d..=j] {
            for (k, c) in [cp.x, cp.y, cp.z].into_iter().enumerate() {
                lo[k] = lo[k].min(c);
                hi[k] = hi[k].max(c);
            }
        }
        for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
            if c < lo[k] - pad || c > hi[k] + pad {
                inside = false;
                break;
            }
        }
        if inside {
            return true;
        }
    }
    false
}

/// If `point` lies on the REAL curve of graph edge `he_idx` (which
/// must be CURVED — straight edges are fully served by the chord
/// test), return its chord-equivalent fraction. Coarse scan +
/// local ternary refinement; interior hits only (endpoint
/// matches are the caller's node-snap job).
pub(in crate::merge) fn point_on_curved_graph_edge(
    graph: &IntersectionGraph,
    he_idx: usize,
    point: &Point3,
    tolerance: f64,
) -> Option<f64> {
    let he = &graph.half_edges[he_idx];
    let curve = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())?;
    let curved = match &curve.geometry {
        popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
        popkern_core::geometry::curve::CurveGeometry::Nurbs(n) => n.degree > 1,
        _ => true,
    };
    if !curved {
        return None;
    }
    // Knot-seeded global inversion: the private 16-sample scan here
    // ranked the ADJACENT COIL TURN's basin above the true one on a
    // 6.5-turn helical seam (0.209 vs 0.219), so the seam × plate
    // pierce point classified OnFaceInterior, the seam never split at
    // the crossing, and the whole wall dropped to a knife-edge sample
    // (PK-449 — the station-aligned deep grind's real mechanism).
    let (ta, tb) = (curve.t_start, curve.t_end);
    // Span-hull pre-gate (PK-496): a point beyond every overlapping
    // span's local control hull cannot be an interior on-curve hit —
    // reject before paying the global inversion. A hit outside the
    // queried [ta, tb] sub-range would clamp to an endpoint fraction
    // and be rejected below anyway, so scoping the gate to
    // overlapping spans preserves the answer exactly.
    if let popkern_core::geometry::curve::CurveGeometry::Nurbs(n) = &curve.geometry
        && !span_hull_within(n, ta, tb, point, tolerance)
    {
        return None;
    }
    let (best_t, best_d) =
        popkern_core::geometry::curve_inversion::nearest_param_on_curve(curve, point);
    if best_d > tolerance {
        return None;
    }
    let f = ((best_t - ta) / (tb - ta)).clamp(0.0, 1.0);
    // Interior only: endpoints belong to node snapping. Clear of the
    // ends by the tolerance as a length, never a parameter fraction
    // (PK-1131): a ball's seam crossing of another ball 2.1e-6 up the
    // seam from its vertex on a box's top lay a millionth of the seam
    // piece's parameter from it, inside or outside by the placement,
    // and where it went unsplit the circle the balls overlap in started
    // at the vertex (PK-1181).
    let ends = [
        graph.nodes[he.origin].position,
        graph.nodes[graph.half_edges[he.twin].origin].position,
    ];
    if f <= 0.0 || f >= 1.0 || ends.iter().any(|e| e.distance_to(point) <= tolerance) {
        return None;
    }
    Some(f)
}
