//! Hole → region nesting for the legacy fragment election.
//!
//! The flat assignment hung every hole on the primary fragment; with
//! CONCENTRIC rings (countersunk bore: the stale top rim sits inside
//! the new countersink circle) the annulus piece lost its hole and
//! the primary kept a stale one, leaving a single-use edge that made
//! the cellular assembly non-manifold — the silent graph-merge
//! fallback then shipped an invalid body. A hole moves to a separate
//! region only when that region strictly contains it (positive,
//! curve-sampled proof), winds WITH the outer, and shares no graph
//! edge with it — the twin-sharing exclusion that keeps the tuned
//! paired-ring tool semantics (same rule as assemble_chart_roles).

use std::collections::HashSet;

use popkern_core::geometry::Vec3;

use super::super::graph::IntersectionGraph;
use super::walk_geometry::{self, signed_area};

/// Reassign the election's holes to their smallest strictly
/// containing separate region. `inner_loops_out` keeps the holes that
/// stay on the primary fragment; the returned vector is indexed like
/// `separate_outers` and carries each region's holes.
pub(super) fn assign_holes_to_regions(
    graph: &IntersectionGraph,
    outer_newell: &Vec3,
    inner_loops_out: &mut Vec<Vec<usize>>,
    separate_outers: &[Vec<usize>],
) -> Vec<Vec<Vec<usize>>> {
    let mut separate_inners: Vec<Vec<Vec<usize>>> = vec![Vec::new(); separate_outers.len()];
    let mut kept: Vec<Vec<usize>> = Vec::new();
    for hole in inner_loops_out.drain(..) {
        let hole_twins: HashSet<usize> = hole
            .iter()
            .map(|&he_idx| graph.half_edges[he_idx].twin)
            .collect();
        let mut best: Option<(f64, usize)> = None;
        for (si, sep) in separate_outers.iter().enumerate() {
            if sep.iter().any(|h| hole_twins.contains(h)) {
                continue;
            }
            let sep_newell = walk_geometry::loop_newell(graph, sep);
            if sep_newell.dot(outer_newell) <= 0.0 {
                continue;
            }
            // Positive-proof containment (curve-sampled, so 1-node
            // closed rims vote with real points).
            if !loop_strictly_inside_outer(graph, &hole, sep, &sep_newell) {
                continue;
            }
            let sep_area = signed_area(graph, sep).abs();
            if best.is_none_or(|(ba, _)| sep_area < ba) {
                best = Some((sep_area, si));
            }
        }
        match best {
            Some((_, si)) => separate_inners[si].push(hole),
            None => kept.push(hole),
        }
    }
    *inner_loops_out = kept;
    separate_inners
}

/// True when `lp` is UNAMBIGUOUSLY inside the candidate outer ring:
/// every sampled point inside and at least three of them. Closed
/// single-half-edge rims are sampled along their CURVE (they carry
/// one node, and the node-quorum of [`loop_outside_outer`] can never
/// clear on them — a countersunk bore's stale rim hole failed the
/// outside test against a faraway notch piece and was reassigned to
/// it, deleting the bore from the face). Used by the nesting pass:
/// holes move fragments only on positive proof.
fn loop_strictly_inside_outer(
    graph: &IntersectionGraph,
    lp: &[usize],
    outer: &[usize],
    outer_newell: &Vec3,
) -> bool {
    let n = *outer_newell;
    let (ax, ay) = if n.x.abs() >= n.y.abs() && n.x.abs() >= n.z.abs() {
        (1, 2)
    } else if n.y.abs() >= n.z.abs() {
        (0, 2)
    } else {
        (0, 1)
    };
    let coord = |p: &popkern_core::geometry::Point3, a: usize| match a {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    };
    let sample_loop = |hes: &[usize]| -> Vec<(f64, f64)> {
        let mut pts = Vec::new();
        for &he_idx in hes {
            let he = &graph.half_edges[he_idx];
            let own = he.curve.is_some();
            let curve = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref());
            let curved = curve.filter(|c| {
                !matches!(
                    &c.geometry,
                    popkern_core::geometry::curve::CurveGeometry::Line { .. }
                )
            });
            if let Some(c) = curved {
                use popkern_core::geometry::curve::CurveEval;
                let (t0, t1) = (c.t_start, c.t_end);
                for i in 0..8 {
                    let f = i as f64 / 8.0;
                    let t = if own {
                        t0 + (t1 - t0) * f
                    } else {
                        t1 - (t1 - t0) * f
                    };
                    let q = c.geometry.point_at(t);
                    pts.push((coord(&q, ax), coord(&q, ay)));
                }
            } else {
                let p = graph.nodes[he.origin].position;
                pts.push((coord(&p, ax), coord(&p, ay)));
            }
        }
        pts
    };
    let ring = sample_loop(outer);
    if ring.len() < 3 {
        return false;
    }
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(x, y) in &ring {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    let band = 1.0e-3 * ((hi.0 - lo.0) + (hi.1 - lo.1)).max(1e-9);
    let mut inside_n = 0usize;
    for (px, py) in sample_loop(lp) {
        let mut min_d = f64::MAX;
        for i in 0..ring.len() {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % ring.len()];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((px - x0) * dx + (py - y0) * dy) / l2).clamp(0.0, 1.0);
            let d = ((px - x0 - dx * t).powi(2) + (py - y0 - dy * t).powi(2)).sqrt();
            min_d = min_d.min(d);
        }
        if min_d < band {
            continue;
        }
        let mut crossings = 0;
        for i in 0..ring.len() {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % ring.len()];
            if (y0 > py) != (y1 > py) {
                let xc = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
                if xc > px {
                    crossings += 1;
                }
            }
        }
        if crossings % 2 == 1 {
            inside_n += 1;
        } else {
            return false; // any outside vote = not contained
        }
    }
    inside_n >= 3
}
