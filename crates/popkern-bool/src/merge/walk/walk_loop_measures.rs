//! Loop measures for walked fragments: sampled Newell vectors,
//! signed-area proxies, and outer-loop containment votes.
//!
//! Split from `walk_geometry.rs` under the structural line ratchet.

use popkern_core::geometry::Vec3;

use crate::merge::graph::IntersectionGraph;

/// Compute a signed area proxy for a loop using Newell's method.
///
/// Positive = CCW (outer), negative = CW (inner/hole).
pub(super) fn signed_area(graph: &IntersectionGraph, loop_hes: &[usize]) -> f64 {
    let normal = loop_newell(graph, loop_hes);
    normal.length() * normal.z.signum()
}

/// The (borrow-signed) Newell vector of a loop: |N| = 2·area, and
/// the DIRECTION encodes winding — two coplanar loops wind opposite
/// exactly when their Newell vectors point opposite ways. This is
/// the orientation-safe comparison (the historical `sign(N.z)` trick
/// degenerates on vertical faces).
pub(crate) fn loop_newell(graph: &IntersectionGraph, loop_hes: &[usize]) -> Vec3 {
    // Sum cross products of consecutive positions projected onto an
    // approximate normal. CURVED half-edges are SAMPLED along their
    // span — origin-only polygons measure a curved loop as its
    // chord polygon, and an ellipse-fronted face lost the outer-loop
    // election to its own seat ring because the chord quad
    // under-measured the true area (the bulge is ~30% of the face).
    let mut positions: Vec<popkern_core::geometry::Point3> = Vec::new();
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        let origin_pos = graph.nodes[he.origin].position;
        let own = he.curve.is_some();
        let curve = he
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[he.twin].curve.as_ref());
        let curved_curve = curve.filter(|c| match &c.geometry {
            popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
            popkern_core::geometry::curve::CurveGeometry::Nurbs(n) => n.degree > 1,
            _ => true,
        });
        if let Some(c) = curved_curve {
            use popkern_core::geometry::curve::CurveEval;
            let (t0, t1) = (c.t_start, c.t_end);
            // Walk the samples in THIS half-edge's traversal
            // direction: a borrowed (twin) curve runs opposite.
            let n_s = 8;
            for i in 0..n_s {
                let f = i as f64 / n_s as f64;
                let t = if own {
                    t0 + (t1 - t0) * f
                } else {
                    t1 - (t1 - t0) * f
                };
                positions.push(c.geometry.point_at(t));
            }
        } else {
            positions.push(origin_pos);
        }
    }

    // Check if all positions are coincident (self-loop structure).
    // If so, sample curve geometry to compute area.
    let all_coincident = positions.len() >= 2
        && positions
            .iter()
            .all(|p| p.distance_to(&positions[0]) < 1e-10);

    if positions.len() < 3 || all_coincident {
        // Self-loop edge (closed curve like a circle): estimate area from
        // the curve geometry by sampling points along the curve.
        // Try each half-edge in the loop until one has curve geometry.
        for &he_idx in loop_hes {
            let he = &graph.half_edges[he_idx];
            // Curve geometry may live on either member of the pair —
            // a rim half-edge without its own curve measured area 0,
            // lost the outer-loop election to its 12-segment socket
            // hole, and the cap fragment came out inside-out. When
            // the curve is BORROWED from the twin, the samples run
            // opposite this loop's traversal, so the SIGN of the
            // area must flip (hole/piece role assignment reads the
            // sign; an unflipped borrow made cap inner loops wind
            // with their outer and reopened the socket seams).
            let own = he.curve.is_some();
            let curve = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref());
            let borrow_sign = if own { 1.0 } else { -1.0 };
            if let Some(curve) = curve {
                use popkern_core::geometry::curve::CurveEval;
                let n_samples = 16;
                let (t0, t1) = (curve.t_start, curve.t_end);
                let mut pts = Vec::with_capacity(n_samples);
                for i in 0..n_samples {
                    let t = t0 + (t1 - t0) * (i as f64) / (n_samples as f64);
                    pts.push(curve.geometry.point_at(t));
                }
                // Newell's method on sampled points. Scale matches the
                // polygon branch below (|N| = 2·area) so closed rims
                // and polygonal loops compete on one measure.
                let mut nx = 0.0;
                let mut ny = 0.0;
                let mut nz = 0.0;
                for i in 0..pts.len() {
                    let cur = &pts[i];
                    let next = &pts[(i + 1) % pts.len()];
                    nx += (cur.y - next.y) * (cur.z + next.z);
                    ny += (cur.z - next.z) * (cur.x + next.x);
                    nz += (cur.x - next.x) * (cur.y + next.y);
                }
                let normal = Vec3::new(nx, ny, nz);
                return normal * borrow_sign;
            }
        }
        return Vec3::new(0.0, 0.0, 0.0);
    }

    // Newell's method for 3D polygon normal (proportional to area).
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;

    let n = positions.len();
    for i in 0..n {
        let cur = &positions[i];
        let next = &positions[(i + 1) % n];
        nx += (cur.y - next.y) * (cur.z + next.z);
        ny += (cur.z - next.z) * (cur.x + next.x);
        nz += (cur.x - next.x) * (cur.y + next.y);
    }

    // The length of the vector is proportional to twice the area;
    // the direction encodes winding.
    Vec3::new(nx, ny, nz)
}

/// True when `lp`'s nodes are (by strict majority) OUTSIDE the outer
/// loop's sampled polygon — such a loop is a disjoint region, never
/// a hole of this fragment. Points within a small band of the
/// boundary are skipped (loops may share junction nodes with the
/// outer ring). Projection: the face plane's two dominant axes from
/// the outer Newell vector.
pub(super) fn loop_outside_outer(
    graph: &IntersectionGraph,
    lp: &[usize],
    outer: &[usize],
    outer_newell: &Vec3,
) -> bool {
    // 2D chart: drop the dominant normal component.
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
    // Sampled outer polygon (curved half-edges sampled like
    // loop_newell does, via node origins + curve stations).
    let mut ring: Vec<(f64, f64)> = Vec::new();
    for &he_idx in outer {
        let he = &graph.half_edges[he_idx];
        let own = he.curve.is_some();
        let curve = he
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[he.twin].curve.as_ref());
        let curved_curve = curve.filter(|c| match &c.geometry {
            popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
            popkern_core::geometry::curve::CurveGeometry::Nurbs(nc) => nc.degree > 1,
            _ => true,
        });
        if let Some(c) = curved_curve {
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
                ring.push((coord(&q, ax), coord(&q, ay)));
            }
        } else {
            let p = graph.nodes[he.origin].position;
            ring.push((coord(&p, ax), coord(&p, ay)));
        }
    }
    if ring.len() < 3 {
        return false;
    }
    // Scale-aware boundary band.
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(x, y) in &ring {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    let band = 1.0e-3 * ((hi.0 - lo.0) + (hi.1 - lo.1)).max(1e-9);
    let mut inside_n = 0usize;
    let mut outside_n = 0usize;
    for &he_idx in lp {
        let p = graph.nodes[graph.half_edges[he_idx].origin].position;
        let (px, py) = (coord(&p, ax), coord(&p, ay));
        // Distance to the ring: skip near-boundary points.
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
            outside_n += 1;
        }
    }
    // Require an UNAMBIGUOUS verdict: every tested node outside and
    // enough of them to trust. A mixed or thin verdict keeps the
    // existing size-rule behavior (measured: a through-slot's ring
    // on a long tapered side face produced noisy outside votes and
    // the demotion blocked the tunnel — χ' 0 → 2).
    inside_n == 0 && outside_n >= 3
}
