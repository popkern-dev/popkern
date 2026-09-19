//! Bowyer–Watson Delaunay triangulation over 2-D sample points —
//! split from `boundary_first.rs` for the 500-line/frozen-ratchet
//! structural gate.

use std::collections::HashMap;

/// Bowyer–Watson Delaunay triangulation over a point set. Returns
/// `None` on degenerate input (all collinear etc.).
pub(super) fn bowyer_watson(points: &[(f64, f64)]) -> Option<Vec<[usize; 3]>> {
    if points.len() < 3 {
        return None;
    }
    // Super-triangle enclosing everything.
    let (mut x0, mut x1, mut y0, mut y1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(x, y) in points {
        x0 = x0.min(x);
        x1 = x1.max(x);
        y0 = y0.min(y);
        y1 = y1.max(y);
    }
    let d = ((x1 - x0).max(y1 - y0)).max(1e-9) * 16.0;
    let cx = 0.5 * (x0 + x1);
    let cy = 0.5 * (y0 + y1);
    let mut pts: Vec<(f64, f64)> = points.to_vec();
    let s0 = pts.len();
    pts.push((cx - d, cy - d * 0.6));
    pts.push((cx + d, cy - d * 0.6));
    pts.push((cx, cy + d));

    let mut tris: Vec<[usize; 3]> = vec![[s0, s0 + 1, s0 + 2]];
    for (pi, &p) in points.iter().enumerate() {
        // Find triangles whose circumcircle contains p.
        let mut bad: Vec<usize> = Vec::new();
        for (ti, t) in tris.iter().enumerate() {
            if in_circumcircle(p, pts[t[0]], pts[t[1]], pts[t[2]]) {
                bad.push(ti);
            }
        }
        // Boundary of the cavity: edges appearing exactly once.
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        let mut edge_dir: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        for &ti in &bad {
            let t = tris[ti];
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let key = (a.min(b), a.max(b));
                *edge_count.entry(key).or_insert(0) += 1;
                edge_dir.insert(key, (a, b));
            }
        }
        // Remove bad triangles (descending index).
        bad.sort_unstable_by(|a, b| b.cmp(a));
        for ti in bad {
            tris.swap_remove(ti);
        }
        // Deterministic cavity fill: HashMap iteration order varies
        // with the hasher's per-thread seed, so filling in map order
        // reordered the emitted triangles thread-to-thread and
        // run-to-run (caught by the PK-277 thread-count identity
        // pin). The triangle SET is order-independent; the emission
        // order must be too.
        let mut cavity_boundary: Vec<(usize, usize)> = edge_count
            .iter()
            .filter(|&(_, &count)| count == 1)
            .map(|(&key, _)| key)
            .collect();
        cavity_boundary.sort_unstable();
        for key in cavity_boundary {
            let (a, b) = edge_dir[&key];
            tris.push([a, b, pi]);
        }
    }
    tris.retain(|t| t.iter().all(|&i| i < s0));
    if tris.is_empty() { None } else { Some(tris) }
}

fn in_circumcircle(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    // EXACT adaptive predicates (Shewchuk, via the `robust` crate) —
    // PK-494. The float determinant with a relative-epsilon guard
    // degraded on multi-wrap wall charts: metric scaling puts the
    // chart at ~776×0.43 with thousands of exactly-collinear
    // chart-edge samples, the super-triangle multiplies coordinate
    // magnitudes by 16×extent, and sliver circumcircle determinants
    // fell inside the epsilon — cavity boundaries went non-manifold,
    // consecutive crease-chain segments vanished from the
    // triangulation, and the coverage gate sent every tool wall face
    // to the uniform-grid fallback (the 41-tool-turn torn serve).
    // Exactly-cocircular stays OUT of the circumcircle, now exactly:
    // any triangulation of a cocircular set is equally Delaunay, and
    // a hole ring sampled from an analytic circle must not cascade
    // arbitrary re-triangulations.
    let orient = robust::orient2d(
        robust::Coord { x: a.0, y: a.1 },
        robust::Coord { x: b.0, y: b.1 },
        robust::Coord { x: c.0, y: c.1 },
    );
    if orient == 0.0 {
        // Degenerate triangle: no circumcircle.
        return false;
    }
    let det = robust::incircle(
        robust::Coord { x: a.0, y: a.1 },
        robust::Coord { x: b.0, y: b.1 },
        robust::Coord { x: c.0, y: c.1 },
        robust::Coord { x: p.0, y: p.1 },
    );
    if orient > 0.0 { det > 0.0 } else { det < 0.0 }
}
