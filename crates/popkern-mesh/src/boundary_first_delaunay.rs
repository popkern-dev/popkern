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
    // Standard incircle determinant with orientation normalization.
    let orient = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    if orient.abs() < 1e-300 {
        return false;
    }
    let (ax, ay) = (a.0 - p.0, a.1 - p.1);
    let (bx, by) = (b.0 - p.0, b.1 - p.1);
    let (cx, cy) = (c.0 - p.0, c.1 - p.1);
    let t1 = (ax * ax + ay * ay) * (bx * cy - by * cx);
    let t2 = (bx * bx + by * by) * (ax * cy - ay * cx);
    let t3 = (cx * cx + cy * cy) * (ax * by - ay * bx);
    let det = t1 - t2 + t3;
    // Relative-epsilon guard: exactly-cocircular points (e.g. hole
    // rings sampled from an analytic circle) give det = ±1e-17 float
    // noise; strict comparison then cascades arbitrary
    // re-triangulations that break boundary segments. Near-zero is
    // OUT of the circumcircle — any triangulation of a cocircular
    // set is equally Delaunay.
    let scale = t1.abs() + t2.abs() + t3.abs();
    let eps = scale * 1.0e-12;
    if orient > 0.0 { det > eps } else { det < -eps }
}
