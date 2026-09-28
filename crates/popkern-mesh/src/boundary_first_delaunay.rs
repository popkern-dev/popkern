//! Bowyer–Watson Delaunay triangulation over 2-D sample points —
//! split from `boundary_first.rs` for the 500-line/frozen-ratchet
//! structural gate.

use std::collections::HashMap;

/// Bowyer–Watson Delaunay triangulation over a point set. Returns
/// `None` on degenerate input (all collinear etc.).
///
/// Each insertion's cavity (every triangle whose circumcircle strictly
/// contains the point) is found by walking to the triangle containing
/// the point and growing across shared edges, not by scanning every
/// triangle: the scan made the whole triangulation O(n²), and at a 6°
/// facet bound an 8k-point coil wall spent 9 s of its 10 s face build
/// here, once per encroachment round (PK-589). The conflict region of
/// a Delaunay triangulation is connected and contains the triangle
/// holding the point, so the grown cavity is the scanned one. Removal
/// and emission follow the scan's exact order, so the output is
/// identical. Whenever the walk or the seed cannot decide (a degenerate
/// triangle, a point on a vertex, a walk that leaves the hull or does
/// not settle), the insertion falls back to the full scan.
pub(super) fn bowyer_watson(points: &[(f64, f64)]) -> Option<Vec<[usize; 3]>> {
    triangulate(points, true)
}

/// Bowyer–Watson with the walked cavity (`walk`) or the reference scan.
fn triangulate(points: &[(f64, f64)], walk: bool) -> Option<Vec<[usize; 3]>> {
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

    let mut mesh = Adjacency::default();
    mesh.push([s0, s0 + 1, s0 + 2]);
    let mut hint = 0usize;
    for (pi, &p) in points.iter().enumerate() {
        let mut bad = walk
            .then(|| mesh.conflict_region(&pts, p, hint))
            .flatten()
            .unwrap_or_else(|| mesh.conflict_scan(&pts, p));
        // Boundary of the cavity: edges appearing exactly once.
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        let mut edge_dir: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        for &ti in &bad {
            let t = mesh.tris[ti];
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let key = (a.min(b), a.max(b));
                *edge_count.entry(key).or_insert(0) += 1;
                edge_dir.insert(key, (a, b));
            }
        }
        // Remove bad triangles (descending index).
        bad.sort_unstable_by(|a, b| b.cmp(a));
        for ti in bad {
            mesh.swap_remove(ti);
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
            mesh.push([a, b, pi]);
        }
        hint = mesh.tris.len().saturating_sub(1);
    }
    let mut tris = mesh.tris;
    tris.retain(|t| t.iter().all(|&i| i < s0));
    if tris.is_empty() { None } else { Some(tris) }
}

/// Unset neighbour slot.
const NONE: usize = usize::MAX;

/// Triangles plus an edge → incident-triangles map, kept current under
/// `push` and `swap_remove`.
#[derive(Default)]
struct Adjacency {
    tris: Vec<[usize; 3]>,
    edges: HashMap<(usize, usize), [usize; 2]>,
}

fn edge_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

impl Adjacency {
    fn link(&mut self, a: usize, b: usize, ti: usize) {
        let slot = self.edges.entry(edge_key(a, b)).or_insert([NONE, NONE]);
        if slot[0] == NONE {
            slot[0] = ti;
        } else {
            slot[1] = ti;
        }
    }

    fn unlink(&mut self, a: usize, b: usize, ti: usize) {
        let key = edge_key(a, b);
        if let Some(slot) = self.edges.get_mut(&key) {
            if slot[0] == ti {
                slot[0] = slot[1];
                slot[1] = NONE;
            } else if slot[1] == ti {
                slot[1] = NONE;
            }
            if slot[0] == NONE {
                self.edges.remove(&key);
            }
        }
    }

    fn relink(&mut self, a: usize, b: usize, from: usize, to: usize) {
        if let Some(slot) = self.edges.get_mut(&edge_key(a, b)) {
            for s in slot.iter_mut() {
                if *s == from {
                    *s = to;
                }
            }
        }
    }

    fn push(&mut self, t: [usize; 3]) {
        let ti = self.tris.len();
        self.tris.push(t);
        for k in 0..3 {
            self.link(t[k], t[(k + 1) % 3], ti);
        }
    }

    /// `Vec::swap_remove` with the edge map following the moved
    /// triangle.
    fn swap_remove(&mut self, ti: usize) {
        let t = self.tris[ti];
        for k in 0..3 {
            self.unlink(t[k], t[(k + 1) % 3], ti);
        }
        let last = self.tris.len() - 1;
        if ti != last {
            let m = self.tris[last];
            for k in 0..3 {
                self.relink(m[k], m[(k + 1) % 3], last, ti);
            }
        }
        self.tris.swap_remove(ti);
    }

    /// The triangle across edge `(a, b)` from `ti`, if any.
    fn across(&self, a: usize, b: usize, ti: usize) -> Option<usize> {
        let slot = self.edges.get(&edge_key(a, b))?;
        let other = if slot[0] == ti { slot[1] } else { slot[0] };
        (other != NONE && other != ti).then_some(other)
    }

    /// Every triangle whose circumcircle strictly contains `p` (the
    /// reference definition).
    fn conflict_scan(&self, pts: &[(f64, f64)], p: (f64, f64)) -> Vec<usize> {
        let mut bad = Vec::new();
        for (ti, t) in self.tris.iter().enumerate() {
            if in_circumcircle(p, pts[t[0]], pts[t[1]], pts[t[2]]) {
                bad.push(ti);
            }
        }
        bad
    }

    /// The conflict region grown from the triangle containing `p`,
    /// located by a visibility walk from `hint`; `None` when the walk
    /// or the seed cannot decide (the caller then scans).
    fn conflict_region(
        &self,
        pts: &[(f64, f64)],
        p: (f64, f64),
        hint: usize,
    ) -> Option<Vec<usize>> {
        let coord = |i: usize| robust::Coord {
            x: pts[i].0,
            y: pts[i].1,
        };
        let q = robust::Coord { x: p.0, y: p.1 };
        let mut ti = hint.min(self.tris.len().checked_sub(1)?);
        let mut settled = false;
        for step in 0..self.tris.len() + 3 {
            let t = self.tris[ti];
            let mut moved = false;
            for r in 0..3 {
                // Rotate the first edge tried so a walk cannot cycle
                // between two triangles on a tie.
                let k = (r + step) % 3;
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let own = robust::orient2d(coord(a), coord(b), coord(c));
                if own == 0.0 {
                    return None;
                }
                let side = robust::orient2d(coord(a), coord(b), q);
                if side != 0.0 && (side > 0.0) != (own > 0.0) {
                    ti = self.across(a, b, ti)?;
                    moved = true;
                    break;
                }
            }
            if !moved {
                settled = true;
                break;
            }
        }
        if !settled {
            return None;
        }
        let seed = self.tris[ti];
        if !in_circumcircle(p, pts[seed[0]], pts[seed[1]], pts[seed[2]]) {
            return None;
        }
        let mut bad = vec![ti];
        let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::from([ti]);
        let mut next = 0;
        while next < bad.len() {
            let t = self.tris[bad[next]];
            let here = bad[next];
            next += 1;
            for k in 0..3 {
                let Some(n) = self.across(t[k], t[(k + 1) % 3], here) else {
                    continue;
                };
                if seen.insert(n) {
                    let m = self.tris[n];
                    if in_circumcircle(p, pts[m[0]], pts[m[1]], pts[m[2]]) {
                        bad.push(n);
                    }
                }
            }
        }
        Some(bad)
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic scatter (xorshift), no external RNG.
    fn scatter(n: usize, seed: u64) -> Vec<(f64, f64)> {
        let mut x = seed;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n).map(|_| (next(), next())).collect()
    }

    /// The walked cavity reproduces the reference scan exactly (same
    /// triangles, same order) on the inputs that stress it: random
    /// scatter, an exactly cocircular lattice (every cavity a tie),
    /// duplicated points, collinear chart-edge rows, and a thin
    /// anisotropic chart like a metric-scaled coil wall.
    #[test]
    fn walked_cavities_match_the_reference_scan() {
        let mut cases: Vec<(&str, Vec<(f64, f64)>)> = vec![
            ("scatter", scatter(600, 0x9e37_79b9_7f4a_7c15)),
            (
                "lattice",
                (0..20)
                    .flat_map(|i| (0..20).map(move |j| (i as f64, j as f64)))
                    .collect(),
            ),
        ];
        let mut dup = scatter(200, 7);
        dup.extend(scatter(50, 7));
        cases.push(("duplicates", dup));
        let mut rows: Vec<(f64, f64)> = (0..200).map(|i| (i as f64 * 0.01, 0.0)).collect();
        rows.extend((0..200).map(|i| (i as f64 * 0.01, 1.0)));
        rows.extend(
            scatter(100, 3)
                .into_iter()
                .map(|(x, y)| (2.0 * x, 0.05 + 0.9 * y)),
        );
        cases.push(("collinear rows", rows));
        cases.push((
            "thin chart",
            scatter(800, 11)
                .into_iter()
                .map(|(x, y)| (776.0 * x, 0.43 * y))
                .collect(),
        ));
        for (name, pts) in cases {
            let walked = triangulate(&pts, true);
            let scanned = triangulate(&pts, false);
            assert!(scanned.is_some(), "{name}");
            assert_eq!(walked, scanned, "{name}");
        }
    }
}
