//! Constrained-Delaunay segment recovery for boundary-first charts.
//!
//! The unconstrained Bowyer–Watson triangulation keeps a boundary
//! segment only when its diametral disk is free of other points. Where
//! two boundary edges meet at a sharp input angle, each edge's ring
//! vertices sit inside the other's disks, and midpoint splitting
//! ping-pongs instead of converging (Ruppert's small-angle case;
//! measured: an M4 socket-cap tip cone ran out of encroachment rounds
//! on a 0.03-wide wedge and fell back with 131 open mesh edges). A
//! segment's crossing edges can instead be flipped away (Sloan's
//! recovery), after which Lawson flips restore the Delaunay property
//! on every edge that is not a segment.
//!
//! Recovery is exact-predicate and all-or-nothing: it reports failure
//! (and the caller keeps the unconstrained triangulation, whose
//! coverage check then requests the split rounds as before) when a
//! point lies on a segment's interior — a T-point only a conforming
//! split resolves — or when segments cross or the flip budget runs out.

use std::collections::{HashMap, HashSet};

/// Orientation of `c` relative to the directed line `a → b` (exact).
fn orient(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.0, y: a.1 },
        robust::Coord { x: b.0, y: b.1 },
        robust::Coord { x: c.0, y: c.1 },
    )
}

/// Whether `d` lies strictly inside the circumcircle of the
/// counter-clockwise triangle `(a, b, c)` (exact).
fn in_circle(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64)) -> bool {
    robust::incircle(
        robust::Coord { x: a.0, y: a.1 },
        robust::Coord { x: b.0, y: b.1 },
        robust::Coord { x: c.0, y: c.1 },
        robust::Coord { x: d.0, y: d.1 },
    ) > 0.0
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Counter-clockwise triangles with edge → triangles and vertex →
/// triangles incidence, both kept current under flips.
struct Mesh {
    pts: Vec<(f64, f64)>,
    tris: Vec<[usize; 3]>,
    edges: HashMap<(usize, usize), Vec<usize>>,
    vertex_tris: Vec<Vec<usize>>,
}

impl Mesh {
    fn new(pts: &[(f64, f64)], tris: &[[usize; 3]]) -> Option<Self> {
        let mut mesh = Mesh {
            pts: pts.to_vec(),
            tris: Vec::with_capacity(tris.len()),
            edges: HashMap::new(),
            vertex_tris: vec![Vec::new(); pts.len()],
        };
        for &t in tris {
            let o = orient(pts[t[0]], pts[t[1]], pts[t[2]]);
            if o == 0.0 {
                return None;
            }
            let t = if o > 0.0 { t } else { [t[0], t[2], t[1]] };
            let ti = mesh.tris.len();
            mesh.tris.push(t);
            for &v in &t {
                mesh.vertex_tris.get_mut(v)?.push(ti);
            }
            for k in 0..3 {
                mesh.edges
                    .entry(key(t[k], t[(k + 1) % 3]))
                    .or_default()
                    .push(ti);
            }
        }
        Some(mesh)
    }

    /// The two triangles of interior edge `(p, q)` as `(p, q, r)` and
    /// `(q, p, s)` (both counter-clockwise) with their indices.
    fn quad(&self, p: usize, q: usize) -> Option<(usize, usize, usize, usize)> {
        let inc = self.edges.get(&key(p, q))?;
        if inc.len() != 2 {
            return None;
        }
        let apex = |ti: usize, from: usize, to: usize| -> Option<usize> {
            let t = self.tris[ti];
            (0..3)
                .find(|&k| t[k] == from && t[(k + 1) % 3] == to)
                .map(|k| t[(k + 2) % 3])
        };
        let (t1, t2) = if apex(inc[0], p, q).is_some() {
            (inc[0], inc[1])
        } else {
            (inc[1], inc[0])
        };
        Some((apex(t1, p, q)?, apex(t2, q, p)?, t1, t2))
    }

    /// Whether the quad around `(p, q)` is strictly convex, so the
    /// diagonal may flip to `(r, s)`.
    fn flippable(&self, p: usize, q: usize, r: usize, s: usize) -> bool {
        let (pp, pq, pr, ps) = (self.pts[p], self.pts[q], self.pts[r], self.pts[s]);
        orient(pr, ps, pp) * orient(pr, ps, pq) < 0.0
            && orient(pp, pq, pr) * orient(pp, pq, ps) < 0.0
    }

    /// Replace diagonal `(p, q)` by `(r, s)`.
    fn flip(&mut self, p: usize, q: usize, r: usize, s: usize, t_pqr: usize, t_qps: usize) {
        for (ti, t) in [(t_pqr, [p, q, r]), (t_qps, [q, p, s])] {
            for k in 0..3 {
                if let Some(inc) = self.edges.get_mut(&key(t[k], t[(k + 1) % 3])) {
                    inc.retain(|&x| x != ti);
                }
            }
        }
        self.edges.remove(&key(p, q));
        // Quad order (counter-clockwise): p, s, q, r. (p, q, r) becomes
        // (p, s, r) and (q, p, s) becomes (s, q, r).
        self.vertex_tris[q].retain(|&x| x != t_pqr);
        self.vertex_tris[s].push(t_pqr);
        self.vertex_tris[p].retain(|&x| x != t_qps);
        self.vertex_tris[r].push(t_qps);
        let new = [(t_pqr, [p, s, r]), (t_qps, [s, q, r])];
        for (ti, t) in new {
            self.tris[ti] = t;
            for k in 0..3 {
                self.edges
                    .entry(key(t[k], t[(k + 1) % 3]))
                    .or_default()
                    .push(ti);
            }
        }
    }

    /// The edges the open segment `(a, b)` crosses, in order, found by
    /// walking the triangles it passes through from `a` (local: the
    /// cost is the crossing count, not the mesh size). `None` when a
    /// vertex lies on the segment's interior (a T-point) or the walk
    /// cannot proceed.
    fn crossing_edges(&self, a: usize, b: usize) -> Option<Vec<(usize, usize)>> {
        let (pa, pb) = (self.pts[a], self.pts[b]);
        let on_open_segment = |v: usize| {
            let p = self.pts[v];
            orient(pa, pb, p) == 0.0 && {
                let t = (p.0 - pa.0) * (pb.0 - pa.0) + (p.1 - pa.1) * (pb.1 - pa.1);
                let l2 = (pb.0 - pa.0).powi(2) + (pb.1 - pa.1).powi(2);
                t > 0.0 && t < l2
            }
        };
        // The triangle at `a` whose opposite edge (p, q) the segment
        // leaves through: p right of a -> b, q left.
        let mut edge = None;
        for &ti in &self.vertex_tris[a] {
            let t = self.tris[ti];
            let k = (0..3).find(|&k| t[k] == a)?;
            let (p, q) = (t[(k + 1) % 3], t[(k + 2) % 3]);
            if on_open_segment(p) || on_open_segment(q) {
                return None;
            }
            if orient(pa, pb, self.pts[p]) < 0.0 && orient(pa, pb, self.pts[q]) > 0.0 {
                edge = Some((p, q));
                break;
            }
        }
        let (mut p, mut q) = edge?;
        let mut crossed = vec![(p, q)];
        for _ in 0..self.tris.len() {
            // The triangle beyond (p, q): its apex s.
            let (_, s, _, _) = self.quad(p, q)?;
            if s == b {
                return Some(crossed);
            }
            if on_open_segment(s) {
                return None;
            }
            if orient(pa, pb, self.pts[s]) < 0.0 {
                p = s;
            } else {
                q = s;
            }
            crossed.push((p, q));
        }
        None
    }

    /// Whether edge `(p, q)` crosses the open segment `(a, b)` at an
    /// interior point of both.
    fn crosses(&self, p: usize, q: usize, a: usize, b: usize) -> bool {
        if p == a || p == b || q == a || q == b {
            return false;
        }
        let (pa, pb, pp, pq) = (self.pts[a], self.pts[b], self.pts[p], self.pts[q]);
        orient(pa, pb, pp) * orient(pa, pb, pq) < 0.0
            && orient(pp, pq, pa) * orient(pp, pq, pb) < 0.0
    }
}

/// Make every `segments` pair an edge of the triangulation `tris` of
/// `pts` by edge flips, then restore the Delaunay property on every
/// non-segment edge. Returns the constrained triangulation
/// (counter-clockwise), or `None` when some segment cannot be
/// recovered (see the module docs); `tris` is never modified.
pub(super) fn recover_segments(
    pts: &[(f64, f64)],
    tris: &[[usize; 3]],
    segments: &[(usize, usize)],
) -> Option<Vec<[usize; 3]>> {
    let mut mesh = Mesh::new(pts, tris)?;
    let constrained: HashSet<(usize, usize)> = segments.iter().map(|&(a, b)| key(a, b)).collect();
    let mut touched: Vec<(usize, usize)> = Vec::new();
    for &(a, b) in segments {
        if a == b || mesh.edges.contains_key(&key(a, b)) {
            continue;
        }
        // A point on the segment's interior can never be flipped
        // around: the caller's conforming split owns it.
        let mut queue: std::collections::VecDeque<(usize, usize)> =
            mesh.crossing_edges(a, b)?.into_iter().collect();
        // Sloan's bound: a crossing edge flips at most once per pass,
        // and a pass with no flip on a non-empty queue cannot happen
        // in exact arithmetic; the budget only guards a broken input
        // (crossing segments).
        let mut budget = 64 * queue.len() + 64;
        while let Some((p, q)) = queue.pop_front() {
            if budget == 0 {
                return None;
            }
            budget -= 1;
            if constrained.contains(&key(p, q)) {
                return None; // two segments cross
            }
            let (r, s, t_pqr, t_qps) = mesh.quad(p, q)?;
            if !mesh.flippable(p, q, r, s) {
                queue.push_back((p, q));
                continue;
            }
            mesh.flip(p, q, r, s, t_pqr, t_qps);
            if mesh.crosses(r, s, a, b) {
                queue.push_back((r, s));
            } else {
                touched.push((r, s));
            }
        }
        if !mesh.edges.contains_key(&key(a, b)) {
            return None;
        }
    }
    // Lawson: restore the empty-circle property around every edge the
    // recovery created, never flipping a segment.
    let mut budget = 64 * touched.len() + 64;
    while let Some((p, q)) = touched.pop() {
        if budget == 0 {
            break;
        }
        budget -= 1;
        if constrained.contains(&key(p, q)) {
            continue;
        }
        let Some((r, s, t_pqr, t_qps)) = mesh.quad(p, q) else {
            continue;
        };
        if !mesh.flippable(p, q, r, s) || !in_circle(pts[p], pts[q], pts[r], pts[s]) {
            continue;
        }
        mesh.flip(p, q, r, s, t_pqr, t_qps);
        touched.extend([(p, s), (s, q), (q, r), (r, p)]);
    }
    Some(mesh.tris)
}

#[path = "boundary_first_cdt_insert.rs"]
mod insert;
pub(super) use insert::ConstrainedMesh;

#[cfg(test)]
#[path = "boundary_first_cdt_tests.rs"]
mod tests;
