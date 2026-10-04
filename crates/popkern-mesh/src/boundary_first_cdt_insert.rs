//! Incremental Steiner insertion into a constrained Delaunay
//! triangulation (Bowyer–Watson with constraints).
//!
//! Deviation refinement adds points round after round; re-triangulating
//! the whole face per round cost a full Bowyer–Watson each time, and on
//! long thin helical charts that is far from linear (a 5.8k-point thread
//! flank took 0.95 s, eight rounds on three flanks per tessellation:
//! 1.4 s -> 26 s for a 24-turn rod; PK-651). Here each point is located
//! by a visibility walk, its cavity (the triangles whose circumcircle
//! holds it, grown without crossing a constrained edge) is removed, and
//! the cavity is re-starred from the point. Constrained edges are never
//! removed, so ring segments stay edges by construction; the result is
//! the constrained Delaunay triangulation of the grown point set.

use std::collections::HashSet;

use super::{Mesh, in_circle, key, orient};

/// A constrained triangulation that accepts new interior points.
pub(in crate::boundary_first) struct ConstrainedMesh {
    mesh: Mesh,
    constrained: HashSet<(usize, usize)>,
    /// Walk start for the next location (the last re-starred triangle).
    hint: usize,
}

impl ConstrainedMesh {
    /// Wrap the triangulation `tris` of `pts` (any orientation) whose
    /// `segments` are edges. `None` when a triangle is degenerate.
    pub(in crate::boundary_first) fn new(
        pts: &[(f64, f64)],
        tris: &[[usize; 3]],
        segments: &[(usize, usize)],
    ) -> Option<Self> {
        Some(ConstrainedMesh {
            mesh: Mesh::new(pts, tris)?,
            constrained: segments.iter().map(|&(a, b)| key(a, b)).collect(),
            hint: 0,
        })
    }

    /// The current triangles (counter-clockwise).
    pub(in crate::boundary_first) fn triangles(&self) -> Vec<[usize; 3]> {
        self.mesh.tris.clone()
    }

    /// Insert `p`; returns its vertex index (the next index in the point
    /// list), or `None`, leaving the mesh unchanged, when `p` lies outside
    /// every triangle, on a constrained edge or on an existing vertex,
    /// or its cavity is not star-shaped from it.
    pub(in crate::boundary_first) fn insert(&mut self, p: (f64, f64)) -> Option<usize> {
        let start = self.locate(p)?;
        let t = self.mesh.tris[start];
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if self.mesh.pts[a] == p {
                return None;
            }
            if orient(self.mesh.pts[a], self.mesh.pts[b], p) == 0.0
                && self.constrained.contains(&key(a, b))
            {
                return None;
            }
        }
        // Cavity: grown across unconstrained edges into triangles whose
        // circumcircle holds p (a point on the start triangle's edge
        // lies inside the neighbour's circumcircle, so it joins too).
        let mut cavity = vec![start];
        let mut seen: HashSet<usize> = HashSet::from([start]);
        let mut next = 0;
        while next < cavity.len() {
            let ti = cavity[next];
            next += 1;
            let t = self.mesh.tris[ti];
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if self.constrained.contains(&key(a, b)) {
                    continue;
                }
                let Some(n) = self.across(a, b, ti) else {
                    continue;
                };
                if seen.insert(n) {
                    let m = self.mesh.tris[n];
                    let (pa, pb, pc) = (
                        self.mesh.pts[m[0]],
                        self.mesh.pts[m[1]],
                        self.mesh.pts[m[2]],
                    );
                    if in_circle(pa, pb, pc, p) {
                        cavity.push(n);
                    }
                }
            }
        }
        // Boundary: directed cavity edges whose other side is outside.
        let members: HashSet<usize> = cavity.iter().copied().collect();
        let mut rim: Vec<(usize, usize)> = Vec::new();
        for &ti in &cavity {
            let t = self.mesh.tris[ti];
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let inside = self.across(a, b, ti).is_some_and(|n| members.contains(&n));
                if !inside {
                    // Star-shaped from p: every new triangle (a, b, p)
                    // is counter-clockwise.
                    if orient(self.mesh.pts[a], self.mesh.pts[b], p) <= 0.0 {
                        return None;
                    }
                    rim.push((a, b));
                }
            }
        }
        let v = self.mesh.pts.len();
        self.mesh.pts.push(p);
        self.mesh.vertex_tris.push(Vec::new());
        for &ti in &cavity {
            self.unlink(ti);
        }
        let mut slots = cavity.into_iter();
        for (a, b) in rim {
            let t = [a, b, v];
            let ti = match slots.next() {
                Some(ti) => {
                    self.mesh.tris[ti] = t;
                    ti
                }
                None => {
                    self.mesh.tris.push(t);
                    self.mesh.tris.len() - 1
                }
            };
            self.link(ti);
            self.hint = ti;
        }
        // A cavity of k triangles has k + 2 rim edges, so every slot was
        // reused; nothing is left to remove.
        Some(v)
    }

    /// The triangle holding `p` (boundary inclusive): a visibility walk
    /// from the hint, then a scan when the walk leaves the mesh (a
    /// concave trim region).
    fn locate(&self, p: (f64, f64)) -> Option<usize> {
        let n = self.mesh.tris.len();
        if n == 0 {
            return None;
        }
        let holds = |ti: usize| {
            let t = self.mesh.tris[ti];
            (0..3).all(|k| orient(self.mesh.pts[t[k]], self.mesh.pts[t[(k + 1) % 3]], p) >= 0.0)
        };
        let mut ti = self.hint.min(n - 1);
        for step in 0..n {
            let t = self.mesh.tris[ti];
            let mut moved = false;
            for r in 0..3 {
                let k = (r + step) % 3;
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if orient(self.mesh.pts[a], self.mesh.pts[b], p) < 0.0 {
                    match self.across(a, b, ti) {
                        Some(next) => {
                            ti = next;
                            moved = true;
                        }
                        None => return (0..n).find(|&i| holds(i)),
                    }
                    break;
                }
            }
            if !moved {
                return Some(ti);
            }
        }
        (0..n).find(|&i| holds(i))
    }

    /// The triangle across edge `(a, b)` from `ti`.
    fn across(&self, a: usize, b: usize, ti: usize) -> Option<usize> {
        self.mesh
            .edges
            .get(&key(a, b))?
            .iter()
            .copied()
            .find(|&x| x != ti)
    }

    /// Drop triangle `ti` from the incidence maps (its slot is reused).
    fn unlink(&mut self, ti: usize) {
        let t = self.mesh.tris[ti];
        for k in 0..3 {
            let e = key(t[k], t[(k + 1) % 3]);
            if let Some(inc) = self.mesh.edges.get_mut(&e) {
                inc.retain(|&x| x != ti);
                if inc.is_empty() {
                    self.mesh.edges.remove(&e);
                }
            }
            self.mesh.vertex_tris[t[k]].retain(|&x| x != ti);
        }
    }

    /// Enter triangle `ti` into the incidence maps.
    fn link(&mut self, ti: usize) {
        let t = self.mesh.tris[ti];
        for k in 0..3 {
            self.mesh
                .edges
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(ti);
            self.mesh.vertex_tris[t[k]].push(ti);
        }
    }
}
