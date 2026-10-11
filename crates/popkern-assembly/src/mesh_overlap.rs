//! Certain overlap from meshes alone, and BVH point queries (PK-723).
//!
//! Crossing meshes do not prove that the bodies overlap: each mesh
//! strays from its surface by up to its deviation, so two bodies a
//! hair apart can still have crossing meshes. The interference check
//! therefore certified a crossing pair on ever finer meshes and built
//! the exact overlap. A swing that drives the upper arm into the base
//! paid that for every bisection step — meshes of a 300 mm casting at a
//! 0.006 mm chord, then a boolean — though the parts were millimetres
//! deep in each other: 400 s for one swing.
//!
//! A mesh VERTEX lies on its body's true surface (a tessellation's
//! vertices are surface points). If one lies inside the other mesh
//! deeper than that mesh's largest deviation, it lies inside the other
//! true body, and a boundary point of one body strictly inside the
//! other means their interiors overlap. That proof needs no finer mesh
//! and no boolean. When no vertex near a crossing is that deep, the
//! check falls back to certification as before.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use popkern_core::geometry::{Point3, Vec3};

use super::{Aabb, WorldMesh};
use crate::tri_geometry::{ray_hit, triangles_cross};

/// Crossing triangle pairs examined for a deep vertex.
const CROSSINGS: usize = 32;

/// Ray directions for the parity vote (fixed and skew to the axes).
const RAYS: [[f64; 3]; 3] = [
    [0.5773, 0.5774, 0.5773],
    [-0.2672, 0.5345, 0.8017],
    [0.8729, -0.2182, 0.4364],
];

/// Whether the ray from `p` along `d` meets `b` (slab test).
fn ray_meets(p: &Point3, d: &Vec3, b: &Aabb) -> bool {
    let (o, v) = ([p.x, p.y, p.z], [d.x, d.y, d.z]);
    let (mut t0, mut t1) = (0.0_f64, f64::INFINITY);
    for i in 0..3 {
        if v[i].abs() < 1e-300 {
            if o[i] < b.lo[i] || o[i] > b.hi[i] {
                return false;
            }
            continue;
        }
        let (a, c) = ((b.lo[i] - o[i]) / v[i], (b.hi[i] - o[i]) / v[i]);
        t0 = t0.max(a.min(c));
        t1 = t1.min(a.max(c));
        if t0 > t1 {
            return false;
        }
    }
    true
}

/// Distance from `p` to box `b` (0 inside).
fn point_box(p: &Point3, b: &Aabb) -> f64 {
    let c = [p.x, p.y, p.z];
    (0..3)
        .map(|i| (b.lo[i] - c[i]).max(0.0).max(c[i] - b.hi[i]))
        .map(|d| d * d)
        .sum::<f64>()
        .sqrt()
}

/// Distance from `p` to triangle `t`.
fn point_triangle(p: &Point3, t: &crate::mesh_bvh::Tri) -> f64 {
    crate::tri_geometry::triangle_distance(&[*p, *p, *p], t).0
}

impl WorldMesh {
    /// Triangles the ray from `p` along `d` hits (BVH-pruned).
    fn ray_hits(&self, p: &Point3, d: &Vec3) -> usize {
        if self.nodes.is_empty() {
            return 0;
        }
        let mut hits = 0;
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let node = &self.nodes[n];
            if !ray_meets(p, d, &node.bounds) {
                continue;
            }
            match node.children {
                Some((l, r)) => stack.extend([l, r]),
                None => {
                    hits += (node.first..node.first + node.count)
                        .filter(|&t| ray_hit(p, d, &self.tris[t]).is_some())
                        .count();
                }
            }
        }
        hits
    }

    /// Whether `p` is inside this closed mesh: the majority of three
    /// ray-parity votes along fixed, axis-skew directions.
    pub fn contains(&self, p: &Point3) -> bool {
        RAYS.iter()
            .filter(|d| self.ray_hits(p, &Vec3::new(d[0], d[1], d[2])) % 2 == 1)
            .count()
            >= 2
    }

    /// Distance from `p` to this mesh.
    pub fn point_distance(&self, p: &Point3) -> f64 {
        let mut best = f64::INFINITY;
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let Some(node) = self.nodes.get(n) else {
                continue;
            };
            if point_box(p, &node.bounds) >= best {
                continue;
            }
            match node.children {
                Some((l, r)) => {
                    let (dl, dr) = (
                        point_box(p, &self.nodes[l].bounds),
                        point_box(p, &self.nodes[r].bounds),
                    );
                    if dl <= dr {
                        stack.extend([r, l]);
                    } else {
                        stack.extend([l, r]);
                    }
                }
                None => {
                    for t in node.first..node.first + node.count {
                        best = best.min(point_triangle(p, &self.tris[t]));
                    }
                }
            }
        }
        best
    }

    /// Up to `limit` crossing triangle pairs (indices here, there).
    fn crossings(&self, other: &WorldMesh, limit: usize) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return out;
        }
        let mut stack = vec![(0usize, 0usize)];
        while let Some((i, j)) = stack.pop() {
            let (a, b) = (&self.nodes[i], &other.nodes[j]);
            if a.bounds.distance(&b.bounds) > 0.0 {
                continue;
            }
            match (a.children, b.children) {
                (None, None) => {
                    for ta in a.first..a.first + a.count {
                        for tb in b.first..b.first + b.count {
                            if triangles_cross(&self.tris[ta], &other.tris[tb]) {
                                out.push((ta, tb));
                                if out.len() >= limit {
                                    return out;
                                }
                            }
                        }
                    }
                }
                (Some((l, r)), None) => stack.extend([(l, j), (r, j)]),
                (None, Some((l, r))) => stack.extend([(i, l), (i, r)]),
                (Some((al, ar)), Some(_)) => stack.extend([(al, j), (ar, j)]),
            }
        }
        out
    }
}

/// Whether the bodies meshed by `a` and `b` certainly overlap: a vertex
/// of a crossing triangle lies inside the other mesh deeper than that
/// mesh's largest deviation plus `margin` (see the module docs).
pub fn certain_overlap(a: &WorldMesh, b: &WorldMesh, margin: f64) -> bool {
    let deep_in =
        |p: &Point3, m: &WorldMesh| m.point_distance(p) > m.max_eps() + margin && m.contains(p);
    a.crossings(b, CROSSINGS).iter().any(|&(ta, tb)| {
        a.tris[ta].iter().any(|p| deep_in(p, b)) || b.tris[tb].iter().any(|p| deep_in(p, a))
    })
}
