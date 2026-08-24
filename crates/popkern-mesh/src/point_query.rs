//! Point↔mesh queries: distance, containment, signed distance.
//!
//! This is the kernel-independent substrate for the offset SDF oracle
//! (design doc §9, PK-323): given a tessellation of a body, classify
//! arbitrary points against it without consulting any BRep machinery —
//! the queries see only triangles, so they cannot inherit a defect
//! from the serving pipeline they are used to judge.
//!
//! Acceleration is a flat AABB tree (median split on centroids), giving
//! logarithmic distance and ray queries; oracle sweeps run hundreds of
//! thousands of queries per fixture.

use popkern_core::geometry::{Point3, Vec3};

use crate::mesh::TriMesh;

/// One triangle stored world-space, with its AABB for pruning.
#[derive(Debug, Clone)]
struct IndexedTriangle {
    a: Point3,
    b: Point3,
    c: Point3,
    lo: Point3,
    hi: Point3,
}

#[derive(Debug, Clone)]
struct Node {
    lo: Point3,
    hi: Point3,
    /// Leaf: `(start, len)` into the triangle permutation.
    /// Internal: `(left_child, right_child)` into `nodes`.
    a: usize,
    b: usize,
    leaf: bool,
}

/// Spatial index over a triangle mesh answering point distance,
/// containment (ray-crossing parity), and signed distance.
///
/// Containment and the sign of [`PointMeshIndex::signed_distance`] are
/// meaningful only for meshes that bound a closed region; distance is
/// meaningful for any mesh. An empty mesh reports infinite distance
/// and contains nothing.
#[derive(Debug, Clone)]
pub struct PointMeshIndex {
    triangles: Vec<IndexedTriangle>,
    /// Permutation of triangle indices grouped by leaf.
    order: Vec<usize>,
    nodes: Vec<Node>,
    root_lo: Point3,
    root_hi: Point3,
}

const LEAF_SIZE: usize = 8;

/// Deterministic skew ray directions for parity counting. A hit near a
/// triangle edge is ambiguous (the neighbouring triangle may double- or
/// zero-count the crossing); on ambiguity the query retries with the
/// next direction.
const PARITY_DIRS: [Vec3; 3] = [
    Vec3 {
        x: 0.577_350_1,
        y: 0.577_351_2,
        z: 0.577_352_3,
    },
    Vec3 {
        x: -0.267_351_9,
        y: 0.802_452_8,
        z: 0.533_651_7,
    },
    Vec3 {
        x: 0.148_253_1,
        y: -0.622_354_9,
        z: 0.768_557_3,
    },
];

impl PointMeshIndex {
    /// Build the index from a combined tessellation.
    pub fn build(mesh: &TriMesh) -> Self {
        let mut triangles = Vec::with_capacity(mesh.triangles.len());
        for t in &mesh.triangles {
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            let lo = Point3::new(
                a.x.min(b.x).min(c.x),
                a.y.min(b.y).min(c.y),
                a.z.min(b.z).min(c.z),
            );
            let hi = Point3::new(
                a.x.max(b.x).max(c.x),
                a.y.max(b.y).max(c.y),
                a.z.max(b.z).max(c.z),
            );
            triangles.push(IndexedTriangle { a, b, c, lo, hi });
        }
        let mut order: Vec<usize> = (0..triangles.len()).collect();
        let mut nodes = Vec::new();
        if !triangles.is_empty() {
            build_node(&triangles, &mut order, 0, triangles.len(), &mut nodes);
        }
        let (root_lo, root_hi) = nodes.first().map_or_else(
            || (Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
            |r| (r.lo, r.hi),
        );
        Self {
            triangles,
            order,
            nodes,
            root_lo,
            root_hi,
        }
    }

    /// Number of triangles indexed.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// Axis-aligned bounds of the whole mesh (`(min, max)`); zero at
    /// the origin for an empty mesh.
    pub fn bounds(&self) -> (Point3, Point3) {
        (self.root_lo, self.root_hi)
    }

    /// Unsigned distance from `p` to the nearest triangle
    /// (branch-and-bound over the AABB tree). Infinite for an empty
    /// mesh.
    pub fn distance(&self, p: Point3) -> f64 {
        if self.nodes.is_empty() {
            return f64::INFINITY;
        }
        let mut best_sq = f64::INFINITY;
        // Depth-first with near-child-first ordering; prune on the
        // squared AABB distance against the current best.
        let mut stack = vec![0_usize];
        while let Some(ni) = stack.pop() {
            let node = &self.nodes[ni];
            if aabb_distance_sq(&node.lo, &node.hi, &p) >= best_sq {
                continue;
            }
            if node.leaf {
                for &ti in &self.order[node.a..node.a + node.b] {
                    let t = &self.triangles[ti];
                    if aabb_distance_sq(&t.lo, &t.hi, &p) >= best_sq {
                        continue;
                    }
                    let d = point_triangle_distance_sq(&p, &t.a, &t.b, &t.c);
                    if d < best_sq {
                        best_sq = d;
                    }
                }
            } else {
                let (l, r) = (node.a, node.b);
                let dl = aabb_distance_sq(&self.nodes[l].lo, &self.nodes[l].hi, &p);
                let dr = aabb_distance_sq(&self.nodes[r].lo, &self.nodes[r].hi, &p);
                // Push the farther child first so the nearer is
                // processed next (tightening the bound early).
                if dl <= dr {
                    stack.push(r);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(r);
                }
            }
        }
        best_sq.sqrt()
    }

    /// Whether `p` lies inside the closed region bounded by the mesh
    /// (ray-crossing parity, orientation-independent). Ambiguous
    /// edge-grazing rays retry along alternate skew directions; if all
    /// retries stay ambiguous the last parity is used.
    pub fn contains(&self, p: Point3) -> bool {
        let mut last = false;
        for dir in &PARITY_DIRS {
            let (crossings, ambiguous) = self.ray_crossings(&p, dir);
            last = crossings % 2 == 1;
            if !ambiguous {
                return last;
            }
        }
        last
    }

    /// Signed distance: negative inside the closed region, positive
    /// outside.
    pub fn signed_distance(&self, p: Point3) -> f64 {
        let d = self.distance(p);
        if self.contains(p) { -d } else { d }
    }

    /// Count ray/triangle crossings from `p` along `dir`; the flag
    /// reports whether any hit landed near a triangle boundary
    /// (parity untrustworthy).
    fn ray_crossings(&self, p: &Point3, dir: &Vec3) -> (usize, bool) {
        if self.nodes.is_empty() {
            return (0, false);
        }
        let inv = Vec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
        let mut crossings = 0_usize;
        let mut ambiguous = false;
        let mut stack = vec![0_usize];
        while let Some(ni) = stack.pop() {
            let node = &self.nodes[ni];
            if !ray_hits_aabb(&node.lo, &node.hi, p, &inv) {
                continue;
            }
            if node.leaf {
                for &ti in &self.order[node.a..node.a + node.b] {
                    let t = &self.triangles[ti];
                    match ray_triangle_hit(p, dir, &t.a, &t.b, &t.c) {
                        RayHit::Interior => crossings += 1,
                        RayHit::NearBoundary => {
                            crossings += 1;
                            ambiguous = true;
                        }
                        RayHit::Miss => {}
                    }
                }
            } else {
                stack.push(node.a);
                stack.push(node.b);
            }
        }
        (crossings, ambiguous)
    }
}

fn build_node(
    triangles: &[IndexedTriangle],
    order: &mut [usize],
    start: usize,
    len: usize,
    nodes: &mut Vec<Node>,
) -> usize {
    let slice = &order[start..start + len];
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &ti in slice {
        let t = &triangles[ti];
        lo = Point3::new(lo.x.min(t.lo.x), lo.y.min(t.lo.y), lo.z.min(t.lo.z));
        hi = Point3::new(hi.x.max(t.hi.x), hi.y.max(t.hi.y), hi.z.max(t.hi.z));
    }
    let this = nodes.len();
    nodes.push(Node {
        lo,
        hi,
        a: start,
        b: len,
        leaf: true,
    });
    if len <= LEAF_SIZE {
        return this;
    }
    // Median split on the longest AABB axis by triangle centroid.
    let ext = hi - lo;
    let axis = if ext.x >= ext.y && ext.x >= ext.z {
        0
    } else if ext.y >= ext.z {
        1
    } else {
        2
    };
    let centroid = |ti: usize| {
        let t = &triangles[ti];
        match axis {
            0 => (t.a.x + t.b.x + t.c.x) / 3.0,
            1 => (t.a.y + t.b.y + t.c.y) / 3.0,
            _ => (t.a.z + t.b.z + t.c.z) / 3.0,
        }
    };
    let mid = len / 2;
    order[start..start + len].select_nth_unstable_by(mid, |&x, &y| {
        centroid(x)
            .partial_cmp(&centroid(y))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let left = build_node(triangles, order, start, mid, nodes);
    let right = build_node(triangles, order, start + mid, len - mid, nodes);
    nodes[this].a = left;
    nodes[this].b = right;
    nodes[this].leaf = false;
    this
}

fn aabb_distance_sq(lo: &Point3, hi: &Point3, p: &Point3) -> f64 {
    let dx = (lo.x - p.x).max(0.0).max(p.x - hi.x);
    let dy = (lo.y - p.y).max(0.0).max(p.y - hi.y);
    let dz = (lo.z - p.z).max(0.0).max(p.z - hi.z);
    dx * dx + dy * dy + dz * dz
}

fn ray_hits_aabb(lo: &Point3, hi: &Point3, origin: &Point3, inv_dir: &Vec3) -> bool {
    let mut tmin = 0.0_f64;
    let mut tmax = f64::INFINITY;
    for (o, i, l, h) in [
        (origin.x, inv_dir.x, lo.x, hi.x),
        (origin.y, inv_dir.y, lo.y, hi.y),
        (origin.z, inv_dir.z, lo.z, hi.z),
    ] {
        let t1 = (l - o) * i;
        let t2 = (h - o) * i;
        tmin = tmin.max(t1.min(t2));
        tmax = tmax.min(t1.max(t2));
        if tmin > tmax {
            return false;
        }
    }
    true
}

enum RayHit {
    Miss,
    Interior,
    /// Hit with a barycentric coordinate near an edge of the triangle
    /// or a ray parameter near zero — parity may double-count.
    NearBoundary,
}

/// Barycentric margin under which a crossing is flagged ambiguous.
const EDGE_MARGIN: f64 = 1.0e-9;

fn ray_triangle_hit(origin: &Point3, dir: &Vec3, a: &Point3, b: &Point3, c: &Point3) -> RayHit {
    // Möller–Trumbore for `origin + t·dir`, `t > 0`.
    let e1 = *b - *a;
    let e2 = *c - *a;
    let p = dir.cross(&e2);
    let det = e1.dot(&p);
    if det.abs() < 1e-15 {
        return RayHit::Miss;
    }
    let inv_det = 1.0 / det;
    let s = *origin - *a;
    let u = s.dot(&p) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return RayHit::Miss;
    }
    let q = s.cross(&e1);
    let v = dir.dot(&q) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return RayHit::Miss;
    }
    let t = e2.dot(&q) * inv_det;
    if t <= 1e-12 {
        return RayHit::Miss;
    }
    let near_edge = u < EDGE_MARGIN
        || v < EDGE_MARGIN
        || u + v > 1.0 - EDGE_MARGIN
        || u > 1.0 - EDGE_MARGIN
        || t < 1e-9;
    if near_edge {
        RayHit::NearBoundary
    } else {
        RayHit::Interior
    }
}

/// Squared distance from `p` to triangle `abc` (Ericson, *Real-Time
/// Collision Detection* §5.1.5).
fn point_triangle_distance_sq(p: &Point3, a: &Point3, b: &Point3, c: &Point3) -> f64 {
    let ab = *b - *a;
    let ac = *c - *a;
    let ap = *p - *a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.dot(&ap);
    }
    let bp = *p - *b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return bp.dot(&bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let denom = d1 - d3;
        let v = if denom.abs() > 0.0 { d1 / denom } else { 0.0 };
        let q = *a + ab * v;
        let dq = *p - q;
        return dq.dot(&dq);
    }
    let cp = *p - *c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return cp.dot(&cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let denom = d2 - d6;
        let w = if denom.abs() > 0.0 { d2 / denom } else { 0.0 };
        let q = *a + ac * w;
        let dq = *p - q;
        return dq.dot(&dq);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let denom = (d4 - d3) + (d5 - d6);
        let w = if denom.abs() > 0.0 {
            (d4 - d3) / denom
        } else {
            0.0
        };
        let q = *b + (*c - *b) * w;
        let dq = *p - q;
        return dq.dot(&dq);
    }
    let denom = va + vb + vc;
    let (v, w) = if denom.abs() > 0.0 {
        (vb / denom, vc / denom)
    } else {
        (0.0, 0.0)
    };
    let q = *a + ab * v + ac * w;
    let dq = *p - q;
    dq.dot(&dq)
}

#[cfg(test)]
#[path = "point_query_tests.rs"]
mod tests;
