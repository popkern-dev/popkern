//! A body's triangles in world coordinates with a bounding-volume
//! hierarchy, for the interference check: branch-and-bound distance
//! between two meshes, crossing detection, and point containment by
//! ray parity.
//!
//! Every triangle carries its MEASURED deviation `eps` from the true
//! surface (the largest distance from its centroid and edge midpoints
//! to the face's surface, see [`deviation`]), so a mesh distance `d`
//! between triangles `ta`, `tb` bounds the true distance of the faces
//! they sample from below by `d − eps(ta) − eps(tb)`.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv_near;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TriMesh;

use crate::placement::Placement;
use crate::tri_geometry::{Tri, triangle_distance, triangles_cross};

/// An axis-aligned box.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    /// Lower corner.
    pub lo: [f64; 3],
    /// Upper corner.
    pub hi: [f64; 3],
}

impl Aabb {
    fn empty() -> Aabb {
        Aabb {
            lo: [f64::INFINITY; 3],
            hi: [f64::NEG_INFINITY; 3],
        }
    }

    fn add(&mut self, p: &Point3) {
        for (k, x) in [p.x, p.y, p.z].into_iter().enumerate() {
            self.lo[k] = self.lo[k].min(x);
            self.hi[k] = self.hi[k].max(x);
        }
    }

    /// Smallest distance between points of the two boxes (0 when they
    /// overlap).
    pub fn distance(&self, o: &Aabb) -> f64 {
        let mut s = 0.0;
        for k in 0..3 {
            let gap = (o.lo[k] - self.hi[k]).max(self.lo[k] - o.hi[k]).max(0.0);
            s += gap * gap;
        }
        s.sqrt()
    }

    /// The diagonal length.
    pub fn diagonal(&self) -> f64 {
        (0..3)
            .map(|k| (self.hi[k] - self.lo[k]).powi(2))
            .sum::<f64>()
            .sqrt()
    }
}

/// The triangles of one definition body, before placement: corners,
/// the face each samples, and each triangle's measured deviation.
#[derive(Debug, Clone)]
pub struct LocalMesh {
    /// Triangle corners.
    pub tris: Vec<Tri>,
    /// The BRep face each triangle samples.
    pub faces: Vec<FaceHandle>,
    /// Each triangle's measured deviation from its face's surface.
    pub eps: Vec<f64>,
    /// Per-corner surface parameters (projection hints).
    pub uvs: Vec<[(f64, f64); 3]>,
    /// The box of every corner.
    pub bounds: Aabb,
    /// The largest triangle deviation.
    pub max_eps: f64,
    /// The hierarchy's triangle order and nodes, built once here; a
    /// placement transforms the triangles and refits the node boxes.
    order: Vec<usize>,
    nodes: Vec<Node>,
}

impl LocalMesh {
    /// From a tessellation of `brep`, measuring each triangle's deviation.
    pub fn new(brep: &BRep, mesh: &TriMesh) -> LocalMesh {
        let mut out = LocalMesh {
            tris: Vec::with_capacity(mesh.triangles.len()),
            faces: Vec::with_capacity(mesh.triangles.len()),
            eps: Vec::with_capacity(mesh.triangles.len()),
            uvs: Vec::with_capacity(mesh.triangles.len()),
            bounds: Aabb::empty(),
            max_eps: 0.0,
            order: Vec::new(),
            nodes: Vec::new(),
        };
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let corners = [
                mesh.vertices[tri[0]].position,
                mesh.vertices[tri[1]].position,
                mesh.vertices[tri[2]].position,
            ];
            let uvs = [
                mesh.vertices[tri[0]].uv,
                mesh.vertices[tri[1]].uv,
                mesh.vertices[tri[2]].uv,
            ];
            let face = mesh.face_map[t];
            let eps = deviation(brep, face, &corners, &uvs);
            out.max_eps = out.max_eps.max(eps);
            for p in &corners {
                for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                    out.bounds.lo[i] = out.bounds.lo[i].min(c);
                    out.bounds.hi[i] = out.bounds.hi[i].max(c);
                }
            }
            out.eps.push(eps);
            out.tris.push(corners);
            out.faces.push(face);
            out.uvs.push(uvs);
        }
        out.order = (0..out.tris.len()).collect();
        build(
            &out.tris,
            &out.eps,
            &mut out.order,
            0,
            out.tris.len(),
            &mut out.nodes,
        );
        out
    }
}

/// The largest distance from the triangle's centroid and edge
/// midpoints to the face's surface (projected from the interpolated
/// parameters), or a conservative stand-in (the triangle's longest
/// edge) when the projection fails.
pub fn deviation(brep: &BRep, face: FaceHandle, tri: &Tri, uvs: &[(f64, f64); 3]) -> f64 {
    let longest = (0..3)
        .map(|i| (tri[i] - tri[(i + 1) % 3]).length())
        .fold(0.0, f64::max);
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return longest;
    };
    let samples: [[f64; 3]; 4] = [
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
    ];
    let mut worst: f64 = 0.0;
    for w in samples {
        let p = Point3::new(
            w[0] * tri[0].x + w[1] * tri[1].x + w[2] * tri[2].x,
            w[0] * tri[0].y + w[1] * tri[1].y + w[2] * tri[2].y,
            w[0] * tri[0].z + w[1] * tri[1].z + w[2] * tri[2].z,
        );
        let hint = (
            w[0] * uvs[0].0 + w[1] * uvs[1].0 + w[2] * uvs[2].0,
            w[0] * uvs[0].1 + w[1] * uvs[1].1 + w[2] * uvs[2].1,
        );
        let Some((u, v)) = project_point_to_uv_near(&surface.geometry, &p, hint.0, hint.1) else {
            return longest;
        };
        worst = worst.max((surface.geometry.point_at(u, v) - p).length());
    }
    worst
}

#[derive(Debug, Clone)]
struct Node {
    bounds: Aabb,
    /// Largest triangle deviation below this node.
    eps: f64,
    /// Children, or `None` for a leaf over `first..first + count`.
    children: Option<(usize, usize)>,
    first: usize,
    count: usize,
}

/// One placed body's triangles in world coordinates with a BVH.
#[derive(Debug, Clone)]
pub struct WorldMesh {
    /// World triangles, in BVH order.
    pub tris: Vec<Tri>,
    /// Index into the local mesh's triangle arrays, in BVH order.
    pub local: Vec<usize>,
    eps: Vec<f64>,
    nodes: Vec<Node>,
    /// The whole mesh's bounds.
    pub bounds: Aabb,
}

const LEAF: usize = 4;

#[path = "mesh_overlap.rs"]
mod overlap;
pub use overlap::certain_overlap;

impl WorldMesh {
    /// Place `local` at `world`: its triangles transformed, in the
    /// hierarchy order `local` built once, with every node box refitted
    /// to them bottom-up (children follow their parent in `nodes`).
    ///
    /// A sweep re-places every moving occurrence at every step, and
    /// rebuilding each one's hierarchy (a sort per level) was most of
    /// a base-yaw swing step (PK-760). The refitted boxes are the exact
    /// world bounds of the same triangles, so every query (closest
    /// pair, crossing) visits a valid hierarchy and returns the same
    /// geometric answer.
    pub fn new(local: &LocalMesh, world: &Placement) -> WorldMesh {
        let tris: Vec<Tri> = local
            .order
            .iter()
            .map(|&i| {
                let t = &local.tris[i];
                [
                    world.apply_point(&t[0]),
                    world.apply_point(&t[1]),
                    world.apply_point(&t[2]),
                ]
            })
            .collect();
        let mut nodes = local.nodes.clone();
        for k in (0..nodes.len()).rev() {
            let mut bounds = Aabb::empty();
            match nodes[k].children {
                None => {
                    for t in &tris[nodes[k].first..nodes[k].first + nodes[k].count] {
                        for p in t {
                            bounds.add(p);
                        }
                    }
                }
                Some((l, r)) => {
                    for c in [l, r] {
                        for axis in 0..3 {
                            bounds.lo[axis] = bounds.lo[axis].min(nodes[c].bounds.lo[axis]);
                            bounds.hi[axis] = bounds.hi[axis].max(nodes[c].bounds.hi[axis]);
                        }
                    }
                }
            }
            nodes[k].bounds = bounds;
        }
        let bounds = nodes.first().map(|n| n.bounds).unwrap_or_else(Aabb::empty);
        WorldMesh {
            tris,
            eps: local.order.iter().map(|&i| local.eps[i]).collect(),
            local: local.order.clone(),
            nodes,
            bounds,
        }
    }

    /// The largest triangle deviation of the whole mesh.
    pub fn max_eps(&self) -> f64 {
        self.nodes.first().map(|n| n.eps).unwrap_or(0.0)
    }

    /// Branch-and-bound closest pair against `other`: `(mesh distance,
    /// certified lower bound, point here, point there, triangle here,
    /// triangle there)` (triangle indices in BVH order). The lower
    /// bound subtracts both triangles' deviations at every pair, so it
    /// bounds the true face distance.
    pub fn closest(&self, other: &WorldMesh) -> Option<Closest> {
        self.closest_within(other, f64::INFINITY)
    }

    /// [`WorldMesh::closest`] where it matters: `None` when the meshes
    /// are certainly farther apart than `cutoff` (every node pair's gap
    /// less both deviations exceeds it). Node pairs beyond the cutoff
    /// are never opened, and nearer child pairs are visited first, so a
    /// pair well clear of a swing's clearance costs a few box tests.
    pub fn closest_within(&self, other: &WorldMesh, cutoff: f64) -> Option<Closest> {
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return None;
        }
        let gap_of = |i: usize, j: usize| self.nodes[i].bounds.distance(&other.nodes[j].bounds);
        let mut best: Option<Closest> = None;
        let mut lower = f64::INFINITY;
        let mut stack = vec![(0usize, 0usize)];
        while let Some((i, j)) = stack.pop() {
            let (a, b) = (&self.nodes[i], &other.nodes[j]);
            let gap = a.bounds.distance(&b.bounds);
            if gap - a.eps - b.eps > cutoff {
                continue;
            }
            let best_d = best.as_ref().map_or(f64::INFINITY, |c| c.distance);
            // A node pair can improve neither the mesh distance nor the
            // certified lower bound once both its gaps exceed them.
            if gap >= best_d && gap - a.eps - b.eps >= lower {
                continue;
            }
            match (a.children, b.children) {
                (None, None) => {
                    for ta in a.first..a.first + a.count {
                        for tb in b.first..b.first + b.count {
                            let (d, pa, pb) = triangle_distance(&self.tris[ta], &other.tris[tb]);
                            lower = lower.min((d - self.eps[ta] - other.eps[tb]).max(0.0));
                            if d < best.as_ref().map_or(f64::INFINITY, |c| c.distance) {
                                best = Some(Closest {
                                    distance: d,
                                    lower_bound: 0.0,
                                    here: pa,
                                    there: pb,
                                    tri_here: ta,
                                    tri_there: tb,
                                });
                            }
                        }
                    }
                }
                (Some(_), _) | (None, Some(_)) => {
                    // Split the larger node; the nearer child pair is
                    // pushed last, so it is opened first.
                    let split_a = match (a.children, b.children) {
                        (Some(_), None) => true,
                        (None, Some(_)) => false,
                        _ => a.bounds.diagonal() >= b.bounds.diagonal(),
                    };
                    let (x, y) = match (split_a, a.children, b.children) {
                        (true, Some((l, r)), _) => ((l, j), (r, j)),
                        (false, _, Some((l, r))) => ((i, l), (i, r)),
                        _ => continue,
                    };
                    if gap_of(x.0, x.1) <= gap_of(y.0, y.1) {
                        stack.extend([y, x]);
                    } else {
                        stack.extend([x, y]);
                    }
                }
            }
        }
        // Every pruned node pair lies beyond the cutoff: the certified
        // bound stands as `lower` while it is within it.
        best.filter(|_| lower <= cutoff).map(|mut c| {
            c.lower_bound = lower.min(c.distance);
            c
        })
    }

    /// Whether any triangle here crosses one of `other`.
    pub fn crosses(&self, other: &WorldMesh) -> bool {
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return false;
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
                                return true;
                            }
                        }
                    }
                }
                (Some((l, r)), None) => stack.extend([(l, j), (r, j)]),
                (None, Some((l, r))) => stack.extend([(i, l), (i, r)]),
                (Some((al, ar)), Some(_)) => stack.extend([(al, j), (ar, j)]),
            }
        }
        false
    }
}

/// The closest pair two meshes found.
#[derive(Debug, Clone)]
pub struct Closest {
    /// Distance between the meshes.
    pub distance: f64,
    /// Certified lower bound on the true distance of the bodies.
    pub lower_bound: f64,
    /// The closest point on the first mesh.
    pub here: Point3,
    /// The closest point on the second mesh.
    pub there: Point3,
    /// Its triangle on the first mesh (BVH order).
    pub tri_here: usize,
    /// Its triangle on the second mesh (BVH order).
    pub tri_there: usize,
}

fn build(
    tris: &[Tri],
    eps: &[f64],
    order: &mut [usize],
    first: usize,
    count: usize,
    nodes: &mut Vec<Node>,
) -> usize {
    let mut bounds = Aabb::empty();
    let mut node_eps: f64 = 0.0;
    for &i in &order[first..first + count] {
        for p in &tris[i] {
            bounds.add(p);
        }
        node_eps = node_eps.max(eps[i]);
    }
    let index = nodes.len();
    nodes.push(Node {
        bounds,
        eps: node_eps,
        children: None,
        first,
        count,
    });
    if count <= LEAF {
        return index;
    }
    let axis = (0..3)
        .max_by(|&a, &b| (bounds.hi[a] - bounds.lo[a]).total_cmp(&(bounds.hi[b] - bounds.lo[b])))
        .unwrap_or(0);
    let centre = |i: usize| {
        let t = &tris[i];
        [
            t[0].x + t[1].x + t[2].x,
            t[0].y + t[1].y + t[2].y,
            t[0].z + t[1].z + t[2].z,
        ][axis]
    };
    order[first..first + count].sort_by(|&a, &b| centre(a).total_cmp(&centre(b)));
    let half = count / 2;
    let left = build(tris, eps, order, first, half, nodes);
    let right = build(tris, eps, order, first + half, count - half, nodes);
    nodes[index].children = Some((left, right));
    index
}
