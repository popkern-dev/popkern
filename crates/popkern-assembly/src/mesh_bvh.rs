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

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv_near;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TriMesh;

use crate::placement::Placement;
use crate::tri_geometry::{Tri, ray_hit, triangle_distance, triangles_cross};

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
}

impl LocalMesh {
    /// From a tessellation of `brep`, measuring each triangle's deviation.
    pub fn new(brep: &BRep, mesh: &TriMesh) -> LocalMesh {
        let mut out = LocalMesh {
            tris: Vec::with_capacity(mesh.triangles.len()),
            faces: Vec::with_capacity(mesh.triangles.len()),
            eps: Vec::with_capacity(mesh.triangles.len()),
            uvs: Vec::with_capacity(mesh.triangles.len()),
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
            out.eps.push(deviation(brep, face, &corners, &uvs));
            out.tris.push(corners);
            out.faces.push(face);
            out.uvs.push(uvs);
        }
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

impl WorldMesh {
    /// Place `local` at `world` and build the hierarchy.
    pub fn new(local: &LocalMesh, world: &Placement) -> WorldMesh {
        let placed: Vec<Tri> = local
            .tris
            .iter()
            .map(|t| {
                [
                    world.apply_point(&t[0]),
                    world.apply_point(&t[1]),
                    world.apply_point(&t[2]),
                ]
            })
            .collect();
        let mut order: Vec<usize> = (0..placed.len()).collect();
        let mut nodes = Vec::new();
        build(&placed, &local.eps, &mut order, 0, placed.len(), &mut nodes);
        let bounds = nodes.first().map(|n| n.bounds).unwrap_or_else(Aabb::empty);
        WorldMesh {
            tris: order.iter().map(|&i| placed[i]).collect(),
            eps: order.iter().map(|&i| local.eps[i]).collect(),
            local: order,
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
        if self.nodes.is_empty() || other.nodes.is_empty() {
            return None;
        }
        let mut best: Option<Closest> = None;
        let mut lower = f64::INFINITY;
        let mut stack = vec![(0usize, 0usize)];
        while let Some((i, j)) = stack.pop() {
            let (a, b) = (&self.nodes[i], &other.nodes[j]);
            let gap = a.bounds.distance(&b.bounds);
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
                (Some((l, r)), None) => stack.extend([(l, j), (r, j)]),
                (None, Some((l, r))) => stack.extend([(i, l), (i, r)]),
                (Some((al, ar)), Some((bl, br))) => {
                    if a.bounds.diagonal() >= b.bounds.diagonal() {
                        stack.extend([(al, j), (ar, j)]);
                    } else {
                        stack.extend([(i, bl), (i, br)]);
                    }
                }
            }
        }
        best.map(|mut c| {
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

    /// Whether `p` is inside this closed mesh: the majority of three
    /// ray-parity votes along fixed, axis-skew directions.
    pub fn contains(&self, p: &Point3) -> bool {
        let dirs = [
            Vec3::new(0.5773, 0.5774, 0.5773),
            Vec3::new(-0.2672, 0.5345, 0.8017),
            Vec3::new(0.8729, -0.2182, 0.4364),
        ];
        let inside = dirs
            .iter()
            .filter(|d| {
                self.tris
                    .iter()
                    .filter(|t| ray_hit(p, d, t).is_some())
                    .count()
                    % 2
                    == 1
            })
            .count();
        inside >= 2
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
