//! Tessellation distances and containment for the generator's own
//! checks (split from `geometry.rs` for the file-size limit).

use std::collections::HashMap;

use popkern_core::geometry::{Point3, Vec3};
use popkern_mesh::TriMesh;

pub fn point_triangle_distance(p: Point3, a: Point3, b: Point3, c: Point3) -> f64 {
    // Ericson, closest point on a triangle.
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(&ap), ac.dot(&ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return p.distance_to(&a);
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(&bp), ac.dot(&bp));
    if d3 >= 0.0 && d4 <= d3 {
        return p.distance_to(&b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let t = d1 / (d1 - d3);
        return p.distance_to(&(a + ab * t));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(&cp), ac.dot(&cp));
    if d6 >= 0.0 && d5 <= d6 {
        return p.distance_to(&c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let t = d2 / (d2 - d6);
        return p.distance_to(&(a + ac * t));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return p.distance_to(&(b + (c - b) * t));
    }
    let denom = 1.0 / (va + vb + vc);
    let (s, t) = (vb * denom, vc * denom);
    p.distance_to(&(a + ab * s + ac * t))
}

/// Triangles of a mesh bucketed on a uniform grid.
pub struct Grid<'a> {
    mesh: &'a TriMesh,
    cell: f64,
    buckets: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl<'a> Grid<'a> {
    pub fn new(mesh: &'a TriMesh, cell: f64) -> Grid<'a> {
        // A cell at least the mesh's longest sampled edge: a cell much
        // smaller than a triangle buckets it into millions of cells (a box
        // face or a cylinder strip spans the whole body).
        let n = mesh.triangles.len().max(1);
        let edge = mesh
            .triangles
            .iter()
            .step_by((n / 512).max(1))
            .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
            .map(|(i, j)| {
                mesh.vertices[i]
                    .position
                    .distance_to(&mesh.vertices[j].position)
            })
            .fold(0.0f64, f64::max);
        let cell = cell.max(edge).max(1e-12);
        let mut buckets: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        let key = |x: f64| (x / cell).floor() as i64;
        for (i, t) in mesh.triangles.iter().enumerate() {
            let ps: Vec<Point3> = t.iter().map(|&k| mesh.vertices[k].position).collect();
            let lo = |f: fn(&Point3) -> f64| ps.iter().map(f).fold(f64::INFINITY, f64::min);
            let hi = |f: fn(&Point3) -> f64| ps.iter().map(f).fold(f64::NEG_INFINITY, f64::max);
            for x in key(lo(|p| p.x))..=key(hi(|p| p.x)) {
                for y in key(lo(|p| p.y))..=key(hi(|p| p.y)) {
                    for z in key(lo(|p| p.z))..=key(hi(|p| p.z)) {
                        buckets.entry((x, y, z)).or_default().push(i);
                    }
                }
            }
        }
        Grid {
            mesh,
            cell,
            buckets,
        }
    }

    /// The distance from `p` to the mesh, if under `reach` (else `None`).
    pub fn distance_within(&self, p: Point3, reach: f64) -> Option<f64> {
        let key = |x: f64| (x / self.cell).floor() as i64;
        let r = (reach / self.cell).ceil() as i64;
        let (cx, cy, cz) = (key(p.x), key(p.y), key(p.z));
        let mut best = f64::INFINITY;
        for x in cx - r..=cx + r {
            for y in cy - r..=cy + r {
                for z in cz - r..=cz + r {
                    for &i in self.buckets.get(&(x, y, z)).into_iter().flatten() {
                        let t = self.mesh.triangles[i];
                        let q = |k: usize| self.mesh.vertices[t[k]].position;
                        best = best.min(point_triangle_distance(p, q(0), q(1), q(2)));
                    }
                }
            }
        }
        (best <= reach).then_some(best)
    }
}

/// The smallest vertex-to-surface distance between two meshes within
/// `reach`, skipping vertices within `exclude` of a point (the contact):
/// `None` when every vertex is farther than `reach`.
pub fn closest_approach(
    a: &TriMesh,
    b: &TriMesh,
    reach: f64,
    exclude: Option<(Point3, f64)>,
) -> Option<f64> {
    closest_approach_except(a, b, reach, exclude.as_slice())
}

/// [`closest_approach`] skipping vertices within any of several balls (a
/// contact curve's tube).
pub fn closest_approach_except(
    a: &TriMesh,
    b: &TriMesh,
    reach: f64,
    exclude: &[(Point3, f64)],
) -> Option<f64> {
    let (ga, gb) = (Grid::new(a, 2.0 * reach), Grid::new(b, 2.0 * reach));
    let skip = |p: &Point3| exclude.iter().any(|(c, r)| p.distance_to(c) <= *r);
    let mut best: Option<f64> = None;
    for (mesh, grid) in [(a, &gb), (b, &ga)] {
        for vtx in &mesh.vertices {
            if skip(&vtx.position) {
                continue;
            }
            if let Some(d) = grid.distance_within(vtx.position, reach) {
                best = Some(best.map_or(d, |x: f64| x.min(d)));
            }
        }
    }
    best
}

/// Whether `p` is inside a closed mesh: ray parity along three unrelated
/// directions, by majority (a ray through a shared vertex or edge, such as
/// a cap fan's centre, can miscount one parity).
pub fn inside(mesh: &TriMesh, p: Point3) -> bool {
    let dirs = [
        Vec3::new(0.5773, 0.5774, 0.5775),
        Vec3::new(-0.2672, 0.8017, 0.5345),
        Vec3::new(0.71, -0.1, 0.7),
    ];
    dirs.iter().filter(|d| ray_parity(mesh, p, **d)).count() >= 2
}

fn ray_parity(mesh: &TriMesh, p: Point3, d: Vec3) -> bool {
    let mut crossings = 0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        // Möller–Trumbore.
        let (e1, e2) = (b - a, c - a);
        let h = d.cross(&e2);
        let det = e1.dot(&h);
        if det.abs() < 1e-18 {
            continue;
        }
        let f = 1.0 / det;
        let s = p - a;
        let u = f * s.dot(&h);
        if !(0.0..=1.0).contains(&u) {
            continue;
        }
        let q = s.cross(&e1);
        let w = f * d.dot(&q);
        if w < 0.0 || u + w > 1.0 {
            continue;
        }
        if f * e2.dot(&q) > 0.0 {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

/// The distance from `p` to a polyline (closed when `closed`).
pub fn polyline_distance(p: Point3, line: &[Point3], closed: bool) -> f64 {
    let n = line.len();
    let segs = if closed { n } else { n.saturating_sub(1) };
    (0..segs)
        .map(|i| {
            let (a, b) = (line[i], line[(i + 1) % n]);
            let ab = b - a;
            let t = ((p - a).dot(&ab) / ab.dot(&ab).max(1e-300)).clamp(0.0, 1.0);
            p.distance_to(&(a + ab * t))
        })
        .fold(
            if n == 1 {
                p.distance_to(&line[0])
            } else {
                f64::INFINITY
            },
            f64::min,
        )
}

/// The axis-aligned box of a mesh's vertices (`popkern_measure`'s
/// `body_bounding_box` reads B-rep vertices only, PK-1258).
pub fn mesh_box(mesh: &TriMesh) -> (Point3, Point3) {
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for v in &mesh.vertices {
        let p = v.position;
        lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    (lo, hi)
}
