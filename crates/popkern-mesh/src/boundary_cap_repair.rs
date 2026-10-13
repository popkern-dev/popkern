//! Local repairs of leaning boundary caps that a chord split alone
//! cannot fix (see the parent module): flips that hand a cap a farther
//! apex, and collapsing a cap-zone apex onto the new boundary vertex.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceGeometry;

use super::{boundary_edges, chart_area, spread_of, tilt_of, triangle_spread, triangle_tilt};
use crate::mesh::{FaceMesh, MeshVertex};
use crate::vertex_normal::singular_without_unique_normal;

/// A cap flip candidate: worst new tilt, the neighbour triangle, the
/// two new triangles and the edge they replace.
pub(super) type FlipChoice = (f64, usize, [usize; 3], [usize; 3], (usize, usize));

/// Placeholder index of the vertex a planned split inserts.
pub(super) const NEW: usize = usize::MAX;

/// Plan collapsing apex `c` of cap `cap` onto a new boundary vertex `m`
/// (normal `n`): every other triangle around `c` takes `m` instead. Kept
/// only when no rewritten triangle flips its chart winding, loses its
/// area or faces against its normals, and the star's worst tilt drops.
pub(super) fn plan_apex_collapse(
    mesh: &FaceMesh,
    surface: &SurfaceGeometry,
    cap: usize,
    c: usize,
    m: &MeshVertex,
    topo: &MeshTopology,
    bound: f64,
) -> Option<Vec<(usize, [usize; 3])>> {
    let uv = mesh.vertices[c].uv;
    if topo.boundary.contains(&c) || singular_without_unique_normal(surface, uv.0, uv.1) {
        return None;
    }
    let star: Vec<usize> = topo
        .incident
        .get(&c)
        .map(|v| v.iter().copied().filter(|&ti| ti != cap).collect())
        .unwrap_or_default();
    let mut old = triangle_tilt(mesh, surface, mesh.triangles[cap])?;
    let (mut new, mut out) = (0.0f64, Vec::with_capacity(star.len()));
    let corner = |i: usize| {
        if i == NEW {
            (m.position, m.normal, m.uv)
        } else {
            let v = &mesh.vertices[i];
            (v.position, v.normal, v.uv)
        }
    };
    for ti in star {
        let t = mesh.triangles[ti];
        old = old.max(triangle_tilt(mesh, surface, t)?);
        let r = t.map(|i| if i == c { NEW } else { i });
        let q = r.map(corner);
        let area = |p: [(f64, f64); 3]| {
            0.5 * ((p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1))
        };
        let (a_old, a_new) = (chart_area(mesh, t), area(q.map(|x| x.2)));
        if a_new.signum() != a_old.signum() || a_new.abs() <= 1e-9 * a_old.abs() {
            return None;
        }
        let tilt = tilt_of(q.map(|x| x.0), q.map(|x| x.1))?;
        if spread_of(q.map(|x| x.1)) > bound.max(triangle_spread(mesh, t)) {
            return None;
        }
        let f = (q[1].0 - q[0].0).cross(&(q[2].0 - q[0].0));
        if q.iter().any(|x| f.dot(&x.1) <= 0.0) {
            return None;
        }
        new = new.max(tilt);
        out.push((ti, r));
    }
    (new < old).then_some(out)
}

/// A face mesh's boundary vertices and vertex -> triangles map, built
/// once per round for a mesh that needs a collapse and patched by every
/// split and collapse the round lands on that mesh.
pub(super) struct MeshTopology {
    pub(super) boundary: std::collections::HashSet<usize>,
    pub(super) incident: HashMap<usize, Vec<usize>>,
}

impl MeshTopology {
    pub(super) fn of(mesh: &FaceMesh) -> Self {
        let boundary = boundary_edges(mesh)
            .into_iter()
            .flat_map(|(a, b, _)| [a, b])
            .collect();
        let mut incident: HashMap<usize, Vec<usize>> = HashMap::new();
        for (ti, t) in mesh.triangles.iter().enumerate() {
            for &v in t {
                incident.entry(v).or_default().push(ti);
            }
        }
        Self { boundary, incident }
    }

    pub(super) fn unlink(&mut self, v: usize, t: usize) {
        if let Some(list) = self.incident.get_mut(&v) {
            list.retain(|&x| x != t);
        }
    }
}

/// Flip the interior edges of leaning boundary facets when the flip
/// lowers the pair's worst tilt: a cap whose apex sits close beside
/// its boundary chord takes the neighbour's farther apex instead (the
/// chord split cannot help it: both halves keep the near apex). Guards:
/// both new triangles keep the old pair's chart winding with area, and
/// face along their vertex normals. Returns whether anything flipped.
pub(super) fn flip_leaning_caps(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    bound: f64,
) -> bool {
    let mut edge_tris: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            edge_tris.entry((a.min(b), a.max(b))).or_default().push(ti);
        }
    }
    let leaning: Vec<usize> = (0..mesh.triangles.len())
        .filter(|&ti| {
            let t = mesh.triangles[ti];
            (0..3).any(|k| {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                edge_tris
                    .get(&(a.min(b), a.max(b)))
                    .is_some_and(|v| v.len() == 1)
            }) && triangle_tilt(mesh, surface, t).is_some_and(|x| x > bound)
        })
        .collect();
    let mut any = false;
    for t0 in leaning {
        let tri0 = mesh.triangles[t0];
        let Some(old0) = triangle_tilt(mesh, surface, tri0) else {
            continue;
        };
        if old0 <= bound {
            continue; // fixed by an earlier flip this pass
        }
        let mut best: Option<FlipChoice> = None;
        for k in 0..3 {
            let (a, b) = (tri0[k], tri0[(k + 1) % 3]);
            let e = (a.min(b), a.max(b));
            let Some(&[x, y]) = edge_tris.get(&e).map(Vec::as_slice) else {
                continue;
            };
            let t1 = if x == t0 { y } else { x };
            let tri1 = mesh.triangles[t1];
            let c = tri0[(k + 2) % 3];
            let Some(d) = tri1.iter().copied().find(|&v| v != a && v != b) else {
                continue;
            };
            if edge_tris.contains_key(&(c.min(d), c.max(d))) {
                continue;
            }
            // tri0 runs a -> b; new pair (a, d, c) and (d, b, c).
            let (n0, n1) = ([a, d, c], [d, b, c]);
            let sense = chart_area(mesh, tri0).signum();
            let span = chart_area(mesh, tri0).abs() + chart_area(mesh, tri1).abs();
            if [n0, n1].iter().any(|&t| {
                let ar = chart_area(mesh, t);
                ar.signum() != sense || ar.abs() <= 1e-9 * span
            }) {
                continue;
            }
            let old = old0.max(triangle_tilt(mesh, surface, tri1).unwrap_or(f64::INFINITY));
            let (Some(x0), Some(x1)) = (
                triangle_tilt(mesh, surface, n0),
                triangle_tilt(mesh, surface, n1),
            ) else {
                continue;
            };
            let new = x0.max(x1);
            let faces_ok = [n0, n1].iter().all(|&t| {
                let p = t.map(|i| mesh.vertices[i].position);
                let f = (p[1] - p[0]).cross(&(p[2] - p[0]));
                t.iter().all(|&i| f.dot(&mesh.vertices[i].normal) > 0.0)
            });
            let spread_ok = triangle_spread(mesh, n0).max(triangle_spread(mesh, n1))
                <= bound.max(triangle_spread(mesh, tri0).max(triangle_spread(mesh, tri1)));
            if faces_ok && spread_ok && new < old && best.is_none_or(|b| new < b.0) {
                best = Some((new, t1, n0, n1, e));
            }
        }
        if let Some((_, t1, n0, n1, e)) = best {
            mesh.triangles[t0] = n0;
            mesh.triangles[t1] = n1;
            // Rebuild the touched adjacency lazily: a full rebuild per
            // flip is quadratic, so later candidates re-check via the
            // tilt guard and the stale-edge checks above.
            edge_tris.remove(&e);
            for (t, tri) in [(t0, n0), (t1, n1)] {
                for k in 0..3 {
                    let (a, b) = (tri[k], tri[(k + 1) % 3]);
                    let list = edge_tris.entry((a.min(b), a.max(b))).or_default();
                    list.retain(|&o| o != t0 && o != t1);
                    list.push(t);
                }
            }
            any = true;
        }
    }
    any
}
