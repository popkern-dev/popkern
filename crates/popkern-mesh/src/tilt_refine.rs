//! Obtuse-facet bisection for the facet-tilt contract (E4 step 2 of
//! docs/architecture/display-tessellation-in-face-folds.md).
//!
//! After the 3D Delaunay flips, the interior facets still leaning past
//! the bound are obtuse triangles at the vertex-normal spread limit
//! (coil tube: smallest angle ~6 deg, spread 11.8-12 deg): the flip that
//! would cure them joins vertices whose normals differ past the bound,
//! so the spread guard refuses it. Bisecting the triangle's 3D-longest
//! edge — the one opposite its obtuse angle — at the surface midpoint
//! removes the obtuse angle (Rivara) without joining far-apart normals.
//! The edge must be interior (boundary edges carry the shared canonical
//! samples) and the split must hold: every new facet faces along its
//! vertex normals, and the four new facets' worst tilt and spread are
//! no worse than the pair's. Splits never cascade (no LEPP): a refused
//! split leaves the facet for the census.

use std::collections::HashMap;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

use crate::crease_flip::{EdgeKey, angle_between, key, spread};
use crate::mesh::{FaceMesh, MeshVertex};
use crate::vertex_normal::{chart_singular, limit_normal, singular_without_unique_normal};

/// Bisection rounds; each splits every leaning interior facet at most
/// once.
const ROUNDS: usize = 4;

fn tilt_of(p: [Point3; 3], n: [Vec3; 3]) -> f64 {
    let f = (p[1] - p[0]).cross(&(p[2] - p[0]));
    n.iter().map(|&v| angle_between(f, v)).fold(0.0, f64::max)
}

fn tilt(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    tilt_of(
        t.map(|i| mesh.vertices[i].position),
        t.map(|i| mesh.vertices[i].normal),
    )
}

/// Bisect the 3D-longest edge of interior facets leaning past `bound`
/// (see the module docs). Returns the triangles created or rewritten,
/// for a local flip pass.
pub(crate) fn bisect_leaning_facets(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    bound: f64,
) -> Vec<usize> {
    if matches!(surface, SurfaceGeometry::Plane { .. }) {
        return Vec::new();
    }
    // The face's normal sense, recovered once (as angular_refine does).
    let sign = mesh
        .vertices
        .iter()
        .find_map(|v| {
            let d = surface.normal_at(v.uv.0, v.uv.1).dot(&v.normal);
            (d.abs() > 1e-6).then_some(if d > 0.0 { 1.0 } else { -1.0 })
        })
        .unwrap_or(1.0);
    let mut touched = Vec::new();
    for _ in 0..ROUNDS {
        let mut edge_tris: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
        for (ti, t) in mesh.triangles.iter().enumerate() {
            for k in 0..3 {
                edge_tris
                    .entry(key(t[k], t[(k + 1) % 3]))
                    .or_default()
                    .push(ti);
            }
        }
        let singular = |mesh: &FaceMesh, v: usize| {
            let uv = mesh.vertices[v].uv;
            singular_without_unique_normal(surface, uv.0, uv.1)
        };
        let leaning: Vec<usize> = (0..mesh.triangles.len())
            .filter(|&ti| {
                let t = mesh.triangles[ti];
                !t.iter().any(|&v| singular(mesh, v)) && tilt(mesh, t) > bound
            })
            .collect();
        let mut done: std::collections::HashSet<usize> = Default::default();
        let mut split_any = false;
        for ti in leaning {
            if done.contains(&ti) {
                continue;
            }
            let t = mesh.triangles[ti];
            let p = |i: usize| mesh.vertices[i].position;
            let (k, _) = (0..3)
                .map(|k| (k, p(t[k]).distance_to(&p(t[(k + 1) % 3]))))
                .fold((0, 0.0), |b, x| if x.1 > b.1 { x } else { b });
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let Some(&[x, y]) = edge_tris.get(&key(a, b)).map(Vec::as_slice) else {
                continue; // boundary edge
            };
            let tn = if x == ti { y } else { x };
            if done.contains(&tn)
                || singular(
                    mesh,
                    mesh.triangles[tn]
                        .iter()
                        .copied()
                        .find(|&v| v != a && v != b)
                        .unwrap_or(a),
                )
            {
                continue;
            }
            let (va, vb) = (&mesh.vertices[a], &mesh.vertices[b]);
            let uv = (0.5 * (va.uv.0 + vb.uv.0), 0.5 * (va.uv.1 + vb.uv.1));
            let position = surface.point_at(uv.0, uv.1);
            let len = va.position.distance_to(&vb.position);
            let chord_mid = Point3::new(
                0.5 * (va.position.x + vb.position.x),
                0.5 * (va.position.y + vb.position.y),
                0.5 * (va.position.z + vb.position.z),
            );
            if chart_singular(surface, uv.0, uv.1) || position.distance_to(&chord_mid) > 0.5 * len {
                continue;
            }
            let Some(n) = limit_normal(surface, uv.0, uv.1) else {
                continue;
            };
            let normal = Vec3::new(sign * n.x, sign * n.y, sign * n.z);
            let m = MeshVertex {
                position,
                normal,
                uv,
            };
            // Plan both halves of both triangles; land only if all hold.
            let mut halves: Vec<(usize, [usize; 3], [usize; 3])> = Vec::new();
            let (mut old_tilt, mut old_spread, mut new_tilt, mut new_spread) =
                (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            let mut ok = true;
            for tj in [ti, tn] {
                let tri = mesh.triangles[tj];
                let Some(r) = (0..3).find(|&q| key(tri[q], tri[(q + 1) % 3]) == key(a, b)) else {
                    ok = false;
                    break;
                };
                let (pp, qq, cc) = (tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]);
                old_tilt = old_tilt.max(tilt(mesh, tri));
                old_spread = old_spread.max(spread(mesh, tri));
                for (h0, h1, h2) in [(pp, usize::MAX, cc), (usize::MAX, qq, cc)] {
                    let corner = |i: usize| {
                        if i == usize::MAX {
                            (m.position, m.normal)
                        } else {
                            (mesh.vertices[i].position, mesh.vertices[i].normal)
                        }
                    };
                    let c = [corner(h0), corner(h1), corner(h2)];
                    let f = (c[1].0 - c[0].0).cross(&(c[2].0 - c[0].0));
                    if c.iter().any(|x| f.dot(&x.1) <= 0.0) {
                        ok = false;
                    }
                    new_tilt = new_tilt.max(tilt_of(c.map(|x| x.0), c.map(|x| x.1)));
                    new_spread = new_spread.max(
                        angle_between(c[0].1, c[1].1)
                            .max(angle_between(c[1].1, c[2].1))
                            .max(angle_between(c[0].1, c[2].1)),
                    );
                }
                halves.push((tj, [pp, usize::MAX, cc], [usize::MAX, qq, cc]));
            }
            if !ok || new_tilt >= old_tilt || new_spread > bound.max(old_spread) {
                continue;
            }
            let mv = mesh.vertices.len();
            mesh.vertices.push(m);
            for (tj, h0, h1) in halves {
                let fill = |h: [usize; 3]| h.map(|i| if i == usize::MAX { mv } else { i });
                mesh.triangles[tj] = fill(h0);
                touched.push(tj);
                touched.push(mesh.triangles.len());
                mesh.triangles.push(fill(h1));
                done.insert(tj);
            }
            split_any = true;
        }
        if !split_any {
            break;
        }
    }
    touched.sort_unstable();
    touched.dedup();
    touched
}
