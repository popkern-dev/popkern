//! Tangential relaxation of leaning facets' interior vertices (E3 of
//! docs/architecture/display-tessellation-in-face-folds.md).
//!
//! Some leaning facets are flat caps made by vertex PLACEMENT: an
//! interior lattice point sits almost on the line between a neighbour
//! and a boundary sample (propeller trailing edge: uv rows at v = 0.890
//! against a boundary sample at 0.891, angles 157-174 deg). Flips join
//! far-apart normals and fans keep the bad point; moving it does not.
//! Each interior, non-singular vertex of a leaning facet moves to the
//! chart centroid of its one-ring (evaluated on the surface, with the
//! face-sense limit normal) when every facet of its star keeps its chart
//! winding with area and faces along its vertex normals, the star's
//! spread stays within the bound (or its old spread), and the star's
//! worst tilt drops. Boundary vertices never move: they are the
//! canonical samples both faces share.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

use crate::crease_flip::{angle_between, key};
use crate::mesh::FaceMesh;
use crate::vertex_normal::{chart_singular, limit_normal, singular_without_unique_normal};

/// Relaxation rounds.
const ROUNDS: usize = 3;

fn tilt_of(p: [Point3; 3], n: [Vec3; 3]) -> f64 {
    let f = (p[1] - p[0]).cross(&(p[2] - p[0]));
    n.iter().map(|&v| angle_between(f, v)).fold(0.0, f64::max)
}

fn spread_of(n: [Vec3; 3]) -> f64 {
    angle_between(n[0], n[1])
        .max(angle_between(n[1], n[2]))
        .max(angle_between(n[0], n[2]))
}

fn chart_area(p: [(f64, f64); 3]) -> f64 {
    (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1)
}

/// Relax the interior vertices of facets leaning past `bound` (see the
/// module docs), starting from the facets `seeds` (`None`: all).
/// `sign` is the face's normal sense. Returns the
/// triangles around the vertices that moved.
pub(crate) fn relax_leaning_vertices(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    sign: f64,
    bound: f64,
    seeds: Option<&[usize]>,
) -> Vec<usize> {
    let mut moved: Vec<usize> = Vec::new();
    // Moving a vertex keeps the topology: build it once.
    let mut uses: HashMap<(usize, usize), usize> = HashMap::new();
    let mut incident: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            *uses.entry(key(t[k], t[(k + 1) % 3])).or_default() += 1;
            incident.entry(t[k]).or_default().push(ti);
        }
    }
    let boundary: HashSet<usize> = uses
        .iter()
        .filter(|(_, n)| **n == 1)
        .flat_map(|(e, _)| [e.0, e.1])
        .collect();
    // Round 0 scans every facet; later rounds the stars of moved vertices.
    let mut scope: Option<Vec<usize>> = seeds.map(<[usize]>::to_vec);
    for _ in 0..ROUNDS {
        let tilt = |mesh: &FaceMesh, t: [usize; 3]| {
            tilt_of(
                t.map(|i| mesh.vertices[i].position),
                t.map(|i| mesh.vertices[i].normal),
            )
        };
        let candidates: Vec<usize> = {
            let tris: Vec<usize> = match &scope {
                None => (0..mesh.triangles.len()).collect(),
                Some(s) => s.clone(),
            };
            let mut c: Vec<usize> = tris
                .into_iter()
                .map(|ti| mesh.triangles[ti])
                .filter(|&t| tilt(mesh, t) > bound)
                .flat_map(|t| t.into_iter())
                .filter(|v| !boundary.contains(v))
                .collect();
            c.sort_unstable();
            c.dedup();
            c
        };
        let mut moved_round = 0usize;
        let mut next: Vec<usize> = Vec::new();
        for v in candidates {
            let Some(star) = incident.get(&v) else {
                continue;
            };
            let uv0 = mesh.vertices[v].uv;
            if singular_without_unique_normal(surface, uv0.0, uv0.1) {
                continue;
            }
            let ring: HashSet<usize> = star
                .iter()
                .flat_map(|&ti| mesh.triangles[ti])
                .filter(|&w| w != v)
                .collect();
            if ring.is_empty() {
                continue;
            }
            let uv = ring.iter().fold((0.0, 0.0), |a, &w| {
                let q = mesh.vertices[w].uv;
                (a.0 + q.0 / ring.len() as f64, a.1 + q.1 / ring.len() as f64)
            });
            if chart_singular(surface, uv.0, uv.1) {
                continue;
            }
            let Some(n) = limit_normal(surface, uv.0, uv.1) else {
                continue;
            };
            let (position, normal) = (
                surface.point_at(uv.0, uv.1),
                Vec3::new(sign * n.x, sign * n.y, sign * n.z),
            );
            let corner = |i: usize| {
                if i == v {
                    (position, normal, uv)
                } else {
                    let w = &mesh.vertices[i];
                    (w.position, w.normal, w.uv)
                }
            };
            let (mut old_tilt, mut old_spread, mut new_tilt, mut new_spread) =
                (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            let ok = star.iter().all(|&ti| {
                let t = mesh.triangles[ti];
                let old = t.map(|i| {
                    (
                        mesh.vertices[i].position,
                        mesh.vertices[i].normal,
                        mesh.vertices[i].uv,
                    )
                });
                let new = t.map(corner);
                old_tilt = old_tilt.max(tilt_of(old.map(|x| x.0), old.map(|x| x.1)));
                old_spread = old_spread.max(spread_of(old.map(|x| x.1)));
                let (a_old, a_new) = (chart_area(old.map(|x| x.2)), chart_area(new.map(|x| x.2)));
                if a_new.signum() != a_old.signum() || a_new.abs() <= 1e-9 * a_old.abs() {
                    return false;
                }
                let f = (new[1].0 - new[0].0).cross(&(new[2].0 - new[0].0));
                if new.iter().any(|x| f.dot(&x.1) <= 0.0) {
                    return false;
                }
                new_tilt = new_tilt.max(tilt_of(new.map(|x| x.0), new.map(|x| x.1)));
                new_spread = new_spread.max(spread_of(new.map(|x| x.1)));
                true
            });
            if ok && new_tilt < old_tilt && new_spread <= bound.max(old_spread) {
                let w = &mut mesh.vertices[v];
                (w.position, w.normal, w.uv) = (position, normal, uv);
                moved_round += 1;
                next.extend(star.iter().copied());
            }
        }
        if moved_round == 0 {
            break;
        }
        next.sort_unstable();
        next.dedup();
        moved.extend(next.iter().copied());
        scope = Some(next);
    }
    moved.sort_unstable();
    moved.dedup();
    moved
}
