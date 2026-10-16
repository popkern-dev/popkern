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

use popkern_math::PortableMath;

use crate::crease_flip::{angle_between, key};
use crate::mesh::FaceMesh;
use crate::vertex_normal::{chart_singular, limit_normal, singular_without_unique_normal};

/// A candidate move: (leaning facets, worst tilt), chart point, and the
/// position and normal there.
type Move = ((usize, f64), (f64, f64), Point3, Vec3);

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
            // The star in one lift around `v`: across a periodic seam a
            // ring vertex can be stored a period away (tilt_refine::lift).
            let periods = surface.periodicity();
            let lifted = |w: usize| crate::tilt_refine::lift(mesh.vertices[w].uv, uv0, periods);
            let centroid = ring.iter().fold((0.0, 0.0), |a, &w| {
                let q = lifted(w);
                (a.0 + q.0 / ring.len() as f64, a.1 + q.1 / ring.len() as f64)
            });
            // The worst tilt and spread of `v`'s star with `v` at chart
            // point `uv`, or `None` when a facet would flip its chart
            // winding, lose its area or face against its normals.
            let evaluate = |uv: (f64, f64)| -> Option<(f64, f64, usize, Point3, Vec3)> {
                if chart_singular(surface, uv.0, uv.1) {
                    return None;
                }
                let n = limit_normal(surface, uv.0, uv.1)?;
                let (position, normal) = (
                    surface.point_at(uv.0, uv.1),
                    Vec3::new(sign * n.x, sign * n.y, sign * n.z),
                );
                let (mut tilt_max, mut spread_max, mut over) = (0.0f64, 0.0f64, 0usize);
                for &ti in star {
                    let t = mesh.triangles[ti];
                    let old_uv = t.map(lifted);
                    let new = t.map(|i| {
                        if i == v {
                            (position, normal, uv)
                        } else {
                            let w = &mesh.vertices[i];
                            (w.position, w.normal, lifted(i))
                        }
                    });
                    let (a_old, a_new) = (chart_area(old_uv), chart_area(new.map(|x| x.2)));
                    if a_new.signum() != a_old.signum() || a_new.abs() <= 1e-9 * a_old.abs() {
                        return None;
                    }
                    let f = (new[1].0 - new[0].0).cross(&(new[2].0 - new[0].0));
                    if new.iter().any(|x| f.dot(&x.1) <= 0.0) {
                        return None;
                    }
                    let x = tilt_of(new.map(|x| x.0), new.map(|x| x.1));
                    tilt_max = tilt_max.max(x);
                    over += usize::from(x > bound);
                    spread_max = spread_max.max(spread_of(new.map(|x| x.1)));
                }
                Some((tilt_max, spread_max, over, position, normal))
            };
            let (old_tilt, old_spread, old_over) =
                star.iter().fold((0.0f64, 0.0f64, 0usize), |acc, &ti| {
                    let t = mesh.triangles[ti];
                    let n = t.map(|i| mesh.vertices[i].normal);
                    let x = tilt_of(t.map(|i| mesh.vertices[i].position), n);
                    (
                        acc.0.max(x),
                        acc.1.max(spread_of(n)),
                        acc.2 + usize::from(x > bound),
                    )
                });
            // Candidates: the ring centroid, and a pattern of eight moves
            // at two radii (a quarter and half the ring's mean chart
            // reach): the centroid is a fixed point where the star's worst
            // facet is shaped by its boundary samples.
            let reach = ring
                .iter()
                .map(|&w| {
                    let q = lifted(w);
                    ((q.0 - uv0.0).powi(2) + (q.1 - uv0.1).powi(2)).sqrt()
                })
                .sum::<f64>()
                / ring.len() as f64;
            let mut candidates = vec![centroid];
            for r in [0.25 * reach, 0.5 * reach] {
                for k in 0..8 {
                    let a = std::f64::consts::FRAC_PI_4 * k as f64;
                    candidates.push((uv0.0 + r * a.pcos(), uv0.1 + r * a.psin()));
                }
            }
            // Fewer leaning facets first, then a lower worst tilt; never
            // more leaning facets than before (a lower worst can still
            // tip more of the star past the bound).
            let mut best: Option<Move> = None;
            for uv in candidates {
                if let Some((t, sp, over, position, normal)) = evaluate(uv)
                    && t < old_tilt
                    && over <= old_over
                    && sp <= bound.max(old_spread)
                    && best.is_none_or(|b| (over, t) < b.0)
                {
                    best = Some(((over, t), uv, position, normal));
                }
            }
            if let Some((_, uv, position, normal)) = best {
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
