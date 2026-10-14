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

/// The facet-tilt finish of one face: 3D Delaunay flips, then the
/// fan bisection of facets still leaning, then flips again.
pub(crate) fn finish_tilt(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    chord_height: f64,
    bound: f64,
) {
    finish_tilt_from(mesh, surface, chord_height, bound, None);
}

/// [`finish_tilt`] around the triangles `seeds` only (`None`: the whole
/// face): after a local change (the cross-face boundary splits), only
/// its neighbourhood needs the passes.
pub(crate) fn finish_tilt_from(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    chord_height: f64,
    bound: f64,
    seeds: Option<&[usize]>,
) {
    use crate::surface_delaunay_flip::delaunay_flips_3d_from;
    delaunay_flips_3d_from(mesh, surface, chord_height, bound, seeds);
    let split = bisect_leaning_facets_from(mesh, surface, bound, seeds);
    if !split.is_empty() {
        delaunay_flips_3d_from(mesh, surface, chord_height, bound, Some(&split));
    }
    let sign = face_sense(mesh, surface);
    let moved = crate::tilt_relax::relax_leaning_vertices(mesh, surface, sign, bound, seeds);
    if !moved.is_empty() {
        delaunay_flips_3d_from(mesh, surface, chord_height, bound, Some(&moved));
    }
}

/// The face's normal sense: the sign relating stored vertex normals to
/// the chart normal.
fn face_sense(mesh: &FaceMesh, surface: &SurfaceGeometry) -> f64 {
    mesh.vertices
        .iter()
        .find_map(|v| {
            let d = surface.normal_at(v.uv.0, v.uv.1).dot(&v.normal);
            (d.abs() > 1e-6).then_some(if d > 0.0 { 1.0 } else { -1.0 })
        })
        .unwrap_or(1.0)
}

/// A planned fan: worst new tilt, the new vertex, and per old triangle
/// its two halves (`usize::MAX` marks the new vertex).
type FanPlan = (f64, MeshVertex, Vec<(usize, [usize; 3], [usize; 3])>);

/// Plan inserting a vertex at chart point `uv` into the quad of
/// triangles `(t0, t1)` sharing edge `(a, b)`, as the fan of the quad's
/// four corners around it. `Some((worst new tilt, vertex, halves))` when
/// every new facet keeps its triangle's chart winding with area and
/// faces along its vertex normals, the spread stays within the bound
/// (or the old pair's), the worst tilt drops and no more facets lean
/// past the bound than before.
fn plan_fan(
    mesh: &FaceMesh,
    surface: &SurfaceGeometry,
    sign: f64,
    bound: f64,
    (t0, t1): (usize, usize),
    (a, b): (usize, usize),
    uv: (f64, f64),
) -> Option<FanPlan> {
    if chart_singular(surface, uv.0, uv.1) {
        return None;
    }
    let n = limit_normal(surface, uv.0, uv.1)?;
    let m = MeshVertex {
        position: surface.point_at(uv.0, uv.1),
        normal: Vec3::new(sign * n.x, sign * n.y, sign * n.z),
        uv,
    };
    // The image must sit near the quad (a chart seam or singularity
    // crossing lands it far away).
    let (va, vb) = (&mesh.vertices[a], &mesh.vertices[b]);
    let len = va.position.distance_to(&vb.position);
    let chord_mid = Point3::new(
        0.5 * (va.position.x + vb.position.x),
        0.5 * (va.position.y + vb.position.y),
        0.5 * (va.position.z + vb.position.z),
    );
    if m.position.distance_to(&chord_mid) > len {
        return None;
    }
    let mut halves = Vec::with_capacity(2);
    let (mut old_tilt, mut old_spread, mut new_tilt, mut new_spread) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut old_over, mut new_over) = (0usize, 0usize);
    for tj in [t0, t1] {
        let tri = mesh.triangles[tj];
        let r = (0..3).find(|&q| key(tri[q], tri[(q + 1) % 3]) == key(a, b))?;
        let (pp, qq, cc) = (tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]);
        old_tilt = old_tilt.max(tilt(mesh, tri));
        old_spread = old_spread.max(spread(mesh, tri));
        old_over += usize::from(tilt(mesh, tri) > bound);
        let sense = {
            let p = tri.map(|i| mesh.vertices[i].uv);
            ((p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1)).signum()
        };
        for h in [[pp, usize::MAX, cc], [usize::MAX, qq, cc]] {
            let corner = |i: usize| {
                if i == usize::MAX {
                    (m.position, m.normal, m.uv)
                } else {
                    let v = &mesh.vertices[i];
                    (v.position, v.normal, v.uv)
                }
            };
            let c = h.map(corner);
            let area = (c[1].2.0 - c[0].2.0) * (c[2].2.1 - c[0].2.1)
                - (c[2].2.0 - c[0].2.0) * (c[1].2.1 - c[0].2.1);
            if area.signum() != sense || area == 0.0 {
                return None;
            }
            let f = (c[1].0 - c[0].0).cross(&(c[2].0 - c[0].0));
            if c.iter().any(|x| f.dot(&x.1) <= 0.0) {
                return None;
            }
            let x = tilt_of(c.map(|x| x.0), c.map(|x| x.1));
            new_tilt = new_tilt.max(x);
            new_over += usize::from(x > bound);
            new_spread = new_spread.max(
                angle_between(c[0].1, c[1].1)
                    .max(angle_between(c[1].1, c[2].1))
                    .max(angle_between(c[0].1, c[2].1)),
            );
        }
        halves.push((tj, [pp, usize::MAX, cc], [usize::MAX, qq, cc]));
    }
    // Lower the worst tilt without adding facets past the bound (four
    // facets replace two: a lower worst could still lean more often).
    (new_tilt < old_tilt && new_over <= old_over && new_spread <= bound.max(old_spread))
        .then_some((new_tilt, m, halves))
}

/// Bisect the 3D-longest edge of interior facets leaning past `bound`
/// (see the module docs). Returns the triangles created or rewritten,
/// for a local flip pass.
pub(crate) fn bisect_leaning_facets_from(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    bound: f64,
    seeds: Option<&[usize]>,
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
    let mut next_scope: Option<Vec<usize>> = seeds.map(<[usize]>::to_vec);
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
        // Round 0 scans every facet; later rounds only the facets the
        // previous round created or rewrote (a facet whose candidates all
        // failed fails again while its quad is unchanged).
        let scope: Vec<usize> = match &next_scope {
            None => (0..mesh.triangles.len()).collect(),
            Some(s) => s.clone(),
        };
        let leaning: Vec<usize> = scope
            .into_iter()
            .filter(|&ti| {
                let t = mesh.triangles[ti];
                !t.iter().any(|&v| singular(mesh, v)) && tilt(mesh, t) > bound
            })
            .collect();
        let round_start = touched.len();
        let mut done: std::collections::HashSet<usize> = Default::default();
        let mut split_any = false;
        for ti in leaning {
            if done.contains(&ti) {
                continue;
            }
            let t = mesh.triangles[ti];
            // Every interior edge of the facet, each with three candidate
            // positions for the new vertex m joined to the quad's four
            // corners: the shared edge's midpoint (bisection), the
            // opposite apexes' midpoint (the flip the spread guard refused,
            // halved), and the quad's centre. The lowest new tilt wins.
            let uv_of = |i: usize| mesh.vertices[i].uv;
            let mid = |p: (f64, f64), q: (f64, f64)| (0.5 * (p.0 + q.0), 0.5 * (p.1 + q.1));
            let mut best: Option<FanPlan> = None;
            for k in 0..3 {
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let Some(&[x, y]) = edge_tris.get(&key(a, b)).map(Vec::as_slice) else {
                    continue; // boundary edge
                };
                let tn = if x == ti { y } else { x };
                let Some(d) = mesh.triangles[tn]
                    .iter()
                    .copied()
                    .find(|&v| v != a && v != b)
                else {
                    continue;
                };
                if done.contains(&tn) || singular(mesh, d) {
                    continue;
                }
                let (pa, pb, pc, pd) = (uv_of(a), uv_of(b), uv_of(c), uv_of(d));
                let centre = (
                    0.25 * (pa.0 + pb.0 + pc.0 + pd.0),
                    0.25 * (pa.1 + pb.1 + pc.1 + pd.1),
                );
                for uv in [mid(pa, pb), mid(pc, pd), centre] {
                    if let Some((score, m, halves)) =
                        plan_fan(mesh, surface, sign, bound, (ti, tn), (a, b), uv)
                        && best.as_ref().is_none_or(|x| score < x.0)
                    {
                        best = Some((score, m, halves));
                    }
                }
            }
            let Some((_, m, halves)) = best else {
                continue;
            };
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
        let mut scope: Vec<usize> = touched[round_start..].to_vec();
        scope.sort_unstable();
        scope.dedup();
        next_scope = Some(scope);
    }
    touched.sort_unstable();
    touched.dedup();
    touched
}
