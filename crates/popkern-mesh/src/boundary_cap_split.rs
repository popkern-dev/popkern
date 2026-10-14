//! Cross-face boundary split for leaning boundary facets (E1 of
//! docs/architecture/display-tessellation-in-face-folds.md).
//!
//! The display contract bounds every facet's tilt: the angle between its
//! flat normal and the surface normal at each of its vertices. A facet
//! hanging off a long boundary chord with its third vertex close beside
//! it (a coil's ground cut, a pyramid rim, a thread runout corner) leans
//! toward the chord's plane however its vertex normals agree, and no
//! interior split or flip can change it: the boundary chord is a
//! canonical edge sample pair both incident faces share. This pass runs
//! on the finished, oriented face meshes and splits such a chord at the
//! B-rep EDGE curve's midpoint in EVERY face mesh that carries it (the
//! copies are found by bit-identical end positions), so the weld holds
//! by construction. A split is transactional: it lands only when no
//! copy's new triangles lean past the bound by more than the old
//! triangle did; afterwards the touched faces re-run the guarded
//! refinement and crease flips locally.

use std::collections::HashMap;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use popkern_math::PortableMath;

use crate::mesh::{FaceMesh, MeshVertex};

#[path = "boundary_cap_curve.rs"]
mod curve;
use curve::{CurveSamples, edge_curve_midpoint};
#[path = "boundary_cap_repair.rs"]
mod repair;
use crate::vertex_normal::{limit_normal, singular_without_unique_normal};
use repair::{MeshTopology, NEW, flip_leaning_caps, plan_apex_collapse};

/// Rounds: each flips the leaning caps it can, or else splits each
/// requested chord once.
const ROUNDS: usize = 12;

/// Tilt a leaning cap's split halves may gain (radians, 2 deg).
const STEP_SLACK: f64 = 2.0 * std::f64::consts::PI / 180.0;

type PosKey = [u64; 3];
type EdgeKey = (PosKey, PosKey);

fn pos_key(p: Point3) -> PosKey {
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

fn edge_key(a: Point3, b: Point3) -> EdgeKey {
    let (ka, kb) = (pos_key(a), pos_key(b));
    if ka <= kb { (ka, kb) } else { (kb, ka) }
}

fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let (la, lb) = (a.length(), b.length());
    if la < 1e-300 || lb < 1e-300 {
        return 0.0;
    }
    (a.dot(&b) / (la * lb)).clamp(-1.0, 1.0).pacos()
}

/// Largest angle between a facet with corners `p` and its vertex
/// normals `n`; `None` for a zero-area facet.
fn tilt_of(p: [Point3; 3], n: [Vec3; 3]) -> Option<f64> {
    let f = (p[1] - p[0]).cross(&(p[2] - p[0]));
    (f.length() > 1e-300).then(|| n.iter().map(|&v| angle_between(f, v)).fold(0.0, f64::max))
}

/// Largest pairwise angle between vertex normals `n` (the display
/// preset's per-triangle spread bound, which every repair must keep).
fn spread_of(n: [Vec3; 3]) -> f64 {
    angle_between(n[0], n[1])
        .max(angle_between(n[1], n[2]))
        .max(angle_between(n[0], n[2]))
}

/// Spread of triangle `t` of `mesh`.
fn triangle_spread(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    spread_of(t.map(|i| mesh.vertices[i].normal))
}

/// Tilt of triangle `t` of `mesh`, `None` when a vertex is chart-
/// singular (no unique normal) or the facet has no area.
fn triangle_tilt(mesh: &FaceMesh, surface: &SurfaceGeometry, t: [usize; 3]) -> Option<f64> {
    if t.iter().any(|&i| {
        let uv = mesh.vertices[i].uv;
        singular_without_unique_normal(surface, uv.0, uv.1)
    }) {
        return None;
    }
    tilt_of(
        t.map(|i| mesh.vertices[i].position),
        t.map(|i| mesh.vertices[i].normal),
    )
}

/// Boundary edges of one face mesh (used by one triangle), as
/// `(local a, local b, triangle)` with `a -> b` in the triangle's winding.
fn boundary_edges(mesh: &FaceMesh) -> Vec<(usize, usize, usize)> {
    let mut count: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *count.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let mut out = Vec::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if count.get(&(a.min(b), a.max(b))) == Some(&1) {
                out.push((a, b, ti));
            }
        }
    }
    out
}

/// Signed chart area of triangle `t` (its uv winding).
fn chart_area(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    let p = t.map(|i| mesh.vertices[i].uv);
    0.5 * ((p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1))
}

/// One planned copy of a split: face mesh index, triangle, the edge's
/// local ends in winding order, and the inserted vertex.
struct Copy {
    mesh: usize,
    triangle: usize,
    a: usize,
    b: usize,
    vertex: MeshVertex,
    /// Collapse the cap's apex onto the new vertex: the apex's other
    /// triangles, rewritten with `NEW` for the apex (the cap goes).
    collapse: Option<Vec<(usize, [usize; 3])>>,
}

/// Split leaning boundary chords across faces until no boundary facet
/// tilts past `bound`, or nothing more can be split. Returns the face
/// meshes it changed. `surfaces[i]` is face mesh `i`'s surface.
pub(crate) fn split_leaning_boundary_caps(
    brep: &BRep,
    meshes: &mut [FaceMesh],
    bound: f64,
) -> Vec<usize> {
    let gate = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let mut changed: Vec<usize> = Vec::new();
    let mut curve_cache = CurveSamples::new();
    for _ in 0..ROUNDS {
        // Every boundary edge copy, keyed by its end positions.
        let mut copies: HashMap<EdgeKey, Vec<(usize, usize, usize, usize)>> = HashMap::new();
        let mut requests: Vec<(f64, EdgeKey)> = Vec::new();
        let mut leaning: Vec<usize> = Vec::new();
        for (mi, fm) in meshes.iter().enumerate() {
            let Some(surface) = brep
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            else {
                continue;
            };
            let edges = boundary_edges(fm);
            // Every boundary edge of a leaning facet is a candidate,
            // longest first: a facet whose longest chord cannot split
            // (its halves would lean more, the hex bolt's runout corner)
            // may split along its other boundary edge.
            let mut any_leaning = false;
            for &(a, b, ti) in &edges {
                let (pa, pb) = (fm.vertices[a].position, fm.vertices[b].position);
                let key = edge_key(pa, pb);
                copies.entry(key).or_default().push((mi, ti, a, b));
                if triangle_tilt(fm, surface, fm.triangles[ti]).is_some_and(|t| t > bound) {
                    requests.push((pa.distance_to(&pb), key));
                    any_leaning = true;
                }
            }
            if any_leaning {
                leaning.push(mi);
            }
        }
        requests.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
        let mut requests: Vec<EdgeKey> = requests.into_iter().map(|(_, k)| k).collect();
        let mut seen = std::collections::HashSet::new();
        requests.retain(|k| seen.insert(*k));
        if requests.is_empty() {
            break;
        }
        // Cap flips first, on the meshes that lean: a flip changes the
        // copies just collected, so a round that flips collects again.
        let mut flipped_any = false;
        for &mi in &leaning {
            let fm = &mut meshes[mi];
            if let Some(surface) = brep
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
                && flip_leaning_caps(fm, surface, bound)
            {
                changed.push(mi);
                flipped_any = true;
            }
        }
        if flipped_any {
            continue;
        }
        let mut touched: std::collections::HashSet<(usize, usize)> = Default::default();
        let topology: std::cell::RefCell<HashMap<usize, MeshTopology>> = Default::default();
        let mut split_any = false;
        for key in requests {
            let Some(list) = copies.get(&key) else {
                continue;
            };
            if list
                .iter()
                .any(|&(mi, ti, _, _)| touched.contains(&(mi, ti)))
            {
                continue;
            }
            let (m0, _, a0, b0) = list[0];
            let (pa, pb) = (
                meshes[m0].vertices[a0].position,
                meshes[m0].vertices[b0].position,
            );
            let faces: Vec<FaceHandle> = list.iter().map(|&(mi, ..)| meshes[mi].face).collect();
            // The apex of the copy leaning most: the split may go at its
            // foot when the chord's midpoint does not hold.
            let apex = list
                .iter()
                .filter_map(|&(mi, ti, a, b)| {
                    let fm = &meshes[mi];
                    let t = fm.triangles[ti];
                    let surface = brep.faces.get(fm.face).and_then(|f| f.surface.as_ref())?;
                    let tilt = triangle_tilt(fm, &surface.geometry, t)?;
                    let c = t.iter().copied().find(|&v| v != a && v != b)?;
                    Some((tilt, fm.vertices[c].position))
                })
                .max_by(|x, y| x.0.total_cmp(&y.0))
                .map(|(_, p)| p);
            let mut landed: Option<Vec<Copy>> = None;
            for toward in [None, apex] {
                let Some(m) =
                    edge_curve_midpoint(brep, &faces, pa, pb, gate, &mut curve_cache, toward)
                else {
                    continue;
                };
                // Plan every copy; land the split only if all of them hold.
                let mut plan: Vec<Copy> = Vec::new();
                let ok = list.iter().all(|&(mi, ti, a, b)| {
                    let fm = &meshes[mi];
                    let Some(surface) = brep
                        .faces
                        .get(fm.face)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry)
                    else {
                        return false;
                    };
                    let (va, vb) = (&fm.vertices[a], &fm.vertices[b]);
                    let hint = (0.5 * (va.uv.0 + vb.uv.0), 0.5 * (va.uv.1 + vb.uv.1));
                    let Some(uv) = popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        surface, &m, hint.0, hint.1,
                    ) else {
                        return false;
                    };
                    if surface.point_at(uv.0, uv.1).distance_to(&m) > gate
                        || singular_without_unique_normal(surface, uv.0, uv.1)
                    {
                        return false;
                    }
                    let Some(n) = limit_normal(surface, uv.0, uv.1) else {
                        return false;
                    };
                    let facing = va.normal + vb.normal;
                    let n = if n.dot(&facing) < 0.0 { -n } else { n };
                    let t = fm.triangles[ti];
                    let Some(r) = (0..3).find(|&k| t[k] == a && t[(k + 1) % 3] == b) else {
                        return false;
                    };
                    let c = t[(r + 2) % 3];
                    let vc = &fm.vertices[c];
                    let old = triangle_tilt(fm, surface, t).unwrap_or(f64::INFINITY);
                    let new = [
                        tilt_of([va.position, m, vc.position], [va.normal, n, vc.normal]),
                        tilt_of([m, vb.position, vc.position], [n, vb.normal, vc.normal]),
                    ];
                    let mut worst = new
                        .iter()
                        .map(|x| x.unwrap_or(f64::INFINITY))
                        .fold(0.0, f64::max);
                    let halves_spread = spread_of([va.normal, n, vc.normal])
                        .max(spread_of([n, vb.normal, vc.normal]));
                    if halves_spread > bound.max(triangle_spread(fm, t)) {
                        worst = f64::INFINITY; // the plain split may not land
                    }
                    let vertex = MeshVertex {
                        position: m,
                        normal: n,
                        uv,
                    };
                    // A chord whose sagitta tilts its cap may need two halvings:
                    // the first may leave the halves slightly worse (a cusp tip
                    // on a cylinder rim, 12.8 -> 13.7 deg), the next one cures
                    // them (sagitta falls with the square). A leaning cap's
                    // split may lean up to STEP_SLACK more; the rounds bound it.
                    // A copy that does not lean may cross the bound by as much:
                    // its face takes the tilt finish after this pass (the
                    // propeller's hub copy of a blade-edge chord: 5.5 -> 12.5).
                    let allowance = old.max(bound) + STEP_SLACK;
                    let collapse = if worst > allowance {
                        // The apex sits in the chord's cap zone: both halves
                        // keep leaning toward it. Remove it instead.
                        let mut topo = topology.borrow_mut();
                        let topo = topo.entry(mi).or_insert_with(|| MeshTopology::of(fm));
                        Some(plan_apex_collapse(fm, surface, ti, c, &vertex, topo, bound))
                    } else {
                        None
                    };
                    if collapse.as_ref().is_some_and(Option::is_none) {
                        return false;
                    }
                    plan.push(Copy {
                        mesh: mi,
                        triangle: ti,
                        a,
                        b,
                        vertex,
                        collapse: collapse.flatten(),
                    });
                    true
                });
                if ok && !plan.is_empty() {
                    landed = Some(plan);
                    break;
                }
            }
            let Some(plan) = landed else {
                continue;
            };
            for cp in plan {
                let fm = &mut meshes[cp.mesh];
                let t = fm.triangles[cp.triangle];
                let Some(r) = (0..3).find(|&k| t[k] == cp.a && t[(k + 1) % 3] == cp.b) else {
                    continue;
                };
                let c = t[(r + 2) % 3];
                let mv = fm.vertices.len();
                fm.vertices.push(cp.vertex);
                let mut topo = topology.borrow_mut();
                let topo = topo.get_mut(&cp.mesh);
                match &cp.collapse {
                    None => {
                        let n = fm.triangles.len();
                        fm.triangles[cp.triangle] = [cp.a, mv, c];
                        fm.triangles.push([mv, cp.b, c]);
                        if let Some(topo) = topo {
                            topo.unlink(cp.b, cp.triangle);
                            topo.incident.entry(cp.b).or_default().push(n);
                            topo.incident.entry(c).or_default().push(n);
                            topo.incident.insert(mv, vec![cp.triangle, n]);
                            topo.boundary.insert(mv);
                        }
                    }
                    Some(star) => {
                        fm.triangles[cp.triangle] = [mv, mv, mv]; // purged below
                        for &(ti, tri) in star {
                            fm.triangles[ti] = tri.map(|i| if i == NEW { mv } else { i });
                            touched.insert((cp.mesh, ti));
                        }
                        if let Some(topo) = topo {
                            for v in [cp.a, cp.b] {
                                topo.unlink(v, cp.triangle);
                            }
                            topo.incident.remove(&c);
                            topo.incident
                                .insert(mv, star.iter().map(|&(ti, _)| ti).collect());
                            topo.boundary.insert(mv);
                        }
                    }
                }
                touched.insert((cp.mesh, cp.triangle));
                changed.push(cp.mesh);
                split_any = true;
            }
        }
        for fm in meshes.iter_mut() {
            fm.triangles
                .retain(|t| t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
        }
        if !split_any {
            break;
        }
    }
    changed.sort_unstable();
    changed.dedup();
    changed
}
