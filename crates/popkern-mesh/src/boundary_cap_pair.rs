//! Paired chord split for leaning caps on a narrow trimmed strip
//! (PK-859; E1 of docs/architecture/display-tessellation-in-face-folds.md).
//!
//! A face only one triangle wide between two trims (a thread-root
//! cylinder strip ~0.0126 across) meshes as a row of caps, each with
//! its chord on one trim and its apex on the other. Such a cap leans by
//! about L^2 / (8 r w) for chord length L across width w: halving its
//! chord alone keeps the apex across the strip (the halves lean MORE),
//! and its apex cannot be collapsed (it is a boundary vertex). Halving
//! the chord AND the opposite trim's chord of the neighbouring cap in
//! one transaction quarters it: the quad the two caps make is
//! re-triangulated through both new vertices. Every copy of either chord
//! in another face mesh is split at the same point, so the weld holds by
//! construction; the pair lands only when no new facet in any copy
//! leans more than the facets it replaces (or the bound).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

use super::curve::{CurveSamples, edge_curve_midpoint};
use super::{
    EdgeKey, STEP_SLACK, boundary_edges, chart_area, edge_key, spread_of, tilt_of, triangle_spread,
    triangle_tilt,
};
use crate::mesh::{FaceMesh, MeshVertex};
use crate::vertex_normal::{limit_normal, singular_without_unique_normal};

#[path = "boundary_cap_corner.rs"]
mod corner;
use corner::plan_chord_with_corners;

/// The face mesh's surface.
fn surface_of<'a>(brep: &'a BRep, fm: &FaceMesh) -> Option<&'a SurfaceGeometry> {
    brep.faces
        .get(fm.face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
}

/// The mesh vertex for edge point `m` between local vertices `a` and
/// `b` of `fm`: its chart position near their midpoint and the surface
/// normal oriented with theirs. `None` off the surface or at a chart
/// singularity.
fn vertex_at(
    fm: &FaceMesh,
    surface: &SurfaceGeometry,
    a: usize,
    b: usize,
    m: Point3,
    gate: f64,
) -> Option<MeshVertex> {
    let (va, vb) = (&fm.vertices[a], &fm.vertices[b]);
    let hint = (0.5 * (va.uv.0 + vb.uv.0), 0.5 * (va.uv.1 + vb.uv.1));
    let raw = popkern_core::geometry::uv_projection::project_point_to_uv_near(
        surface, &m, hint.0, hint.1,
    )?;
    // The chord's own branch of a periodic chart: a projection lands in
    // the principal range, and a helical strip's chart is unwrapped
    // over many turns (the new vertex sat whole turns from its chord).
    let (pu, pv) = surface.periodicity();
    let align = |x: f64, near: f64, period: Option<f64>| match period {
        Some(p) if p > 0.0 => x - ((x - near) / p).round() * p,
        _ => x,
    };
    let uv = (align(raw.0, hint.0, pu), align(raw.1, hint.1, pv));
    if surface.point_at(uv.0, uv.1).distance_to(&m) > gate
        || singular_without_unique_normal(surface, uv.0, uv.1)
    {
        return None;
    }
    let n = limit_normal(surface, uv.0, uv.1)?;
    let facing = va.normal + vb.normal;
    let normal = if n.dot(&facing) < 0.0 { -n } else { n };
    Some(MeshVertex {
        position: m,
        normal,
        uv,
    })
}

/// A corner of a planned triangle: an existing vertex or a new one.
#[derive(Clone, Copy)]
enum Corner {
    Old(usize),
    New(usize),
}

/// One mesh's part of a planned pair split: the triangles it replaces,
/// the new vertices and the triangles that take their place.
struct MeshPlan {
    mesh: usize,
    remove: Vec<usize>,
    vertices: Vec<MeshVertex>,
    triangles: Vec<[Corner; 3]>,
}

/// Position, normal and chart point of corner `c` under `plan`.
fn corner(
    fm: &FaceMesh,
    vertices: &[MeshVertex],
    c: Corner,
) -> (Point3, popkern_core::geometry::Vec3, (f64, f64)) {
    let v = match c {
        Corner::Old(i) => &fm.vertices[i],
        Corner::New(i) => &vertices[i],
    };
    (v.position, v.normal, v.uv)
}

/// Worst tilt of the planned triangles, or `None` when one loses its
/// chart winding (`sense`), its area, faces against its normals, or
/// widens the spread past `spread_cap`.
fn plan_tilt(
    fm: &FaceMesh,
    vertices: &[MeshVertex],
    tris: &[[Corner; 3]],
    sense: f64,
    spread_cap: f64,
) -> Option<f64> {
    let mut worst = 0.0f64;
    for t in tris {
        let q = t.map(|c| corner(fm, vertices, c));
        let p = q.map(|x| x.2);
        let area =
            0.5 * ((p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1));
        if area.signum() != sense || area.abs() <= 1e-15 {
            return None;
        }
        let f = (q[1].0 - q[0].0).cross(&(q[2].0 - q[0].0));
        if q.iter().any(|x| f.dot(&x.1) <= 0.0) || spread_of(q.map(|x| x.1)) > spread_cap {
            return None;
        }
        worst = worst.max(tilt_of(q.map(|x| x.0), q.map(|x| x.1))?);
    }
    Some(worst)
}

/// The plain split of another mesh's copy `(mi, ti, a, b)` of a chord at
/// `m`: its triangle halves, each leaning no more than the old one (or
/// the bound) plus the step slack.
fn plan_copy(
    brep: &BRep,
    meshes: &[FaceMesh],
    copy: (usize, usize, usize, usize),
    m: Point3,
    gate: f64,
    bound: f64,
) -> Option<MeshPlan> {
    let (mi, ti, a, b) = copy;
    let fm = &meshes[mi];
    let surface = surface_of(brep, fm)?;
    let t = fm.triangles[ti];
    let r = (0..3).find(|&k| t[k] == a && t[(k + 1) % 3] == b)?;
    let c = t[(r + 2) % 3];
    let vertex = vertex_at(fm, surface, a, b, m, gate)?;
    let tris = vec![
        [Corner::Old(a), Corner::New(0), Corner::Old(c)],
        [Corner::New(0), Corner::Old(b), Corner::Old(c)],
    ];
    let old = triangle_tilt(fm, surface, t).unwrap_or(f64::INFINITY);
    let sense = chart_area(fm, t).signum();
    let worst = plan_tilt(
        fm,
        std::slice::from_ref(&vertex),
        &tris,
        sense,
        bound.max(triangle_spread(fm, t)),
    )?;
    (worst <= old.max(bound) + STEP_SLACK).then(|| MeshPlan {
        mesh: mi,
        remove: vec![ti],
        vertices: vec![vertex],
        triangles: tris,
    })
}

/// Edge -> triangles of one mesh.
type EdgeMap = HashMap<(usize, usize), Vec<usize>>;

/// Each face mesh's edge map, built once per round on first use and
/// dropped when a plan lands on that mesh (the coil's 36k-triangle tube
/// rebuilt it per cap: +8% display time).
type EdgeMaps = HashMap<usize, EdgeMap>;

fn edge_triangles(fm: &FaceMesh) -> EdgeMap {
    let mut map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ti, t) in fm.triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            map.entry((a.min(b), a.max(b))).or_default().push(ti);
        }
    }
    map
}

/// Plan the pair split of leaning cap `t1` of mesh `mi`, whose chord is
/// boundary edge `p0 -> p1`. The neighbour across one of the cap's
/// apex edges must itself be a cap on the opposite trim; the quad
/// `p0 p1 p2 p3` (chords `p0 p1` and `p2 p3`, in winding order) becomes
/// four triangles through both chords' midpoints, the better of the
/// diagonal choices on each side.
#[allow(clippy::too_many_arguments)]
fn plan_pair(
    brep: &BRep,
    meshes: &[FaceMesh],
    copies: &HashMap<EdgeKey, Vec<(usize, usize, usize, usize)>>,
    mi: usize,
    t1: usize,
    (p0, p1): (usize, usize),
    gate: f64,
    bound: f64,
    cache: &mut CurveSamples,
    maps: &mut EdgeMaps,
) -> Option<Vec<MeshPlan>> {
    let fm = &meshes[mi];
    let surface = surface_of(brep, fm)?;
    let tri1 = fm.triangles[t1];
    let r = (0..3).find(|&k| tri1[k] == p0 && tri1[(k + 1) % 3] == p1)?;
    let c = tri1[(r + 2) % 3];
    let edges: &EdgeMap = maps.entry(mi).or_insert_with(|| edge_triangles(fm));
    let boundary = |x: usize, y: usize| {
        edges
            .get(&(x.min(y), x.max(y)))
            .is_some_and(|v| v.len() == 1)
    };
    let old1 = triangle_tilt(fm, surface, tri1)?;
    let sense = chart_area(fm, tri1).signum();
    let mut best: Option<(f64, Vec<MeshPlan>)> = None;
    // The neighbour across p1-c (quad p0 p1 d c) or across c-p0 (quad
    // p0 p1 c e): its chord on the other trim is d -> c or c -> e.
    for across in [(p1, c), (c, p0)] {
        let Some(&[x, y]) = edges
            .get(&(across.0.min(across.1), across.0.max(across.1)))
            .map(Vec::as_slice)
        else {
            continue;
        };
        let t2 = if x == t1 { y } else { x };
        let tri2 = fm.triangles[t2];
        let Some(o) = tri2
            .iter()
            .copied()
            .find(|&v| v != across.0 && v != across.1)
        else {
            continue;
        };
        let quad = if across == (p1, c) {
            [p0, p1, o, c]
        } else {
            [p0, p1, c, o]
        };
        if !boundary(quad[2], quad[3]) {
            continue;
        }
        let Some(old2) = triangle_tilt(fm, surface, tri2) else {
            continue;
        };
        let old = old1.max(old2);
        // Both chords' edge-curve midpoints, from every face that carries them.
        let chord = |a: usize, b: usize| edge_key(fm.vertices[a].position, fm.vertices[b].position);
        let (ka, kb) = (chord(quad[0], quad[1]), chord(quad[2], quad[3]));
        let (Some(list_a), Some(list_b)) = (copies.get(&ka), copies.get(&kb)) else {
            continue;
        };
        let faces = |list: &[(usize, usize, usize, usize)]| {
            list.iter()
                .map(|&(m, ..)| meshes[m].face)
                .collect::<Vec<FaceHandle>>()
        };
        let pos = |i: usize| fm.vertices[i].position;
        let Some(ma) = edge_curve_midpoint(
            brep,
            &faces(list_a),
            pos(quad[0]),
            pos(quad[1]),
            gate,
            cache,
            None,
        ) else {
            continue;
        };
        let Some(mb) = edge_curve_midpoint(
            brep,
            &faces(list_b),
            pos(quad[2]),
            pos(quad[3]),
            gate,
            cache,
            None,
        ) else {
            continue;
        };
        let (Some(va), Some(vb)) = (
            vertex_at(fm, surface, quad[0], quad[1], ma, gate),
            vertex_at(fm, surface, quad[2], quad[3], mb, gate),
        ) else {
            continue;
        };
        let vertices = vec![va, vb];
        let (q0, q1, q2, q3) = (
            Corner::Old(quad[0]),
            Corner::Old(quad[1]),
            Corner::Old(quad[2]),
            Corner::Old(quad[3]),
        );
        let (m, n) = (Corner::New(0), Corner::New(1));
        let spread_cap = bound
            .max(triangle_spread(fm, tri1))
            .max(triangle_spread(fm, tri2));
        // Hexagon q0 m q1 q2 n q3: left quad q0 m n q3, right quad m q1 q2 n.
        let mut here: Option<(f64, Vec<[Corner; 3]>)> = None;
        for left in [[[q0, m, n], [q0, n, q3]], [[q0, m, q3], [m, n, q3]]] {
            for right in [[[m, q1, q2], [m, q2, n]], [[m, q1, n], [q1, q2, n]]] {
                let tris: Vec<[Corner; 3]> = left.iter().chain(right.iter()).copied().collect();
                if let Some(w) = plan_tilt(fm, &vertices, &tris, sense, spread_cap)
                    && here.as_ref().is_none_or(|h| w < h.0)
                {
                    here = Some((w, tris));
                }
            }
        }
        let Some((worst, tris)) = here else {
            continue;
        };
        if worst >= old || best.as_ref().is_some_and(|b| worst >= b.0) {
            continue;
        }
        // Every other copy of either chord takes the plain split.
        let mut plans = vec![MeshPlan {
            mesh: mi,
            remove: vec![t1, t2],
            vertices,
            triangles: tris,
        }];
        let others = list_a
            .iter()
            .map(|&cp| (cp, ma))
            .chain(list_b.iter().map(|&cp| (cp, mb)))
            .filter(|&((m2, t, ..), _)| !(m2 == mi && (t == t1 || t == t2)));
        let mut ok = true;
        for (cp, point) in others {
            match plan_copy(brep, meshes, cp, point, gate, bound) {
                Some(p) => plans.push(p),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            best = Some((worst, plans));
        }
    }
    best.map(|(_, plans)| plans)
}

/// One round of pair (and corner) splits over the face meshes in
/// `leaning`: every leaning boundary cap whose plan holds lands (one per
/// triangle per round). Returns the face meshes changed.
pub(super) fn split_leaning_pairs(
    brep: &BRep,
    meshes: &mut [FaceMesh],
    leaning: &[usize],
    copies: &HashMap<EdgeKey, Vec<(usize, usize, usize, usize)>>,
    bound: f64,
    gate: f64,
    cache: &mut CurveSamples,
) -> Vec<usize> {
    let mut caps: Vec<(f64, usize, usize, (usize, usize))> = Vec::new();
    // Only the meshes the round found leaning boundary caps on (scanning
    // every mesh's boundary per round cost the gallery ~5% display time).
    for &mi in leaning {
        let fm = &meshes[mi];
        let Some(surface) = surface_of(brep, fm) else {
            continue;
        };
        for (a, b, ti) in boundary_edges(fm) {
            if let Some(t) = triangle_tilt(fm, surface, fm.triangles[ti])
                && t > bound
            {
                caps.push((t, mi, ti, (a, b)));
            }
        }
    }
    caps.sort_by(|x, y| y.0.total_cmp(&x.0));
    let mut touched: std::collections::HashSet<(usize, usize)> = Default::default();
    let mut changed = Vec::new();
    let mut maps: EdgeMaps = HashMap::new();
    for (_, mi, ti, chord) in caps {
        if touched.contains(&(mi, ti)) {
            continue;
        }
        let Some(plans) = plan_pair(
            brep, meshes, copies, mi, ti, chord, gate, bound, cache, &mut maps,
        )
        .or_else(|| {
            plan_chord_with_corners(
                brep,
                meshes,
                copies,
                (mi, ti),
                chord,
                gate,
                bound,
                cache,
                &mut maps,
            )
        }) else {
            continue;
        };
        if plans
            .iter()
            .any(|p| p.remove.iter().any(|&t| touched.contains(&(p.mesh, t))))
        {
            continue;
        }
        // One transaction rewrites each triangle once (a neighbour a
        // corner plan takes may be another copy's triangle).
        let mut claimed = std::collections::HashSet::new();
        if !plans
            .iter()
            .all(|p| p.remove.iter().all(|&t| claimed.insert((p.mesh, t))))
        {
            continue;
        }
        for plan in plans {
            maps.remove(&plan.mesh);
            let fm = &mut meshes[plan.mesh];
            let base = fm.vertices.len();
            fm.vertices.extend(plan.vertices);
            let map = |c: Corner| match c {
                Corner::Old(i) => i,
                Corner::New(i) => base + i,
            };
            let mut slots = plan.remove.clone();
            for t in &plan.triangles {
                let tri = t.map(map);
                match slots.pop() {
                    Some(s) => fm.triangles[s] = tri,
                    None => {
                        fm.triangles.push(tri);
                        let n = fm.triangles.len() - 1;
                        touched.insert((plan.mesh, n));
                    }
                }
            }
            for &t in &plan.remove {
                touched.insert((plan.mesh, t));
            }
            changed.push(plan.mesh);
        }
    }
    changed
}
