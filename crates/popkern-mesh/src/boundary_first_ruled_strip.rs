//! Structured tensor tessellation for four-sided ruled charts — the
//! fillet/chamfer strip class (PK-218).
//!
//! A blend strip's chart is a UV rectangle: two rails at constant u
//! and two cross chains (end arcs) at constant v, with the surface
//! straight along v (ruled). Delaunay-with-protection is the wrong
//! tool for these charts: the protecting disks of the long rail
//! spans evacuate the interior lattice (the chart is centimetres in
//! v but only a couple of millimetres of metric width in u), and the
//! boundary-only triangulation that survives chords rail-to-rail
//! across the arc (measured: plate-rim strips read 0.5% low in mesh
//! volume while every vertex lies on the cylinder).
//!
//! The structured mesh is exact by construction instead: boundary
//! vertices are the canonical ring VERBATIM (watertight with the
//! neighbours), interior vertices sit on the cross-station × rail-
//! station lattice evaluated on the surface, and the chord error is
//! the adaptive edge sampling's own chord tolerance. No protection
//! filters, no encroachment rounds.
//!
//! Any decomposition failure declines (with a `PK_BF_TRACE` reason)
//! and the caller continues on the Delaunay path.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::FaceHandle;

use super::BoundaryRing;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

/// One boundary chain of the rectangle decomposition: ring indices
/// (into the single ring's point list) ordered along the chain,
/// inclusive of both corner endpoints.
struct Chain {
    ring_indices: Vec<usize>,
    /// True when the chain holds u constant (a rail); false when it
    /// holds v constant (a cross chain).
    u_const: bool,
}

fn trace(reason: &str) {
    if std::env::var_os("PK_BF_TRACE").is_some() {
        eprintln!("    [bf] ruled-strip declined: {reason}");
    }
}

/// Attempt the structured tensor tessellation. Returns `None` (after
/// recording why under `PK_BF_TRACE`) whenever the face is not a
/// clean four-sided ruled chart; the caller then proceeds with the
/// Delaunay path unchanged.
pub(super) fn try_tessellate(
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    orientation: bool,
    rings: &[BoundaryRing],
    options: &TessellationOptions,
) -> Option<FaceMesh> {
    if std::env::var("PK_BF_RULED_STRIP").as_deref() == Ok("0") {
        return None;
    }
    if matches!(surface, SurfaceGeometry::Plane { .. }) {
        // Planar charts have zero chord error on the Delaunay path
        // already; keep their behaviour untouched.
        return None;
    }
    let [ring] = rings else {
        trace("inner loops present");
        return None;
    };
    let n = ring.uv.len();
    if n < 4 {
        return None;
    }

    // Chart bbox and metric scale (mm per parameter unit in each
    // direction, sampled at the chart centre) — corner detection and
    // constancy tests must run in metric space, not raw UV, because
    // the chart mixes units (radians × millimetres on a cylinder).
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &ring.uv {
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    let (su, sv) = {
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let p = surface.point_at(um, vm);
        let eps_u = ((u1 - u0).abs() * 1e-4).max(1e-9);
        let eps_v = ((v1 - v0).abs() * 1e-4).max(1e-9);
        let du = surface.point_at(um + eps_u, vm).distance_to(&p) / eps_u;
        let dv = surface.point_at(um, vm + eps_v).distance_to(&p) / eps_v;
        (du.max(1e-9), dv.max(1e-9))
    };
    let diag = (su * (u1 - u0)).hypot(sv * (v1 - v0)).max(1e-12);

    // Ruled along v: the v-midline must sit on the chord within the
    // chord tolerance (same test as the caller's ruled path).
    {
        let u_mid = 0.5 * (u0 + u1);
        let pa = surface.point_at(u_mid, v0);
        let pb = surface.point_at(u_mid, v1);
        let pm = surface.point_at(u_mid, 0.5 * (v0 + v1));
        let chord_mid = Point3::new(
            0.5 * (pa.x + pb.x),
            0.5 * (pa.y + pb.y),
            0.5 * (pa.z + pb.z),
        );
        if pm.distance_to(&chord_mid) > options.chord_height {
            return None; // doubly curved — not this path's class
        }
    }

    // Corner detection on the metric polyline: a corner is a turn
    // sharper than 45°. A clean rectangle has exactly four.
    let m: Vec<(f64, f64)> = ring.uv.iter().map(|&(u, v)| (u * su, v * sv)).collect();
    let mut corners: Vec<usize> = Vec::new();
    for i in 0..n {
        let a = m[(i + n - 1) % n];
        let b = m[i];
        let c = m[(i + 1) % n];
        let (d1x, d1y) = (b.0 - a.0, b.1 - a.1);
        let (d2x, d2y) = (c.0 - b.0, c.1 - b.1);
        let l1 = (d1x * d1x + d1y * d1y).sqrt();
        let l2 = (d2x * d2x + d2y * d2y).sqrt();
        if l1 < 1e-300 || l2 < 1e-300 {
            trace("zero-length metric segment");
            return None;
        }
        let cos = (d1x * d2x + d1y * d2y) / (l1 * l2);
        if cos < std::f64::consts::FRAC_1_SQRT_2 {
            corners.push(i);
        }
    }
    if corners.len() != 4 {
        trace(&format!("{} corners, need 4", corners.len()));
        return None;
    }

    // Split the ring into four chains between consecutive corners
    // (each chain includes both corner endpoints) and classify each
    // as u-constant (rail) or v-constant (cross chain).
    let eps_const = 1e-6 * diag;
    let mut chains: Vec<Chain> = Vec::with_capacity(4);
    for k in 0..4 {
        let start = corners[k];
        let end = corners[(k + 1) % 4];
        let mut ring_indices = vec![start];
        let mut i = start;
        while i != end {
            i = (i + 1) % n;
            ring_indices.push(i);
        }
        let xs: Vec<f64> = ring_indices.iter().map(|&i| m[i].0).collect();
        let ys: Vec<f64> = ring_indices.iter().map(|&i| m[i].1).collect();
        let x_span = xs.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
            - xs.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let y_span = ys.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
            - ys.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let u_const = x_span <= eps_const;
        let v_const = y_span <= eps_const;
        if u_const == v_const {
            trace("chain neither u-constant nor v-constant");
            return None;
        }
        chains.push(Chain {
            ring_indices,
            u_const,
        });
    }
    if !(chains[0].u_const != chains[1].u_const && chains[1].u_const != chains[2].u_const) {
        trace("chains do not alternate rail/cross");
        return None;
    }

    // Sort the two rails by u and the two cross chains by v; orient
    // every chain ascending in its running coordinate.
    let mut rails: Vec<Vec<usize>> = Vec::new();
    let mut crosses: Vec<Vec<usize>> = Vec::new();
    for chain in &chains {
        let mut idx = chain.ring_indices.clone();
        let coord = |i: usize| {
            if chain.u_const {
                ring.uv[i].1
            } else {
                ring.uv[i].0
            }
        };
        let ascending = coord(*idx.last()?) >= coord(idx[0]);
        if !ascending {
            idx.reverse();
        }
        // The running coordinate must be strictly monotone or the
        // chain cannot index a tensor row/column.
        for w in idx.windows(2) {
            if coord(w[1]) <= coord(w[0]) {
                trace("chain not strictly monotone");
                return None;
            }
        }
        if chain.u_const {
            rails.push(idx);
        } else {
            crosses.push(idx);
        }
    }
    rails.sort_by(|a, b| ring.uv[a[0]].0.total_cmp(&ring.uv[b[0]].0));
    crosses.sort_by(|a, b| ring.uv[a[0]].1.total_cmp(&ring.uv[b[0]].1));
    let (rail_l, rail_r) = (&rails[0], &rails[1]);
    let (cross_a, cross_b) = (&crosses[0], &crosses[1]);
    if su * (ring.uv[rail_r[0]].0 - ring.uv[rail_l[0]].0) <= eps_const
        || sv * (ring.uv[cross_b[0]].1 - ring.uv[cross_a[0]].1) <= eps_const
    {
        trace("degenerate rectangle extent");
        return None;
    }

    // Station matching: the tensor grid pairs cross-chain station i
    // on both ends and rail station j on both sides. The chains must
    // agree station-for-station within half the smallest gap, or the
    // columns/rows would shear.
    let nu = cross_a.len();
    let nv = rail_l.len();
    if cross_b.len() != nu || rail_r.len() != nv {
        trace("station counts differ between opposite chains");
        return None;
    }
    let min_gap = |idx: &[usize], pick_u: bool, scale: f64| -> f64 {
        idx.windows(2)
            .map(|w| {
                let a = if pick_u {
                    ring.uv[w[0]].0
                } else {
                    ring.uv[w[0]].1
                };
                let b = if pick_u {
                    ring.uv[w[1]].0
                } else {
                    ring.uv[w[1]].1
                };
                scale * (b - a)
            })
            .fold(f64::INFINITY, f64::min)
    };
    let gap_u = min_gap(cross_a, true, su).min(min_gap(cross_b, true, su));
    let gap_v = min_gap(rail_l, false, sv).min(min_gap(rail_r, false, sv));
    for i in 0..nu {
        if su * (ring.uv[cross_a[i]].0 - ring.uv[cross_b[i]].0).abs() > 0.5 * gap_u {
            trace("cross-chain stations misaligned");
            return None;
        }
    }
    for j in 0..nv {
        if sv * (ring.uv[rail_l[j]].1 - ring.uv[rail_r[j]].1).abs() > 0.5 * gap_v {
            trace("rail stations misaligned");
            return None;
        }
    }

    // Assemble: ring points verbatim first (canonical positions →
    // watertight with the neighbours), then the interior lattice.
    let mut uv_pts: Vec<(f64, f64)> = ring.uv.clone();
    let mut pos_pts: Vec<Option<Point3>> = ring.pos.iter().copied().map(Some).collect();
    let interior_base = uv_pts.len();
    for i in 1..nu.saturating_sub(1) {
        let u = 0.5 * (ring.uv[cross_a[i]].0 + ring.uv[cross_b[i]].0);
        for j in 1..nv.saturating_sub(1) {
            let v = 0.5 * (ring.uv[rail_l[j]].1 + ring.uv[rail_r[j]].1);
            uv_pts.push((u, v));
            pos_pts.push(None);
        }
    }
    let vertex_at = |i: usize, j: usize| -> usize {
        if j == 0 {
            cross_a[i]
        } else if j == nv - 1 {
            cross_b[i]
        } else if i == 0 {
            rail_l[j]
        } else if i == nu - 1 {
            rail_r[j]
        } else {
            interior_base + (i - 1) * (nv - 2) + (j - 1)
        }
    };
    let mut triangles: Vec<[usize; 3]> = Vec::with_capacity(2 * (nu - 1) * (nv - 1));
    for i in 0..nu - 1 {
        for j in 0..nv - 1 {
            let (p00, p10, p11, p01) = (
                vertex_at(i, j),
                vertex_at(i + 1, j),
                vertex_at(i + 1, j + 1),
                vertex_at(i, j + 1),
            );
            triangles.push([p00, p10, p11]);
            triangles.push([p00, p11, p01]);
        }
    }
    // Enforce positive UV winding (the quads are built ascending in
    // both axes, but the surface chart may run either handedness).
    for t in &mut triangles {
        let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        let signed = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if signed < 0.0 {
            t.swap(1, 2);
        }
    }

    // Materialize exactly like the Delaunay tail: canonical boundary
    // positions where recorded, surface evaluation for the interior,
    // normals from the chart with the stored orientation flag (the
    // global winding-coherence pass in `tessellate` remains the
    // final authority).
    let mut vertices = Vec::with_capacity(uv_pts.len());
    for (k, &(u, v)) in uv_pts.iter().enumerate() {
        let position = pos_pts[k].unwrap_or_else(|| surface.point_at(u, v));
        let mut normal = surface.normal_at(u, v);
        let len = normal.length();
        if len > 1e-300 {
            normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
        }
        if !orientation {
            normal = -normal;
        }
        vertices.push(MeshVertex {
            position,
            normal,
            uv: (u, v),
        });
    }
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    if std::env::var_os("PK_BF_RING_DUMP").is_some() {
        eprintln!(
            "    [bf-ruled-strip] face {face_handle:?}: {nu}x{nv} stations, {} tris",
            triangles.len()
        );
    }
    Some(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}
