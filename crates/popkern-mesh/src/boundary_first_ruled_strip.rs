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
//! Any decomposition failure declines (with a the `boundary-first` trace topic reason)
//! and the caller continues on the Delaunay path.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, PortableMath, Vec3};
use popkern_core::topology::FaceHandle;

use super::BoundaryRing;
#[path = "boundary_first_ruled_strip_interior.rs"]
mod interior;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};
use interior::resolves_interior;

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
    if popkern_core::trace::trace_enabled("boundary-first") {
        eprintln!("    [bf] ruled-strip declined: {reason}");
    }
}

/// Attempt the structured tensor tessellation. Returns `None` (after
/// recording why under the `boundary-first` trace topic) whenever the face is not a
/// clean four-sided ruled chart; the caller then proceeds with the
/// Delaunay path unchanged.
pub(super) fn try_tessellate(
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    orientation: bool,
    rings: &[BoundaryRing],
    options: &TessellationOptions,
) -> Option<FaceMesh> {
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

    // Ruling axis: the chart direction along which the surface is
    // straight. Strips ruled along v (cylinder-like charts) keep the
    // chart as-is; strips ruled along u (a blend whose chart runs u
    // along the spine and v across the section — variable-radius
    // rounds) are TRANSPOSED so the rest of this builder always sees
    // the ruling as its second coordinate. Before the transpose,
    // u-ruled strips declined into Delaunay, whose isotropic clearance
    // band (sized by the strip's long direction) evacuated the whole
    // interior: a tapered round meshed as 12 rail-to-rail triangles
    // with 90° of normal turn each and 7x the chord tolerance.
    let bbox = |pts: &[(f64, f64)]| {
        pts.iter().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(a0, a1, b0, b1), &(a, b)| (a0.min(a), a1.max(a), b0.min(b), b1.max(b)),
        )
    };
    let straight_along = |swap: bool| -> bool {
        let (u0, u1, v0, v1) = bbox(&ring.uv);
        let at = |a: f64, b: f64| {
            if swap {
                surface.point_at(b, a)
            } else {
                surface.point_at(a, b)
            }
        };
        let (a0, a1, b0, b1) = if swap {
            (v0, v1, u0, u1)
        } else {
            (u0, u1, v0, v1)
        };
        let a_mid = 0.5 * (a0 + a1);
        let pa = at(a_mid, b0);
        let pb = at(a_mid, b1);
        let pm = at(a_mid, 0.5 * (b0 + b1));
        let chord_mid = Point3::new(
            0.5 * (pa.x + pb.x),
            0.5 * (pa.y + pb.y),
            0.5 * (pa.z + pb.z),
        );
        pm.distance_to(&chord_mid) <= options.chord_height
    };
    // A DOUBLY curved four-sided chart (a slit torus fillet band)
    // takes the same tensor lattice: protecting disks evacuated a
    // 662-long 1.26-wide band (facets 0.21 off a 0.8 tube) and
    // isotropic refinement cost 5x the triangles (PK-650).
    let doubly_curved = !straight_along(false) && !straight_along(true);
    // Analytic only: NURBS thread flanks' folded seam charts tore
    // (79-202 open edges on fasteners; PK-651 owns NURBS strips).
    if doubly_curved
        && !matches!(
            surface,
            SurfaceGeometry::Torus { .. } | SurfaceGeometry::Sphere { .. }
        )
    {
        return None;
    }
    let swap = !doubly_curved && !straight_along(false);
    // Chart coordinates (a, b): b runs along the ruling.
    let chart: Vec<(f64, f64)> = ring
        .uv
        .iter()
        .map(|&(u, v)| if swap { (v, u) } else { (u, v) })
        .collect();
    let eval = |a: f64, b: f64| {
        if swap {
            surface.point_at(b, a)
        } else {
            surface.point_at(a, b)
        }
    };

    // Chart bbox and metric scale (mm per parameter unit in each
    // direction, sampled at the chart centre) — corner detection and
    // constancy tests must run in metric space, not raw UV, because
    // the chart mixes units (radians × millimetres on a cylinder).
    let (u0, u1, v0, v1) = bbox(&chart);
    let (su, sv) = {
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let p = eval(um, vm);
        let eps_u = ((u1 - u0).abs() * 1e-4).max(1e-9);
        let eps_v = ((v1 - v0).abs() * 1e-4).max(1e-9);
        let du = eval(um + eps_u, vm).distance_to(&p) / eps_u;
        let dv = eval(um, vm + eps_v).distance_to(&p) / eps_v;
        (du.max(1e-9), dv.max(1e-9))
    };
    let diag = (su * (u1 - u0)).phypot(sv * (v1 - v0)).max(1e-12);
    // Doubly curved: anisotropic only (a near-square wavy face meshed
    // 1.03% short on tensor diagonals vs 0.73% on Delaunay; PK-650).
    let (ext_u, ext_v) = (su * (u1 - u0), sv * (v1 - v0));
    if doubly_curved && ext_u.max(ext_v) < 8.0 * ext_u.min(ext_v) {
        trace("doubly curved chart is not anisotropic");
        return None;
    }

    // Corner detection on the metric polyline: a corner is a turn
    // sharper than 45°. A clean rectangle has exactly four.
    let m: Vec<(f64, f64)> = chart.iter().map(|&(u, v)| (u * su, v * sv)).collect();
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
                chart[i].1
            } else {
                chart[i].0
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
    rails.sort_by(|a, b| chart[a[0]].0.total_cmp(&chart[b[0]].0));
    crosses.sort_by(|a, b| chart[a[0]].1.total_cmp(&chart[b[0]].1));
    let (rail_l, rail_r) = (&rails[0], &rails[1]);
    let (cross_a, cross_b) = (&crosses[0], &crosses[1]);
    if su * (chart[rail_r[0]].0 - chart[rail_l[0]].0) <= eps_const
        || sv * (chart[cross_b[0]].1 - chart[cross_a[0]].1) <= eps_const
    {
        trace("degenerate rectangle extent");
        return None;
    }

    // Station matching. Rails pair station-for-station (the rows of
    // the lattice); cross chains pair station-for-station when they
    // align (the exact tensor lattice), and otherwise ZIP against a
    // union-station lattice — a variable-radius round samples its two
    // end sections at different densities (the chord tolerance meets
    // a 0.18 arc with more samples than a 0.04 one), which declined
    // the whole strip into Delaunay before.
    let nv = rail_l.len();
    if rail_r.len() != nv {
        trace("station counts differ between the rails");
        return None;
    }
    let min_gap = |idx: &[usize], pick_u: bool, scale: f64| -> f64 {
        idx.windows(2)
            .map(|w| {
                let (a, b) = if pick_u {
                    (chart[w[0]].0, chart[w[1]].0)
                } else {
                    (chart[w[0]].1, chart[w[1]].1)
                };
                scale * (b - a)
            })
            .fold(f64::INFINITY, f64::min)
    };
    let gap_u = min_gap(cross_a, true, su).min(min_gap(cross_b, true, su));
    let gap_v = min_gap(rail_l, false, sv).min(min_gap(rail_r, false, sv));
    for j in 0..nv {
        if sv * (chart[rail_l[j]].1 - chart[rail_r[j]].1).abs() > 0.5 * gap_v {
            trace("rail stations misaligned");
            return None;
        }
    }
    let aligned = cross_a.len() == cross_b.len()
        && (0..cross_a.len())
            .all(|i| su * (chart[cross_a[i]].0 - chart[cross_b[i]].0).abs() <= 0.5 * gap_u);
    if !aligned && !swap && !doubly_curved {
        // Chart-ruled-along-v strips with misaligned end chains keep
        // the Delaunay path: its encroachment rounds equalize the
        // SHARED rim samples (a drafted cone's two circles), which the
        // zip cannot — zipping them left the coarser rim coarse
        // (measured: default cone volume 1.2% low where Delaunay reads
        // 0.8%). Interpolated equalization requests from here diverged
        // (hole-rim strips: 12 → 84 split spans per round, 28–62 open
        // edges). The zip below serves only TRANSPOSED strips, whose
        // Delaunay alternative evacuates the interior outright.
        trace("end chains misaligned on a v-ruled chart");
        return None;
    }

    // Boundary stations must resolve the surface mid-chart too (a wavy
    // plate's waves are interior: the lattice flattened them; PK-650).
    let crosses = [cross_a.as_slice(), cross_b.as_slice()];
    if doubly_curved && !resolves_interior(&chart, crosses, rail_l, &eval, options.chord_height) {
        trace("boundary stations under-sample a doubly curved interior");
        return None;
    }

    // Assemble: ring points verbatim first (canonical positions →
    // watertight with the neighbours), then the interior lattice.
    let mut uv_pts: Vec<(f64, f64)> = chart.clone();
    let mut pos_pts: Vec<Option<Point3>> = ring.pos.iter().copied().map(Some).collect();
    let interior_base = uv_pts.len();
    let mut triangles: Vec<[usize; 3]>;
    let nu;
    if aligned {
        nu = cross_a.len();
        for i in 1..nu.saturating_sub(1) {
            let u = 0.5 * (chart[cross_a[i]].0 + chart[cross_b[i]].0);
            for j in 1..nv.saturating_sub(1) {
                let v = 0.5 * (chart[rail_l[j]].1 + chart[rail_r[j]].1);
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
        triangles = Vec::with_capacity(2 * (nu - 1) * (nv - 1));
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
    } else {
        // Union stations: every interior cross-chain station of either
        // end, near-coincident pairs (within half the smallest gap)
        // merged to their mean so no interior column shears a sliver.
        let mut stations: Vec<f64> = cross_a[1..cross_a.len() - 1]
            .iter()
            .chain(cross_b[1..cross_b.len() - 1].iter())
            .map(|&k| chart[k].0)
            .collect();
        stations.sort_by(f64::total_cmp);
        let mut merged: Vec<(f64, usize)> = Vec::with_capacity(stations.len());
        for s_a in stations {
            match merged.last_mut() {
                Some((mean, count)) if su * (s_a - *mean / *count as f64) <= 0.5 * gap_u => {
                    *mean += s_a;
                    *count += 1;
                }
                _ => merged.push((s_a, 1)),
            }
        }
        let stations: Vec<f64> = merged.iter().map(|&(m, c)| m / c as f64).collect();
        // Rows along the ruling: the two end chains verbatim, interior
        // rows on the union stations between their rail vertices.
        let mut rows: Vec<Vec<usize>> = Vec::with_capacity(nv);
        rows.push(cross_a.clone());
        for j in 1..nv - 1 {
            let v = 0.5 * (chart[rail_l[j]].1 + chart[rail_r[j]].1);
            let mut row = vec![rail_l[j]];
            for &u in &stations {
                row.push(uv_pts.len());
                uv_pts.push((u, v));
                pos_pts.push(None);
            }
            row.push(rail_r[j]);
            rows.push(row);
        }
        rows.push(cross_b.clone());
        nu = stations.len() + 2;
        triangles = Vec::new();
        for pair in rows.windows(2) {
            zip_rows(&pair[0], &pair[1], &uv_pts, &mut triangles);
        }
    }
    // Back to the surface chart (undo the transpose), then enforce
    // positive UV winding (the quads are built ascending in both
    // axes, but the surface chart may run either handedness).
    if swap {
        for p in &mut uv_pts {
            *p = (p.1, p.0);
        }
    }
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
    if popkern_core::trace::trace_enabled("boundary-first-ring") {
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

/// Triangulate the band between two lattice rows (vertex indices
/// ascending in the cross coordinate `a`, sharing their first and
/// last `a` at the rails): advance along whichever row's next station
/// comes first, so every triangle spans one station gap of one row.
fn zip_rows(lo: &[usize], hi: &[usize], chart: &[(f64, f64)], triangles: &mut Vec<[usize; 3]>) {
    let (mut i, mut k) = (0usize, 0usize);
    while i + 1 < lo.len() || k + 1 < hi.len() {
        let advance_lo = if k + 1 >= hi.len() {
            true
        } else if i + 1 >= lo.len() {
            false
        } else {
            chart[lo[i + 1]].0 <= chart[hi[k + 1]].0
        };
        if advance_lo {
            triangles.push([lo[i], lo[i + 1], hi[k]]);
            i += 1;
        } else {
            triangles.push([lo[i], hi[k + 1], hi[k]]);
            k += 1;
        }
    }
}
