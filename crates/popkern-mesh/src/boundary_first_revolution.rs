//! Structured tessellation for full-wrap apex charts (PK-336).
//!
//! A face whose outer ring nets a FULL chart period in u cannot be
//! charted by Delaunay in one branch — the classic case is a cone
//! lateral: the base circle sweeps 2π and "returns" through the apex
//! point, whose azimuth is undefined. These faces used to decline to
//! the periodic uniform grid, whose boundary stations stagger against
//! canonical neighbours; the meshes stayed watertight only through
//! the crack zipper's parameter merge (measured: 64 repairs on the
//! primitive cone, 112 on the shelled cone).
//!
//! The structured revolution mesh serves the class by construction:
//! the base cross chain is the canonical ring VERBATIM (bit-identical
//! weld with the neighbouring cap), the apex is one vertex, and the
//! interior is a lattice of latitude rows at the base chain's u
//! stations, row count chosen by the meridian chord criterion. The
//! decomposition is strict — one full-period non-degenerate cross
//! chain at constant v plus one chart-degenerate apex run — and
//! anything else declines (with a the `boundary-first` trace topic reason) to the
//! pre-existing paths.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::FaceHandle;

use super::BoundaryRing;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

fn trace(reason: &str) {
    if popkern_core::trace::trace_enabled("boundary-first") {
        eprintln!("    [bf] full-wrap revolution declined: {reason}");
    }
}

/// Attempt the structured revolution tessellation of a full-u-wrap
/// ring. `u_degenerate` marks ring samples whose azimuth is
/// chart-degenerate (computed by the caller's unwrap preamble).
pub(super) fn try_tessellate(
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    orientation: bool,
    ring: &BoundaryRing,
    u_degenerate: &[bool],
    period_u: f64,
    options: &TessellationOptions,
) -> Option<FaceMesh> {
    let n = ring.uv.len();
    if n < 4 || u_degenerate.len() != n {
        return None;
    }
    // The apex run and the cross chain must each be one contiguous
    // cyclic block.
    let transitions = (0..n)
        .filter(|&i| u_degenerate[i] != u_degenerate[(i + 1) % n])
        .count();
    let apex_count = u_degenerate.iter().filter(|&&d| d).count();
    if transitions == 4 && apex_count > 0 && apex_count < n {
        // Two chart-degenerate runs: the two-apex class (full spheres,
        // PK-339).
        return super::boundary_first_revolution_bipolar::try_tessellate(
            face_handle,
            surface,
            orientation,
            ring,
            u_degenerate,
            period_u,
            options,
        );
    }
    if apex_count == 0 || apex_count == n || transitions != 2 {
        trace("no single apex run");
        return None;
    }
    // Apex: every degenerate sample must be the SAME point.
    let apex_indices: Vec<usize> = (0..n).filter(|&i| u_degenerate[i]).collect();
    let apex_pos = ring.pos[apex_indices[0]];
    let scale = ring
        .pos
        .iter()
        .zip(ring.pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max)
        .max(1e-9);
    if apex_indices
        .iter()
        .any(|&i| ring.pos[i].distance_to(&apex_pos) > 1e-7 * scale)
    {
        trace("degenerate run is not a single apex point");
        return None;
    }
    let v_apex =
        apex_indices.iter().map(|&i| ring.uv[i].1).sum::<f64>() / apex_indices.len() as f64;

    // Cross chain: the non-degenerate block in ring order, starting
    // just after the apex run.
    let start = (0..n).find(|&i| u_degenerate[i] && !u_degenerate[(i + 1) % n])? + 1;
    let chain: Vec<usize> = (0..n)
        .map(|k| (start + k) % n)
        .filter(|&i| !u_degenerate[i])
        .collect();
    if chain.len() < 4 {
        return None;
    }
    // PK-339: a CURVED meridian seam chord-refines into interior
    // stations, so the chain may open and close with meridian LADDER
    // runs (u-constant columns between apex and base) around the
    // constant-v circle block. Straight seams (cone laterals) carry no
    // interior stations — empty ladders keep the original shape (and
    // the original row generation, bit-identical).
    let v_base = ring.uv[chain[chain.len() / 2]].1;
    let eps_v = 1e-6 * (v_base - v_apex).abs().max(1e-9);
    let eps_u = 1e-6 * period_u;
    let off_base = |i: usize| (ring.uv[i].1 - v_base).abs() > eps_v;
    let lead: Vec<usize> = chain.iter().copied().take_while(|&i| off_base(i)).collect();
    let trail_rev: Vec<usize> = chain
        .iter()
        .rev()
        .copied()
        .take_while(|&i| off_base(i))
        .collect();
    if lead.len() + trail_rev.len() >= chain.len() {
        trace("no constant-v circle block between the meridian ladders");
        return None;
    }
    let mut middle: Vec<usize> = chain[lead.len()..chain.len() - trail_rev.len()].to_vec();
    if middle.len() < 4 {
        trace("circle block too short");
        return None;
    }
    if middle.iter().any(|&i| off_base(i)) {
        trace("cross chain not at constant v");
        return None;
    }
    // A rail split into multiple ARC edges (a strip-anchored rim,
    // PK-344) repeats its shared anchor position at each arc joint.
    // Dedupe only BIT-IDENTICAL repeats: near-identical pairs must BOTH
    // stay so the position-welded boundary matches the neighbor's
    // samples (a distance-based dedupe dropped one side of a 1e-9
    // anchor pair and left an unpaired triangle fan at the joint).
    middle.dedup_by(|&mut b, &mut a| {
        let (pa, pb) = (ring.pos[a], ring.pos[b]);
        pa.x.to_bits() == pb.x.to_bits()
            && pa.y.to_bits() == pb.y.to_bits()
            && pa.z.to_bits() == pb.z.to_bits()
    });
    // Chain-continuous unwrap of the block's azimuths: multi-arc rails
    // sample mixed period branches (measured: two half-rails at
    // (−π, π] put a +2π jump mid-block and the whole-wrap arm declined,
    // dropping the cap to a fallback that left 77 open cracks along
    // the rail). Identity for jump-free single-edge rails, so pre-existing
    // charts stay bit-identical.
    let u_unwrapped: Vec<f64> = {
        let mut us = Vec::with_capacity(middle.len());
        let mut prev = ring.uv[middle[0]].0;
        us.push(prev);
        for &i in middle.iter().skip(1) {
            let raw = ring.uv[i].0;
            let mut delta = (raw - prev).rem_euclid(period_u);
            if delta > period_u * 0.5 {
                delta -= period_u;
            }
            prev += delta;
            us.push(prev);
        }
        us
    };
    let u_first = u_unwrapped[0];
    let u_last = *u_unwrapped.last()?;
    // Ladders: u-constant at their seam azimuth, mirrored v stations,
    // v monotone from the apex toward the base.
    if lead.len() != trail_rev.len() {
        trace("meridian ladders differ in station count");
        return None;
    }
    // Azimuth comparisons are periodic: the unwrap phase may leave a
    // ladder on the far branch (u = period while its circle end reads
    // 0 — the same seam azimuth).
    let periodic_du = |a: f64, b: f64| {
        let d = (a - b).rem_euclid(period_u);
        d.min(period_u - d)
    };
    for &i in &lead {
        if periodic_du(ring.uv[i].0, u_first) > eps_u {
            trace("lead meridian ladder leaves the seam azimuth");
            return None;
        }
    }
    for &i in &trail_rev {
        if periodic_du(ring.uv[i].0, u_last) > eps_u {
            trace("trail meridian ladder leaves the seam azimuth");
            return None;
        }
    }
    for (a, b) in lead.iter().zip(trail_rev.iter()) {
        if (ring.uv[*a].1 - ring.uv[*b].1).abs() > eps_v {
            trace("meridian ladders not mirrored");
            return None;
        }
    }
    let toward_base = (v_base - v_apex).signum();
    for w in lead.windows(2) {
        if toward_base * (ring.uv[w[1]].1 - ring.uv[w[0]].1) <= 0.0 {
            trace("meridian ladder v not monotone");
            return None;
        }
    }
    // Circle block: u monotone spanning one full period (checked on
    // the unwrapped azimuths; a surviving near-duplicate anchor pair
    // steps by ~0 and is not a violation).
    let dir = if u_last > u_first { 1.0 } else { -1.0 };
    for w in u_unwrapped.windows(2) {
        if dir * (w[1] - w[0]) < -eps_u {
            trace("cross chain u not monotone");
            return None;
        }
    }
    if ((u_last - u_first).abs() - period_u).abs() > 0.05 * period_u {
        trace("cross chain does not span one period");
        return None;
    }
    // Endpoints of the circle block must be weld-coincident (the seam
    // pair at the base latitude).
    if ring.pos[middle[0]].distance_to(&ring.pos[*middle.last()?]) > 1e-7 * scale {
        trace("chain endpoints not weld-coincident");
        return None;
    }

    // Interior latitude rows. Meridian ladders present: their v
    // stations ARE the rows (the seam chain is already chord-refined;
    // reusing it keeps the face boundary exactly equal to the ring —
    // PK-339). Straight seams: row count from the meridian chord
    // criterion at the chart's u midpoint, doubled until the midpoint
    // sag passes (the original arm, bit-identical).
    let u_mid = 0.5 * (u_first + u_last);
    // Each interior row carries its v and, when it came from the ring,
    // the lead/trail ladder index whose positions the seam columns
    // adopt verbatim.
    let interior_vs: Vec<(f64, Option<usize>)> = if lead.is_empty() {
        let mut rows = 1usize;
        while rows < options.max_subdivisions.max(1) {
            let mut worst = 0.0_f64;
            for j in 0..rows {
                let va = v_apex + (v_base - v_apex) * j as f64 / rows as f64;
                let vb = v_apex + (v_base - v_apex) * (j + 1) as f64 / rows as f64;
                let pa = surface.point_at(u_mid, va);
                let pb = surface.point_at(u_mid, vb);
                let pm = surface.point_at(u_mid, 0.5 * (va + vb));
                let chord_mid = popkern_core::geometry::Point3::new(
                    0.5 * (pa.x + pb.x),
                    0.5 * (pa.y + pb.y),
                    0.5 * (pa.z + pb.z),
                );
                worst = worst.max(pm.distance_to(&chord_mid));
            }
            if worst <= options.chord_height {
                break;
            }
            rows *= 2;
        }
        (1..rows)
            .map(|j| (v_apex + (v_base - v_apex) * j as f64 / rows as f64, None))
            .collect()
    } else {
        // Ladder rows, sag-refined: the seam chain bounds the CURVE's
        // chord error but a coarse ladder can under-resolve the
        // surface between rows (measured: the sphere-capped rod's cap
        // volume fell 1.3% below its analytic pin). Bisect any band
        // whose mid-band sag at u_mid exceeds the chord height;
        // inserted rows are surface evaluations on BOTH seam sides,
        // so welds stay closed.
        let mut rows: Vec<(f64, Option<usize>)> = lead
            .iter()
            .enumerate()
            .map(|(j, &i)| (ring.uv[i].1, Some(j)))
            .collect();
        let sag = |va: f64, vb: f64| {
            let pa = surface.point_at(u_mid, va);
            let pb = surface.point_at(u_mid, vb);
            let pm = surface.point_at(u_mid, 0.5 * (va + vb));
            let chord_mid = popkern_core::geometry::Point3::new(
                0.5 * (pa.x + pb.x),
                0.5 * (pa.y + pb.y),
                0.5 * (pa.z + pb.z),
            );
            pm.distance_to(&chord_mid)
        };
        for _ in 0..options.max_subdivisions.max(1) {
            let mut refined: Vec<(f64, Option<usize>)> = Vec::with_capacity(rows.len() * 2);
            let mut inserted = false;
            let mut prev_v = v_apex;
            for &(v, src) in &rows {
                if sag(prev_v, v) > options.chord_height {
                    refined.push((0.5 * (prev_v + v), None));
                    inserted = true;
                }
                refined.push((v, src));
                prev_v = v;
            }
            if sag(prev_v, v_base) > options.chord_height {
                refined.push((0.5 * (prev_v + v_base), None));
                inserted = true;
            }
            rows = refined;
            if !inserted {
                break;
            }
        }
        rows
    };
    let rows = interior_vs.len() + 1;
    // Circumferential chord error is owned by the base chain's own
    // canonical sampling; interior rows shrink toward the apex, so
    // reusing its u stations only ever refines.

    // Assemble: apex vertex, interior latitude rows at the circle
    // block's u stations, base row = ring points verbatim. Seam-column
    // interior points take the ring's LADDER positions verbatim so the
    // face boundary equals the ring exactly.
    let k = middle.len();
    let mut uv_pts: Vec<(f64, f64)> = Vec::with_capacity(1 + (rows - 1) * k + k);
    let mut vertices: Vec<MeshVertex> = Vec::new();
    let push = |uv_pts: &mut Vec<(f64, f64)>,
                vertices: &mut Vec<MeshVertex>,
                u: f64,
                v: f64,
                pos: Option<popkern_core::geometry::Point3>| {
        let position = pos.unwrap_or_else(|| surface.point_at(u, v));
        let mut normal = surface.normal_at(u, v);
        let len = normal.length();
        if len > 1e-300 {
            normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
        }
        if !orientation {
            normal = -normal;
        }
        uv_pts.push((u, v));
        vertices.push(MeshVertex {
            position,
            normal,
            uv: (u, v),
        });
    };
    // Apex normal is chart-degenerate; evaluate slightly off-apex.
    push(
        &mut uv_pts,
        &mut vertices,
        u_mid,
        v_apex + (v_base - v_apex) * 1e-6,
        Some(apex_pos),
    );
    for &(v, src) in &interior_vs {
        for (idx, &u) in u_unwrapped.iter().enumerate() {
            let seam_pos = match src {
                Some(j) if idx == 0 => Some(ring.pos[lead[j]]),
                Some(j) if idx == k - 1 => Some(ring.pos[trail_rev[j]]),
                _ => None,
            };
            push(&mut uv_pts, &mut vertices, u, v, seam_pos);
        }
    }
    let base_row = vertices.len();
    // Base row: ring positions VERBATIM (the weld contract), but the
    // stored uv uses the UNWRAPPED azimuths — the per-triangle UV
    // winding pass below compares base and interior rows, and a
    // mixed-branch base u at the arc joint made one band triangle's
    // signed area lie, flipping it into an overlapped diagonal
    // (PK-348: the last unpaired sliver at the strip-anchored rim).
    for (idx, &i) in middle.iter().enumerate() {
        push(
            &mut uv_pts,
            &mut vertices,
            u_unwrapped[idx],
            ring.uv[i].1,
            Some(ring.pos[i]),
        );
    }
    let row_at = |j: usize, i: usize| -> usize {
        if j == 0 {
            0
        } else if j == rows {
            base_row + i
        } else {
            1 + (j - 1) * k + i
        }
    };
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    // Apex fan to the first ring of stations.
    for i in 0..k - 1 {
        triangles.push([0, row_at(1, i), row_at(1, i + 1)]);
    }
    // Latitude bands (k stations, k−1 segments — the chain's first and
    // last stations are the weld-coincident seam pair, so no extra
    // wrap segment is needed).
    for j in 1..rows {
        for i in 0..k - 1 {
            let (p00, p10, p11, p01) = (
                row_at(j, i),
                row_at(j, i + 1),
                row_at(j + 1, i + 1),
                row_at(j + 1, i),
            );
            triangles.push([p00, p10, p11]);
            triangles.push([p00, p11, p01]);
        }
    }
    // Positive UV winding, then the stored orientation flag (the
    // global winding-coherence pass remains the final authority).
    for t in &mut triangles {
        let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        let signed = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if signed < 0.0 {
            t.swap(1, 2);
        }
    }
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    if popkern_core::trace::trace_enabled("boundary-first") {
        eprintln!(
            "    [bf] face {face_handle:?} served as full-wrap revolution: {k} stations x {rows} rows"
        );
    }
    Some(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}
