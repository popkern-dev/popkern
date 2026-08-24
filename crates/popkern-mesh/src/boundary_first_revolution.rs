//! Structured tessellation for full-wrap apex charts (PK-336).
//!
//! A face whose outer ring nets a FULL chart period in u cannot be
//! charted by Delaunay in one branch — the classic case is a cone
//! lateral: the base circle sweeps 2π and "returns" through the apex
//! point, whose azimuth is undefined. These faces used to decline to
//! the legacy periodic grid, whose boundary stations stagger against
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
//! anything else declines (with a `PK_BF_TRACE` reason) to the
//! pre-existing paths.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::FaceHandle;

use super::BoundaryRing;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

fn trace(reason: &str) {
    if std::env::var_os("PK_BF_TRACE").is_some() {
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
    // Constant v within chord-scale slack; u strictly monotone
    // spanning one full period.
    let v_base = ring.uv[chain[0]].1;
    let eps_v = 1e-6 * (v_base - v_apex).abs().max(1e-9);
    if chain.iter().any(|&i| (ring.uv[i].1 - v_base).abs() > eps_v) {
        trace("cross chain not at constant v");
        return None;
    }
    let u_first = ring.uv[chain[0]].0;
    let u_last = ring.uv[*chain.last()?].0;
    let dir = if u_last > u_first { 1.0 } else { -1.0 };
    for w in chain.windows(2) {
        if dir * (ring.uv[w[1]].0 - ring.uv[w[0]].0) <= 0.0 {
            trace("cross chain u not monotone");
            return None;
        }
    }
    if ((u_last - u_first).abs() - period_u).abs() > 0.05 * period_u {
        trace("cross chain does not span one period");
        return None;
    }
    // Endpoints of the chain must be weld-coincident (the seam pair).
    if ring.pos[chain[0]].distance_to(&ring.pos[*chain.last()?]) > 1e-7 * scale {
        trace("chain endpoints not weld-coincident");
        return None;
    }

    // Latitude row count from the meridian chord criterion at the
    // chart's u midpoint, doubled until the midpoint sag passes.
    let u_mid = 0.5 * (u_first + u_last);
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
    // Circumferential chord error is owned by the base chain's own
    // canonical sampling; interior rows shrink toward the apex, so
    // reusing its u stations only ever refines.

    // Assemble: apex vertex, interior latitude rows at the chain's u
    // stations, base row = ring points verbatim.
    let k = chain.len();
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
    for j in 1..rows {
        let v = v_apex + (v_base - v_apex) * j as f64 / rows as f64;
        for &i in &chain {
            push(&mut uv_pts, &mut vertices, ring.uv[i].0, v, None);
        }
    }
    let base_row = vertices.len();
    for &i in &chain {
        push(
            &mut uv_pts,
            &mut vertices,
            ring.uv[i].0,
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
    if std::env::var_os("PK_BF_TRACE").is_some() {
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
