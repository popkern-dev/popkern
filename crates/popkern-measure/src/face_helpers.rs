//! Private helpers for face measurement.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle};

/// Compute the UV extent of a face's outer loop from its pcurves.
///
/// Returns `Some(((u0, u1), (v0, v1)))` if every half-edge in the outer loop
/// has a pcurve. Handles periodic surfaces (cylinder u, sphere u, torus u/v)
/// by unwrapping consecutive UV samples so a face that crosses the seam is
/// represented as a contiguous range rather than two fragments.
///
/// Returns `None` if any half-edge lacks a pcurve.
pub(crate) fn pcurve_uv_range(brep: &BRep, face: FaceHandle) -> Option<((f64, f64), (f64, f64))> {
    let face_data = brep.faces.get(face)?;
    let surface_on_face = face_data.surface.as_ref()?;
    let (period_u, period_v) = surface_on_face.geometry.periodicity();

    let loop_ = brep.loops.get(face_data.outer_loop)?;
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return None;
    }

    const SAMPLES: usize = 16;
    let mut all_u: Vec<f64> = Vec::new();
    let mut all_v: Vec<f64> = Vec::new();

    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pcurve = he.pcurve.as_ref()?;
        let (t0, t1) = pcurve.parameter_range();
        for i in 0..=SAMPLES {
            let t = t0 + (t1 - t0) * (i as f64 / SAMPLES as f64);
            let pt = pcurve.point_at(t);
            all_u.push(pt.x);
            all_v.push(pt.y);
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if all_u.is_empty() {
        return None;
    }

    unwrap_periodic(&mut all_u, period_u);
    unwrap_periodic(&mut all_v, period_v);

    let u_min = all_u.iter().cloned().fold(f64::INFINITY, f64::min);
    let u_max = all_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let v_min = all_v.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = all_v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if !u_min.is_finite() || !u_max.is_finite() || !v_min.is_finite() || !v_max.is_finite() {
        return None;
    }

    // On a fully-closed periodic surface, all boundary UV samples lie on the
    // seam so Δu ≈ 0 after unwrapping; in that case the face spans the full period.
    let u_range = match period_u {
        Some(p) if (u_max - u_min) < 1e-6 * p => (u_min, u_min + p),
        _ => (u_min, u_max),
    };
    let v_range = match period_v {
        Some(p) if (v_max - v_min) < 1e-6 * p => (v_min, v_min + p),
        _ => (v_min, v_max),
    };

    Some((u_range, v_range))
}

/// Unwrap periodic values in-place so consecutive entries never jump by more
/// than half the period across the seam.
pub(crate) fn unwrap_periodic(values: &mut [f64], period: Option<f64>) {
    let Some(p) = period else { return };
    if p <= 0.0 || values.len() < 2 {
        return;
    }
    let half = p * 0.5;
    for i in 1..values.len() {
        let diff = values[i] - values[i - 1];
        if diff > half {
            values[i] -= p;
        } else if diff < -half {
            values[i] += p;
        }
    }
}

/// Collect 3-D vertex positions for every half-edge in a loop.
pub(crate) fn loop_positions(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Vec<Point3> {
    let Some(loop_) = brep.loops.get(loop_handle) else {
        return Vec::new();
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }
    let mut positions = Vec::new();
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}
