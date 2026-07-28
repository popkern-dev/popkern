//! Outer-loop pcurve UV-range derivation with periodic unwrapping —
//! split from `tessellate.rs` for the 500-line/frozen-ratchet
//! structural gate.

use popkern_core::topology::{BRep, FaceHandle};

/// Compute UV parameter range from the outer loop's pcurves.
///
/// Returns `Some((u_range, v_range))` if every half-edge in the outer loop
/// has a pcurve; `None` otherwise, signaling the caller to fall back to
/// 3D projection.
///
/// On periodic surfaces, UV values are unwrapped relative to the first
/// sample to prevent seam crossings from inflating the range to the
/// full period.
pub(super) fn pcurve_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Option<((f64, f64), (f64, f64))> {
    let face = brep.faces.get(face_handle)?;
    let surface = face.surface.as_ref()?;
    let (period_u, period_v) = surface.geometry.periodicity();

    let loop_ = brep.loops.get(face.outer_loop)?;
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return None;
    }

    let samples = 16_usize;
    let mut all_u = Vec::new();
    let mut all_v = Vec::new();

    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pcurve = he.pcurve.as_ref()?;
        let (t0, t1) = pcurve.parameter_range();
        for i in 0..=samples {
            let t = t0 + (t1 - t0) * (i as f64 / samples as f64);
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

    // Unwrap periodic values relative to the first sample.
    unwrap_periodic(&mut all_u, period_u);
    unwrap_periodic(&mut all_v, period_v);

    let u_min = all_u.iter().cloned().fold(f64::INFINITY, f64::min);
    let u_max = all_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let v_min = all_v.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = all_v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        Some(((u_min, u_max), (v_min, v_max)))
    } else {
        None
    }
}

/// Unwrap a sequence of values on a periodic domain so consecutive
/// values don't jump across the seam. Adjusts values in-place by
/// adding/subtracting the period to stay within half-period of the
/// previous value.
pub(super) fn unwrap_periodic(values: &mut [f64], period: Option<f64>) {
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
