//! Per-face material-sense resolution via the TRAVERSAL-COHERENT
//! normal (PK-328).
//!
//! Served boolean outputs carry MIXED per-face folded normals
//! (measured: the sphere-capped-rod union folds its cylinder wall and
//! bottom cap outward but its sphere cap center-ward, every face
//! winding-honest), so no component-wide multiplier over flag-folded
//! normals exists — the previous model moved that body's wall skin
//! OUTWARD under an "Inward" shell and served the mix structurally
//! clean at small thickness. Ray-parity probes are no better: the
//! multi-ray point classifier returns false verdicts near filleted
//! corners (PK-345, measured 4 of 8 identical corner octants).
//!
//! The robust invariant is TOPOLOGICAL: on a coherently-traversed
//! manifold, the per-face normal for which the outer loop winds CCW
//! around the face interior (the traversal-coherent normal) is
//! globally on ONE side for the whole component. Measure it per face
//! from the pcurve polygon (is the interior LEFT of traversal in
//! UV?), then fix the component's global side with the wedge integral
//! `∮ p·n_trav dA = ±3V` — positive means the traversal-coherent
//! normals face outward.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use crate::body_distance::FaceTrimEntry;

/// Wedge-integral grid resolution per axis.
const GRID: usize = 6;

/// Standard even-odd ray cast over one polygon (pnpoly).
fn polygon_contains(polygon: &[(f64, f64)], u: f64, v: f64) -> bool {
    let mut inside = false;
    let n = polygon.len();
    let mut j = n - 1;
    for i in 0..n {
        let (ui, vi) = polygon[i];
        let (uj, vj) = polygon[j];
        if ((vi > v) != (vj > v)) && (u < (uj - ui) * (v - vi) / (vj - vi) + ui) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Per-face traversal-coherence data.
pub(crate) struct FaceSense {
    /// +1 when the outer loop winds CCW in UV around the interior
    /// (traversal-coherent normal = `Su×Sv`), −1 when CW.
    pub trav_sign: f64,
    /// Un-folded wedge contribution `∮ p·(Su×Sv) du dv` over the trim.
    pub wedge: f64,
}

/// Measure one face's traversal sense and wedge contribution.
/// `None` when the left-side vote is inconclusive (degenerate pcurve
/// polygons) or the face has no usable trim entry.
pub(crate) fn face_sense(brep: &BRep, fh: FaceHandle) -> Option<FaceSense> {
    let sof = brep.faces.get(fh)?.surface.as_ref()?;
    let entry = FaceTrimEntry::build(brep, fh)?;
    let (u0, v0, u1, v1) = entry.uv_bbox()?;
    let delta = ((u1 - u0).abs() + (v1 - v0).abs()).max(1e-9) * 5e-4;

    // Left-side vote over the OUTER polygon's segments: the interior
    // is LEFT of traversal exactly when stepping along the CCW-rotated
    // segment direction stays in the trim while the CW step leaves it.
    let outer = entry.outer_polygon()?;
    let n = outer.len();
    if n < 3 {
        return None;
    }
    // Periodic-branch unwrap: a seam-straddling chart's pcurves carry
    // mixed period branches (measured: one filleted-cube corner octant
    // spans u ∈ −1.37..4.71 — width 2π — and its left-side votes
    // invert on the broken polygon). Re-branch each sample to within
    // half a period of its predecessor so the polygon is contiguous,
    // then vote against THIS polygon.
    let (period_u, period_v) = sof.geometry.parameter_periodicity();
    let unwrap = |prev: f64, x: f64, period: Option<f64>| -> f64 {
        match period {
            Some(p) if p > 0.0 => {
                let mut x = x;
                while x - prev > p * 0.5 {
                    x -= p;
                }
                while prev - x > p * 0.5 {
                    x += p;
                }
                x
            }
            _ => x,
        }
    };
    let mut outer_unwrapped: Vec<(f64, f64)> = Vec::with_capacity(n);
    for &(u, v) in outer {
        let (pu, pv) = outer_unwrapped.last().copied().unwrap_or((u, v));
        outer_unwrapped.push((unwrap(pu, u, period_u), unwrap(pv, v, period_v)));
    }
    let outer = &outer_unwrapped[..];
    // Chart-jump filter: consecutive samples of one pcurve step
    // smoothly, but the joints BETWEEN half-edges of a periodic chart
    // (seam re-entry, pole meridians) jump across the domain; a
    // left-side vote on a jump segment reads the wrong side. Skip
    // segments far above the median step.
    let mut seg_lens: Vec<f64> = (0..n)
        .map(|k| {
            let (a, b) = (outer[k], outer[(k + 1) % n]);
            ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt()
        })
        .collect();
    let mut sorted = seg_lens.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = sorted[n / 2].max(1e-12);
    let mut votes = 0i32;
    let mut conclusive = 0usize;
    for k in 0..n {
        let len = seg_lens[k];
        if len < 1e-12 || len > median * 3.0 {
            continue;
        }
        let (a, b) = (outer[k], outer[(k + 1) % n]);
        let d = (b.0 - a.0, b.1 - a.1);
        let mid = ((a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5);
        let left = (-d.1 / len * delta, d.0 / len * delta);
        let in_left = polygon_contains(outer, mid.0 + left.0, mid.1 + left.1);
        let in_right = polygon_contains(outer, mid.0 - left.0, mid.1 - left.1);
        if in_left != in_right {
            votes += if in_left { 1 } else { -1 };
            conclusive += 1;
        }
    }
    seg_lens.clear();
    // Demand a strong consensus: a chart whose polygon self-overlaps
    // (wrapped seams) votes both ways — that face falls back rather
    // than guessing.
    let agree = votes.unsigned_abs() as usize;
    if std::env::var_os("PK_SENSE_TRACE").is_some() {
        eprintln!(
            "[sense] {fh:?} votes={votes} conclusive={conclusive} outer_n={n} \
             bbox=({u0:.3},{v0:.3})..({u1:.3},{v1:.3})"
        );
    }
    if conclusive < 3 || agree * 5 < conclusive * 4 {
        return None;
    }
    let trav_sign = if votes > 0 { 1.0 } else { -1.0 };

    // Un-folded wedge over the trim (same divergence-theorem form as
    // the component resolver, before any orientation folding).
    let du = (u1 - u0) / GRID as f64;
    let dv = (v1 - v0) / GRID as f64;
    if du.abs() < f64::EPSILON || dv.abs() < f64::EPSILON {
        return None;
    }
    let mut wedge = 0.0;
    for i in 0..GRID {
        for j in 0..GRID {
            let u = u0 + du * (i as f64 + 0.5);
            let v = v0 + dv * (j as f64 + 0.5);
            if !entry.contains(u, v) {
                continue;
            }
            let p = sof.geometry.point_at(u, v);
            let su = sof.geometry.partial_u(u, v);
            let sv = sof.geometry.partial_v(u, v);
            let nvec = su.cross(&sv);
            wedge += (p.x * nvec.x + p.y * nvec.y + p.z * nvec.z) * du * dv;
        }
    }
    Some(FaceSense { trav_sign, wedge })
}

/// Resolve the per-face coefficient on the STORED surface normal
/// (`Su×Sv`) that yields the out-of-material direction, for every face
/// of the component containing `faces` (the caller groups components).
/// Returns `None` when the component's wedge is too small to call.
pub(crate) fn component_stored_normal_out_signs(
    brep: &BRep,
    component: &[FaceHandle],
    senses: &HashMap<FaceHandle, FaceSense>,
) -> Option<HashMap<FaceHandle, f64>> {
    let mut total = 0.0;
    let mut measured = 0usize;
    for fh in component {
        if let Some(s) = senses.get(fh) {
            total += s.trav_sign * s.wedge;
            measured += 1;
        }
    }
    if measured == 0 {
        return None;
    }
    let scale = brep.tolerance.linear() * 100.0;
    if total.abs() < scale * scale * scale {
        return None;
    }
    let global = if total > 0.0 { 1.0 } else { -1.0 };
    let mut out = HashMap::new();
    for fh in component {
        if let Some(s) = senses.get(fh) {
            out.insert(*fh, global * s.trav_sign);
        }
        // Unmeasured faces (degenerate polygons) fall back to the
        // caller's legacy component verdict.
    }
    Some(out)
}
