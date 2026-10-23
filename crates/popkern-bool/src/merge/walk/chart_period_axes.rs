//! Chart axes of doubly-periodic faces that wrap one period: the
//! v-bounded band window (PK-774) and the transposed walk of v-wrapped
//! u-bounded bands (PK-922).
//!
//! A torus face wrapping the tube (v) through its own v-seam but bounded
//! in u (either half of a torus, a band round the tube) is, for the
//! walk, a v-closed wall: exactly a u-closed wall with the chart axes
//! swapped. The face's embeddings are transposed so the tuned u-seam
//! machinery owns it, σ flips with the transposition, and the fragment
//! charts the torus sampler reads are swapped back (`fragment_chart`).
//! Measured before: such a face got no cut-domain corner (no u-seam),
//! its sampled v window read the full-period wrap as a band, and drills
//! through the merged half torus failed validation (InnerLoopWinding,
//! LoopPeriodWinding).

use popkern_core::topology::{BRep, FaceHandle};

use super::chart_subdivision::FaceChartSubdivision;

/// The u window `(lo, hi)` of `face` (periods `(pu, pv)`) when it wraps
/// v through its own v-seam, carries no u-seam, and is bounded in u
/// strictly inside one period.
pub(super) fn v_wrapped_band(
    brep: &BRep,
    face: FaceHandle,
    (pu, pv): (f64, f64),
) -> Option<(f64, f64)> {
    let seams = super::chart_seam_axes::face_seam_axes(brep, face, (pu, pv));
    if seams.u.is_some() || seams.v.is_none() {
        return None;
    }
    let ((u0, u1), _) = brep.face_trim_uv_window(face)?;
    (u1 - u0 < pu * (1.0 - 1.0e-6)).then_some((u0, u1))
}

/// Swap the chart axes of every embedding of `sub` and flip σ. The
/// face's u (now the second axis, which the walk does not treat as
/// periodic) is folded by whole periods `pu` into the band's u
/// `window`: a record near the band's middle can embed a period off
/// (measured: a radial drill's hole record at u = 2pi on a band of
/// [-pi/2, pi/2], whose orbit then spanned the whole period).
pub(super) fn transpose(sub: &mut FaceChartSubdivision, window: (f64, f64), pu: f64) {
    let swap = |q: &mut (f64, f64)| std::mem::swap(&mut q.0, &mut q.1);
    let centre = 0.5 * (window.0 + window.1);
    for emb in sub.embeddings.values_mut() {
        emb.uv.iter_mut().for_each(swap);
        emb.refined_first.iter_mut().for_each(swap);
        emb.refined_last.iter_mut().for_each(swap);
        let mid = emb.uv[emb.uv.len() / 2].1;
        let shift = ((centre - mid) / pu).round() * pu;
        if shift != 0.0 {
            for q in emb.uv.iter_mut() {
                q.1 += shift;
            }
            for q in [&mut emb.refined_first, &mut emb.refined_last]
                .into_iter()
                .flatten()
            {
                q.1 += shift;
            }
        }
    }
    sub.sigma = -sub.sigma;
    sub.transposed = true;
}

/// The face's v window `(lo, hi)` when its boundary pcurves span
/// strictly less than the surface's v-period: the face is a v-bounded
/// band of a doubly-periodic surface, and its chart never needs the
/// v-seam. The window is the smallest CIRCULAR arc covering every
/// boundary sample (the complement of the largest gap between them,
/// mod the period), so a band whose boundary lies on the v-seam, its
/// samples split between 2π − ε and +ε by roundoff, still reads as one
/// band (PK-774).
pub(super) fn face_v_window(
    brep: &popkern_core::topology::BRep,
    face: popkern_core::topology::FaceHandle,
    v_period: Option<f64>,
) -> Option<(f64, f64)> {
    let period = v_period?;
    let boundary = crate::clip::face_uv_boundary(brep, face)?;
    let mut vs: Vec<f64> = Vec::new();
    for contour in
        std::iter::once(boundary).chain(crate::clip::face_uv_inner_boundaries(brep, face))
    {
        for (_he, pcurve) in &contour {
            let (t0, t1) = pcurve.parameter_range();
            for i in 0..=16 {
                let q = pcurve.point_at(t0 + (t1 - t0) * i as f64 / 16.0);
                if q.y.is_finite() {
                    vs.push(q.y.rem_euclid(period));
                }
            }
        }
    }
    vs.sort_by(f64::total_cmp);
    let (first, last) = (*vs.first()?, *vs.last()?);
    // The largest gap between consecutive samples, the wrap included.
    let (mut gap, mut start) = (first + period - last, first);
    for w in vs.windows(2) {
        if w[1] - w[0] > gap {
            gap = w[1] - w[0];
            start = w[1];
        }
    }
    let extent = period - gap;
    (extent < period * (1.0 - 1.0e-6)).then_some((start, start + extent))
}
