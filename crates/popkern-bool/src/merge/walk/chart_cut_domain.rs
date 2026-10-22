//! Placement of a doubly-wrapped face's chart runs in its cut domain
//! (PK-459 stage 2; docs/architecture/torus-cut-domain-chart-walk.md
//! §3.3).
//!
//! The face's cut domain is D = [u₀, u₀ + Pu] × [v₀, v₀ + Pv], its
//! seams as borders. This pass owns the v axis: v-seam halves go to
//! the low or high v border by the transposed material-left rule, and
//! every other run is shifted by the whole v-period that brings its
//! INTERIOR into D. Endpoint anchoring is wrong on a seamed axis: a
//! node on the v-seam has two chart copies, and the band fold's node
//! projection picks one arbitrarily (measured 2026-10-22: a drill
//! ring's lower half crossing the outer equator kept both endpoints on
//! the low copy while its interior lay a period away, so its run lay
//! along the seam). The u axis keeps the u-periodic machinery
//! (`anchor_periodic`), which sees v-seam halves as ordinary runs.

use std::collections::BTreeSet;

use super::chart_seam_axes::{ChartAxis, constant_axis};
use super::chart_subdivision::FaceChartSubdivision;

/// A face walked on its cut domain: the v window and its v-seam
/// halves (the u axis is carried by the subdivision's `u_period`).
#[derive(Debug, Clone, Default)]
pub(crate) struct CutDomain {
    /// The v-seam height v₀; D spans `[v0, v0 + v_period]` in v.
    pub v0: f64,
    /// The chart v-period.
    pub v_period: f64,
    /// Registered half-edges lying on the v-seam (border copies).
    pub v_seam: BTreeSet<usize>,
}

/// Fold every embedding of `sub` into the cut domain's v window and
/// place the v-seam halves (seam-like runs constant in v) on their
/// borders. `sigma` is the face's final σ.
pub(crate) fn fold_into_cut_domain(
    sub: &mut FaceChartSubdivision,
    v0: f64,
    v_period: f64,
    u_period: f64,
    sigma: f64,
) -> CutDomain {
    let mut domain = CutDomain {
        v0,
        v_period,
        v_seam: BTreeSet::new(),
    };
    let (eps_u, eps_v) = (u_period * 1.0e-6, v_period * 1.0e-6);
    let trace = popkern_core::trace::trace_enabled("seam");
    for (&h, emb) in sub.embeddings.iter_mut() {
        let (first, last) = (emb.uv[0], emb.uv[emb.uv.len() - 1]);
        let on_v_seam = emb.seam_like && constant_axis(&emb.uv, eps_u, eps_v) == Some(ChartAxis::V);
        let shift = if on_v_seam {
            // σ = +1: the interior (v > v₀) is on the left of
            // u-increasing travel at the low border; σ = −1 mirrors.
            let low = (last.0 - first.0 > 0.0) == (sigma > 0.0);
            let target = if low { v0 } else { v0 + v_period };
            if trace {
                eprintln!(
                    "  [seam] v-seam he {h} du {:+.3} sigma {sigma:+.0} -> {}",
                    last.0 - first.0,
                    if low { "low" } else { "HIGH" }
                );
            }
            domain.v_seam.insert(h);
            ((target - first.1) / v_period).round() * v_period
        } else {
            // The sample farthest from the v-seam (mod the period)
            // anchors the run: it is the one an endpoint's roundoff
            // cannot move across a border.
            let depth = |v: f64| {
                let f = ((v - v0) / v_period).rem_euclid(1.0);
                f.min(1.0 - f)
            };
            let anchor = emb
                .uv
                .iter()
                .map(|q| q.1)
                .max_by(|a, b| depth(*a).total_cmp(&depth(*b)))
                .unwrap_or(first.1);
            let target = v0 + (anchor - v0).rem_euclid(v_period);
            ((target - anchor) / v_period).round() * v_period
        };
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
    domain
}
