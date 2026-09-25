//! Chart-orientation honesty for bidegree-9 binonic rings (PK-249).
//!
//! The exact-injection ring builder authors its near-boundary rows
//! with no orientation constraint (boundary-driven frames with zero
//! corner twist, Coons interior, spoke-row L-M against coplanarity
//! only); the only fold protection downstream is the arc-length
//! self-intersection screen. On lopsided rings the distorted
//! sector's rows land past the boundary rails: mildly at 1:20–1:50
//! asymmetry (chart handedness-reversal bands the screen cannot
//! see), catastrophically at 1:100 (a real fold, candidate
//! rejected). This module owns both the band metric and the
//! transactional repair: an L-M solve over the strict-interior CPs
//! with the outer boundary bit-fixed and the spoke curves moving
//! only as shared rows with pinned endpoints (ring G0 exact),
//! driven by an orientation hinge with the spoke-coplanarity field
//! as an in-solve no-regression pressure, and gated on the
//! promotion's commercial gates. (The G4 twin CP identity does NOT
//! hold on the rings this repairs — the exact-boundary spoke-row
//! L-M has already moved `row8`/`col1` independently — so the
//! repair preserves G1 the same way that optimizer does: by
//! measured coplanarity, never by re-welding.)
//! Callers keep their current candidate whenever the repair
//! declines — it can never regress a serve.
//!
//! Design: `docs/architecture/chart-orientation-honest-binonic-rings.md`.

use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;

use super::three_quad_kp_promote::worst_ring_normal_sine;

#[path = "three_quad_chart_orientation_g2_jacobian.rs"]
mod three_quad_chart_orientation_g2_jacobian;
#[path = "three_quad_chart_orientation_metrics.rs"]
mod three_quad_chart_orientation_metrics;
#[path = "three_quad_chart_orientation_polish.rs"]
mod three_quad_chart_orientation_polish;
use crate::construct::kp_spoke_normal_sine::{spoke_fold_count, spoke_senses};
use three_quad_chart_orientation_metrics::{
    any_self_intersects, characteristic_length, worst_outer_fillet_sin,
};
use three_quad_chart_orientation_polish::{boundary_tfi_ring, opposed_proximity_weights};
#[path = "three_quad_chart_orientation_stages.rs"]
mod three_quad_chart_orientation_stages;
use three_quad_chart_orientation_stages::RepairContext;

/// Count chart handedness-reversal bands over candidate patches:
/// contiguous runs of boundary-band samples whose natural normal
/// opposes the patch's mid-face normal. A healthy chart counts zero.
pub(in crate::construct) fn chart_reversal_band_count(patches: &[GregoryPatch]) -> usize {
    patches
        .iter()
        .map(|patch| match patch {
            GregoryPatch::RectangularBiNonic(binonic) => patch_band_count(binonic),
            _ => 0,
        })
        .sum()
}

/// Band count over a raw binonic ring (pre-commit candidates).
pub(in crate::construct) fn binonic_reversal_band_count(
    ring: &[RectangularGregoryBiNonic],
) -> usize {
    ring.iter().map(patch_band_count).sum()
}

fn patch_band_count(binonic: &RectangularGregoryBiNonic) -> usize {
    let normal = |u: f64, v: f64| binonic.partial_u(u, v).cross(&binonic.partial_v(u, v));
    let reference = normal(0.5, 0.5);
    let mut bands = 0usize;
    for band in 0..4 {
        let mut previous_opposed = false;
        for k in 0..=24 {
            let t = 0.02 + 0.96 * (k as f64) / 24.0;
            let (u, v) = match band {
                0 => (0.02, t),
                1 => (0.98, t),
                2 => (t, 0.02),
                _ => (t, 0.98),
            };
            let opposed = normal(u, v).dot(&reference) < 0.0;
            if opposed && !previous_opposed {
                bands += 1;
            }
            previous_opposed = opposed;
        }
    }
    bands
}

/// Compact description of where the ring's opposed rail runs sit —
/// the decline observability that localizes an unrepairable class
/// (and the raw input to any structural argument about it).
fn describe_opposed_regions(ring: &[RectangularGregoryBiNonic]) -> String {
    let mut parts: Vec<String> = Vec::new();
    for (p, binonic) in ring.iter().enumerate() {
        let normal = |u: f64, v: f64| binonic.partial_u(u, v).cross(&binonic.partial_v(u, v));
        let reference = normal(0.5, 0.5);
        for (rail, name) in [(0usize, "u0"), (1, "u1"), (2, "v0"), (3, "v1")] {
            let mut start: Option<f64> = None;
            let mut runs: Vec<String> = Vec::new();
            for k in 0..=24 {
                let t = 0.02 + 0.96 * (k as f64) / 24.0;
                let (u, v) = match rail {
                    0 => (0.02, t),
                    1 => (0.98, t),
                    2 => (t, 0.02),
                    _ => (t, 0.98),
                };
                let opposed = normal(u, v).dot(&reference) < 0.0;
                match (opposed, start) {
                    (true, None) => start = Some(t),
                    (false, Some(s)) => {
                        runs.push(format!("{s:.2}-{t:.2}"));
                        start = None;
                    }
                    _ => {}
                }
            }
            if let Some(s) = start {
                runs.push(format!("{s:.2}-0.98"));
            }
            if !runs.is_empty() {
                parts.push(format!("p{p}:{name}[{}]", runs.join(",")));
            }
        }
    }
    parts.join("|")
}

/// Opposed samples on a coarse strict-interior lattice — the repair
/// must fix the orientation field, not push the violation inward.
fn interior_opposed_count(ring: &[RectangularGregoryBiNonic]) -> usize {
    let mut opposed = 0usize;
    for binonic in ring {
        let normal = |u: f64, v: f64| binonic.partial_u(u, v).cross(&binonic.partial_v(u, v));
        let reference = normal(0.5, 0.5);
        for i in 0..9 {
            let u = 0.1 + 0.8 * (i as f64) / 8.0;
            for j in 0..9 {
                let v = 0.1 + 0.8 * (j as f64) / 8.0;
                if normal(u, v).dot(&reference) < 0.0 {
                    opposed += 1;
                }
            }
        }
    }
    opposed
}

/// Accepted repair: the orientation-clean ring plus the measured
/// before/after quality numbers for the diagnostics fact.
pub(in crate::construct) struct ChartRepairOutcome {
    pub ring: Vec<RectangularGregoryBiNonic>,
    pub pre_bands: usize,
    pub pre_g1: f64,
    pub post_g1: f64,
    pub post_g2: f64,
}

/// Repair a binonic ring's chart orientation with the outer
/// boundary held bit-fixed. `Err` carries the decline reason for
/// the diagnostics fact; on `Err` the caller keeps its candidate.
///
/// `outer_fillets`, when provided, are the per-sub-patch support
/// surfaces along each patch's `u=0` boundary; the repair holds the
/// cap-vs-fillet G1 there under a no-regression gate because the
/// row-1 frame is a degree of freedom.
pub(in crate::construct) fn repair_chart_orientation(
    ring: &[RectangularGregoryBiNonic],
    outer_fillets: Option<&[SurfaceGeometry]>,
) -> Result<ChartRepairOutcome, String> {
    let n = ring.len();
    if n == 0 {
        return Err("empty_ring".to_string());
    }
    let pre_bands = binonic_reversal_band_count(ring);
    let pre_interior = interior_opposed_count(ring);
    let pre_folded = any_self_intersects(ring, 1.0e-6);
    let char_len = characteristic_length(ring);
    let mut references = Vec::with_capacity(n);
    for binonic in ring {
        let reference = binonic
            .partial_u(0.5, 0.5)
            .cross(&binonic.partial_v(0.5, 0.5));
        let len = reference.length();
        // NaN-hostile comparison: a NaN length must decline, so the
        // accept path requires a definite `>`.
        if len.partial_cmp(&(1e-9 * char_len * char_len)) != Some(std::cmp::Ordering::Greater) {
            return Err("degenerate_reference_normal".to_string());
        }
        references.push((reference * (1.0 / len), len));
    }
    // Folds across a spoke are violations like in-patch reversals
    // (PK-544). The orientation-free spoke sine reads a 180° fold as
    // G1, so the sense each spoke must keep comes from the patches'
    // reference normals.
    let senses = spoke_senses(&references);
    let pre_spoke_folds = spoke_fold_count(ring, &senses);
    if pre_bands == 0 && pre_interior == 0 && pre_spoke_folds == 0 && !pre_folded {
        return Err("already_clean".to_string());
    }

    let ring_gp: Vec<GregoryPatch> = ring
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect();
    let pre_g1 = worst_ring_normal_sine(&ring_gp);
    let pre_outer = outer_fillets.map(|f| worst_outer_fillet_sin(ring, f));

    // ── Stage 1: constructive unfold ─────────────────────────
    // The transfinite Coons fill of the four (healthy, fixed)
    // boundary rows is the orientation-clean anchor: blend every
    // strict-interior CP from the current candidate toward it just
    // far enough to clear all violations. Deterministic; when even
    // the full fill stays opposed, the boundary itself forces the
    // fold and the decline records where — the structural-argument
    // measurement.
    let anchor = boundary_tfi_ring(ring);
    let in_patch = |candidate: &[RectangularGregoryBiNonic]| {
        binonic_reversal_band_count(candidate) + interior_opposed_count(candidate)
    };
    if in_patch(&anchor) > 0 {
        return Err(format!(
            "boundary_tfi_still_opposed(at={})",
            describe_opposed_regions(&anchor)
        ));
    }
    if spoke_fold_count(&anchor, &senses) > 0 {
        return Err("boundary_tfi_spoke_fold".to_string());
    }
    // Locality: the violations live in narrow bands, and a uniform
    // pull toward the TFI anchor destroys ring-wide coplanarity
    // (measured 7.6e-4 → 7e-1 at 1:20). Weight each CP's pull by
    // its Greville point's proximity to the opposed samples, so
    // clean patches and clean regions stay put.
    let weights = opposed_proximity_weights(ring, &references);
    let context = RepairContext {
        ring,
        anchor: &anchor,
        weights: &weights,
        references: &references,
        senses: &senses,
        outer_fillets,
        pre_outer,
        char_len,
    };
    // Two unfold strategies, each polished and gated in full (PK-544).
    // The minimal unfold clears only in-patch reversals and leaves
    // spoke folds to the polish's oriented objective. The pull toward
    // the anchor costs continuity, so the smallest pull serves best
    // (cube_extreme 1:40 spoke sine 3.9e-2, against 4.7e-2 after a
    // spoke-clearing pull). The spoke-clearing unfold also clears
    // spoke folds, and runs when the minimal one declines.
    let minimal = context
        .unfold(&|c| in_patch(c))
        .and_then(|u| context.polish_and_gate(u));
    let outcome = match minimal {
        Ok(outcome) => outcome,
        Err(minimal_reason) => context
            .unfold(&|c| in_patch(c) + spoke_fold_count(c, &senses))
            .and_then(|u| context.polish_and_gate(u))
            .map_err(|reason| format!("minimal[{minimal_reason}];spoke_clearing[{reason}]"))?,
    };
    Ok(ChartRepairOutcome {
        ring: outcome.ring,
        pre_bands,
        pre_g1,
        post_g1: outcome.post_g1,
        post_g2: outcome.post_g2,
    })
}

#[cfg(test)]
#[path = "three_quad_chart_orientation_tests.rs"]
mod tests;
