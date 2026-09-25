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

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;

use super::super::kp_lm_solver::{LmConfig, run_lm};
use super::three_quad_kp_promote::{worst_ring_g0_gap, worst_ring_normal_sine};

#[path = "three_quad_chart_orientation_g2_jacobian.rs"]
mod three_quad_chart_orientation_g2_jacobian;
#[path = "three_quad_chart_orientation_metrics.rs"]
mod three_quad_chart_orientation_metrics;
#[path = "three_quad_chart_orientation_polish.rs"]
mod three_quad_chart_orientation_polish;
use three_quad_chart_orientation_metrics::{
    any_self_intersects, characteristic_length, worst_outer_fillet_sin, worst_ring_folded_g2,
};
use three_quad_chart_orientation_polish::{
    SeamPolishEval, boundary_tfi_ring, opposed_proximity_weights, pack_seam_rows, unpack_seam_rows,
};

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
    if pre_bands == 0 && pre_interior == 0 && !pre_folded {
        return Err("already_clean".to_string());
    }
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
    let violations = |candidate: &[RectangularGregoryBiNonic]| {
        binonic_reversal_band_count(candidate) + interior_opposed_count(candidate)
    };
    // Locality: the violations live in narrow bands, and a uniform
    // pull toward the TFI anchor destroys ring-wide coplanarity
    // (measured 7.6e-4 → 7e-1 at 1:20). Weight each CP's pull by
    // its Greville point's proximity to the opposed samples, so
    // clean patches and clean regions stay put.
    let weights = opposed_proximity_weights(ring, &references);
    let blend_to = |lambda: f64, weights: &[[[f64; 10]; 10]]| {
        let mut candidate = ring.to_vec();
        for ((patch, tfi), weight) in candidate.iter_mut().zip(anchor.iter()).zip(weights.iter()) {
            #[allow(clippy::needless_range_loop)]
            for i in 1..=8 {
                for j in 1..=8 {
                    let pull = lambda * weight[i][j];
                    if pull <= 0.0 {
                        continue;
                    }
                    let a = patch.grid[i][j];
                    let b = tfi.grid[i][j];
                    patch.grid[i][j] = Point3::new(
                        a.x + pull * (b.x - a.x),
                        a.y + pull * (b.y - a.y),
                        a.z + pull * (b.z - a.z),
                    );
                }
            }
        }
        candidate
    };
    if violations(&anchor) > 0 {
        return Err(format!(
            "boundary_tfi_still_opposed(at={})",
            describe_opposed_regions(&anchor)
        ));
    }
    // Minimal clearing blend: coarse ladder + bisection, first over
    // the localized weight field, then (only if locality cannot
    // clear) over the uniform pull.
    let uniform = [[1.0_f64; 10]; 10];
    let uniform_weights = vec![uniform; ring.len()];
    let mut chosen: Option<(f64, &Vec<[[f64; 10]; 10]>)> = None;
    for field in [&weights, &uniform_weights] {
        let mut lo = 0.0_f64;
        let mut hi = f64::NAN;
        for step in 1..=10 {
            let lambda = step as f64 / 10.0;
            if violations(&blend_to(lambda, field)) == 0 {
                hi = lambda;
                lo = lambda - 0.1;
                break;
            }
        }
        if hi.is_nan() {
            continue;
        }
        for _ in 0..5 {
            let mid = 0.5 * (lo + hi);
            if violations(&blend_to(mid, field)) == 0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        chosen = Some((hi, field));
        break;
    }
    let Some((lambda_min, field)) = chosen else {
        return Err(format!(
            "blend_cannot_clear(at={})",
            describe_opposed_regions(&blend_to(1.0, &uniform_weights))
        ));
    };
    // Start the polish with headroom past the minimal clearing
    // blend: at the feasibility boundary every descent step bumps
    // into the fold_safe rejection and the damping schedule stalls.
    let lambda = (1.25 * lambda_min).min(1.0);
    let unfolded = blend_to(lambda, field);

    // ── Stage 2: G1 polish inside the clean region ───────────
    // The blend costs spoke coplanarity and outer tangency; recover
    // them with an L-M over the near-seam rows AND the shared spoke
    // rows, whose step acceptance REQUIRES staying violation-free —
    // feasible start, feasibility maintained, so the polish can
    // never trade the orientation back for continuity. Restarted
    // while it keeps improving (the fold_safe rejections stall the
    // damping schedule well before the basin floor).
    let mut repaired = unfolded.clone();
    let mut x = pack_seam_rows(&unfolded);
    let mut previous_g1 = f64::INFINITY;
    for _round in 0..3 {
        let eval = SeamPolishEval {
            base: unfolded.clone(),
            references: references.clone(),
            // Tether to the ORIGINAL rows: the pre-repair candidate
            // had the good continuity, so pulling back toward it
            // (inside the fold_safe feasible region) is exactly the
            // recovery pressure the polish needs.
            x0: pack_seam_rows(ring),
            outer_fillets: outer_fillets.map(|f| f.to_vec()),
            char_len,
            selector_m: None,
        };
        let lm = run_lm(
            eval,
            x,
            &LmConfig {
                fd_eps: 1e-6 * char_len.max(1e-6),
                lambda_init: 1e-3,
                lambda_max: 1e9,
                max_iter: 40,
                target_residual: 1e-12,
                require_feasible_start: true,
            },
            "chart-orientation-polish",
        );
        x = lm.x;
        let mut candidate = unfolded.clone();
        unpack_seam_rows(&x, &mut candidate);
        // The solver's accepted iterates are all fold-safe, but
        // assert the invariant rather than assume it.
        if violations(&candidate) > 0 {
            break;
        }
        repaired = candidate;
        let candidate_gp: Vec<GregoryPatch> = repaired
            .iter()
            .cloned()
            .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
            .collect();
        let g1 = worst_ring_normal_sine(&candidate_gp);
        if g1 <= 5.0e-2 || g1 >= 0.95 * previous_g1 {
            break;
        }
        previous_g1 = g1;
    }

    // ── Transactional gates ─────────────────────────────────
    let post_bands = binonic_reversal_band_count(&repaired);
    let post_interior = interior_opposed_count(&repaired);
    if post_bands > 0 || post_interior > 0 {
        return Err(format!(
            "orientation_residual(bands={post_bands},interior={post_interior},at={})",
            describe_opposed_regions(&repaired)
        ));
    }
    if any_self_intersects(&repaired, 1.0e-7) {
        return Err("fold_margin".to_string());
    }
    let repaired_gp: Vec<GregoryPatch> = repaired
        .iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect();
    let post_g0 = worst_ring_g0_gap(&repaired_gp);
    if post_g0 > 1.0e-10 {
        return Err(format!("g0({post_g0:.3e})"));
    }
    // The acceptance contract is the commercial gate itself: a
    // repaired candidate must still pass 5e-2, and the recorded
    // pre/post pair keeps the quality trade attributable. A
    // relative no-regression clause is deliberately absent — any
    // orientation repair of a 1e-3-class ring costs some
    // coplanarity, and the pre-repair number is not reachable
    // inside the clean region.
    let post_g1 = worst_ring_normal_sine(&repaired_gp);
    if post_g1 > 5.0e-2 {
        return Err(format!("g1({post_g1:.3e})"));
    }
    let post_g2 = worst_ring_folded_g2(&repaired);
    if post_g2 > 5.0e-2 {
        return Err(format!("g2({post_g2:.3e})"));
    }
    if let (Some(fillets), Some(pre)) = (outer_fillets, pre_outer) {
        let post_outer = worst_outer_fillet_sin(&repaired, fillets);
        if post_outer > (1.2 * pre).max(5.0e-3) {
            return Err(format!(
                "outer_g1_regression(pre={pre:.3e},post={post_outer:.3e})"
            ));
        }
    }
    Ok(ChartRepairOutcome {
        ring: repaired,
        pre_bands,
        pre_g1,
        post_g1,
        post_g2,
    })
}

#[cfg(test)]
#[path = "three_quad_chart_orientation_tests.rs"]
mod tests;
