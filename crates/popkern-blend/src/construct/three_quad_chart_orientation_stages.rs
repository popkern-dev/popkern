//! The chart-orientation repair's stages: the constructive unfold
//! toward the boundary TFI anchor, the G1 polish inside the clean
//! region, and the transactional gates (PK-249, restructured under
//! PK-544 so one context can run two unfold strategies).

use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;

use super::super::super::kp_lm_solver::{LmConfig, run_lm};
use super::super::three_quad_kp_promote::{worst_ring_g0_gap, worst_ring_normal_sine};
use super::three_quad_chart_orientation_metrics::{
    any_self_intersects, worst_outer_fillet_sin, worst_ring_folded_g2,
};
use super::three_quad_chart_orientation_polish::{
    SeamPolishEval, pack_seam_rows, unpack_seam_rows,
};
use super::{binonic_reversal_band_count, describe_opposed_regions, interior_opposed_count};
use crate::construct::kp_spoke_normal_sine::spoke_fold_count;

/// Per-patch pull weights toward the anchor, one per control point.
type PullWeights = Vec<[[f64; 10]; 10]>;

/// Everything the repair stages read, fixed for the whole repair.
pub(super) struct RepairContext<'a> {
    /// The candidate ring being repaired (its outer boundary is held).
    pub ring: &'a [RectangularGregoryBiNonic],
    /// The boundary TFI fill: orientation-clean by precondition.
    pub anchor: &'a [RectangularGregoryBiNonic],
    /// Localized pull weights (proximity to the opposed samples).
    pub weights: &'a PullWeights,
    /// Per patch: fixed mid-face unit normal and its magnitude.
    pub references: &'a [(Vec3, f64)],
    /// Per spoke: the relative chart orientation it must keep.
    pub senses: &'a [bool],
    /// Support surfaces along each patch's `u = 0` boundary.
    pub outer_fillets: Option<&'a [SurfaceGeometry]>,
    /// Worst cap-vs-fillet sine before the repair.
    pub pre_outer: Option<f64>,
    /// Characteristic ring length (tether and FD scales).
    pub char_len: f64,
}

/// A repaired ring that passed every gate.
pub(super) struct GatedRepair {
    pub ring: Vec<RectangularGregoryBiNonic>,
    pub post_g1: f64,
    pub post_g2: f64,
}

impl RepairContext<'_> {
    fn blend_to(&self, lambda: f64, weights: &PullWeights) -> Vec<RectangularGregoryBiNonic> {
        let mut candidate = self.ring.to_vec();
        for ((patch, tfi), weight) in candidate.iter_mut().zip(self.anchor).zip(weights) {
            #[allow(clippy::needless_range_loop)]
            for i in 1..=8 {
                for j in 1..=8 {
                    let pull = lambda * weight[i][j];
                    if pull <= 0.0 {
                        continue;
                    }
                    let (a, b) = (patch.grid[i][j], tfi.grid[i][j]);
                    patch.grid[i][j] = Point3::new(
                        a.x + pull * (b.x - a.x),
                        a.y + pull * (b.y - a.y),
                        a.z + pull * (b.z - a.z),
                    );
                }
            }
        }
        candidate
    }

    /// The minimal blend toward the anchor that clears `violations`:
    /// a coarse ladder plus bisection, first over the localized weight
    /// field, then (only if locality cannot clear) over the uniform
    /// pull. The start carries headroom past the minimal clearing
    /// blend: at the feasibility boundary every descent step bumps
    /// into the polish's fold_safe rejection and its damping stalls.
    pub fn unfold(
        &self,
        violations: &dyn Fn(&[RectangularGregoryBiNonic]) -> usize,
    ) -> Result<Vec<RectangularGregoryBiNonic>, String> {
        let uniform: PullWeights = vec![[[1.0_f64; 10]; 10]; self.ring.len()];
        for field in [self.weights, &uniform] {
            let mut lo = 0.0_f64;
            let mut hi = f64::NAN;
            for step in 1..=10 {
                let lambda = step as f64 / 10.0;
                if violations(&self.blend_to(lambda, field)) == 0 {
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
                if violations(&self.blend_to(mid, field)) == 0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return Ok(self.blend_to((1.25 * hi).min(1.0), field));
        }
        Err(format!(
            "blend_cannot_clear(at={})",
            describe_opposed_regions(&self.blend_to(1.0, &uniform))
        ))
    }

    /// Stage 2, the G1 polish inside the clean region, then the gates.
    ///
    /// The blend costs spoke continuity and outer tangency. An L-M over
    /// the near-seam rows and the shared spoke rows recovers them. The
    /// objective is the oriented spoke normal chord, which pushes a
    /// fold toward G1, never deeper (PK-544). Step acceptance keeps
    /// every in-patch violation cleared and never adds a spoke fold, so
    /// the polish cannot trade orientation back for continuity. It is
    /// restarted while it keeps improving, since fold_safe rejections
    /// stall the damping schedule well before the basin floor.
    pub fn polish_and_gate(
        &self,
        unfolded: Vec<RectangularGregoryBiNonic>,
    ) -> Result<GatedRepair, String> {
        let start_spoke_folds = spoke_fold_count(&unfolded, self.senses);
        let mut repaired = unfolded.clone();
        let mut x = pack_seam_rows(&unfolded);
        let mut previous = f64::INFINITY;
        for _round in 0..3 {
            let eval = SeamPolishEval {
                base: unfolded.clone(),
                references: self.references.to_vec(),
                senses: self.senses.to_vec(),
                max_spoke_folds: start_spoke_folds,
                // Tether to the ORIGINAL rows: the pre-repair candidate
                // had the good continuity, so pulling back toward it
                // (inside the feasible region) is the recovery pressure
                // the polish needs.
                x0: pack_seam_rows(self.ring),
                outer_fillets: self.outer_fillets.map(|f| f.to_vec()),
                char_len: self.char_len,
                selector_m: None,
            };
            let lm = run_lm(
                eval,
                x,
                &LmConfig {
                    fd_eps: 1e-6 * self.char_len.max(1e-6),
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
            if binonic_reversal_band_count(&candidate) + interior_opposed_count(&candidate) > 0 {
                break;
            }
            repaired = candidate;
            let folds = spoke_fold_count(&repaired, self.senses);
            let g1 = worst_ring_normal_sine(&as_patches(&repaired));
            let measure = if folds > 0 { 1.0 + folds as f64 } else { g1 };
            if measure <= 5.0e-2 || measure >= 0.95 * previous {
                break;
            }
            previous = measure;
        }
        self.gate(repaired)
    }

    /// The transactional gates: the promotion's commercial gates on the
    /// repaired ring, plus zero spoke folds against the reference
    /// senses.
    fn gate(&self, repaired: Vec<RectangularGregoryBiNonic>) -> Result<GatedRepair, String> {
        let post_bands = binonic_reversal_band_count(&repaired);
        let post_interior = interior_opposed_count(&repaired);
        if post_bands > 0 || post_interior > 0 {
            return Err(format!(
                "orientation_residual(bands={post_bands},interior={post_interior},at={})",
                describe_opposed_regions(&repaired)
            ));
        }
        let post_spoke_folds = spoke_fold_count(&repaired, self.senses);
        if post_spoke_folds > 0 {
            return Err(format!("spoke_fold(samples={post_spoke_folds})"));
        }
        if any_self_intersects(&repaired, 1.0e-7) {
            return Err("fold_margin".to_string());
        }
        let repaired_gp = as_patches(&repaired);
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
        if let (Some(fillets), Some(pre)) = (self.outer_fillets, self.pre_outer) {
            let post_outer = worst_outer_fillet_sin(&repaired, fillets);
            if post_outer > (1.2 * pre).max(5.0e-3) {
                return Err(format!(
                    "outer_g1_regression(pre={pre:.3e},post={post_outer:.3e})"
                ));
            }
        }
        Ok(GatedRepair {
            ring: repaired,
            post_g1,
            post_g2,
        })
    }
}

fn as_patches(ring: &[RectangularGregoryBiNonic]) -> Vec<GregoryPatch> {
    ring.iter()
        .cloned()
        .map(|p| GregoryPatch::RectangularBiNonic(Box::new(p)))
        .collect()
}
