//! Phase H.8 (Option C step 3) — bi-nonic guide interior-CP
//! optimizer.
//!
//! # Problem
//!
//! The bi-quintic optimizer (phase 2) converges the 16-dim interior
//! CP block per sub-patch but cannot express cross-curvature
//! asymmetries that need degree-6+ bases along each boundary. On
//! cube_extreme (1:10 radius ratio), phase 2 bottoms out near
//! 2.19e-1 worst-case G1 coplanarity — still well above the 5e-2
//! commercial gate.
//!
//! # Approach
//!
//! Third optimization phase layered on top of phase 2's converged
//! bi-quintic guides, run as Levenberg-Marquardt nonlinear
//! least-squares over the per-sample G1 coplanarity residuals:
//!
//! 1. Run phase 1 (bicubic C-twist) and phase 2 (bi-quintic
//!    interior) via `optimize_biquintic_and_build_binonic`.
//! 2. If phase-2 residual is already below the target, return.
//! 3. Degree-elevate each bi-quintic guide to bi-nonic (pointwise-
//!    identical, so no regression baseline).
//! 4. Run L-M on the **64 interior CPs** of each bi-nonic guide
//!    (`grid[i][j]` for `i, j ∈ 1..=8`). 64 × 3 × N = 192N DoF.
//!    FD Jacobian, adaptive λ, fold-safe rejection.
//! 5. Return the final bi-nonic ring.

use super::kp_biquintic_optimizer::optimize_biquintic_and_build_binonic;
use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, KpGuidePatchBiNonic,
    RectangularGregoryBiNonic, build_binonic_ring_from_binonic_guides,
};

pub(super) struct KpBiNonicOptimizerResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    #[allow(dead_code)]
    pub final_residual: f64,
}

/// Three-phase optimizer: delegates phases 1–2 to the bi-quintic
/// optimizer, then runs phase 3 (bi-nonic interior CPs via L-M)
/// if residual is still above target.
pub(super) fn optimize_binonic_and_build_binonic(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Option<KpBiNonicOptimizerResult> {
    // Phase-2 residual above this triggers phase-3 L-M. Below it,
    // phase 2 is already in or near the commercial-gate basin (per-
    // sample coplanarity sum ≤ 1e-2 roughly corresponds to worst-
    // case G1 ≤ ~1e-2); skip the expensive L-M to keep common-case
    // test wall-time manageable.
    const PHASE3_ENGAGE_RESIDUAL: f64 = 1e-2;
    const LM_TARGET_RESIDUAL: f64 = 1e-6;

    let bq_result = optimize_biquintic_and_build_binonic(input, options)?;
    if bq_result.final_residual < PHASE3_ENGAGE_RESIDUAL {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        });
    }

    let n = bq_result.bq_guides.len();
    if n == 0 {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        });
    }

    // Elevate phase-2 converged bi-quintic guides to bi-nonic
    // (pointwise-identical at start, so no regression baseline).
    let bn_guides: Vec<KpGuidePatchBiNonic> = bq_result
        .bq_guides
        .iter()
        .map(KpGuidePatchBiNonic::from_biquintic_elevated)
        .collect();

    let initial = build_binonic_ring_from_binonic_guides(&bn_guides, input).ok()?;
    if any_self_intersects(&initial) {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        });
    }

    // Pack the 64 interior CPs × 3 axes × N patches into a flat
    // parameter vector. Order: sub-patch k, slot s ∈ 0..64, axis a.
    let x0 = pack_guides(&bn_guides);

    let eval = BinonicEval {
        base: bn_guides.clone(),
        input,
    };
    let lm_result = run_lm(
        eval,
        x0,
        &LmConfig {
            fd_eps: 1e-4,
            lambda_init: 1e-3,
            lambda_max: 1e6,
            max_iter: 15,
            target_residual: LM_TARGET_RESIDUAL,
        },
        "binonic-phase3",
    );

    // Unpack optimal parameter back into guides; rebuild ring.
    let mut best_guides = bn_guides;
    unpack_guides(&lm_result.x, &mut best_guides);
    let best_binonic = match build_binonic_ring_from_binonic_guides(&best_guides, input) {
        Ok(ring) if !any_self_intersects(&ring) => ring,
        _ => {
            return Some(KpBiNonicOptimizerResult {
                binonic: bq_result.binonic,
                final_residual: bq_result.final_residual,
            });
        }
    };

    // Guard against regression: if L-M ended worse than phase 2
    // (shouldn't happen with a well-behaved Gauss-Newton, but
    // defensive in case the sample-vector sum_sq doesn't exactly
    // match the legacy `measure_ring_residual` on the phase-2
    // output), return phase-2 result.
    let post = measure_ring_residual(&best_binonic);
    if post < bq_result.final_residual {
        Some(KpBiNonicOptimizerResult {
            binonic: best_binonic,
            final_residual: post,
        })
    } else {
        Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        })
    }
}

/// L-M residual callback for phase 3. Keeps a base copy of the
/// guide ring and a reference to the input; each `eval` clones
/// the base, overwrites interior CPs from the parameter vector,
/// rebuilds the binonic ring, and returns a per-sample coplanarity
/// vector.
struct BinonicEval<'a> {
    base: Vec<KpGuidePatchBiNonic>,
    input: &'a DirectBoundaryInput,
}

impl<'a> LmEval for BinonicEval<'a> {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        let mut guides = self.base.clone();
        unpack_guides(x, &mut guides);
        let binonic = build_binonic_ring_from_binonic_guides(&guides, self.input).ok()?;
        Some(residual_vector(&binonic))
    }

    fn fold_safe(&mut self, x: &[f64]) -> bool {
        let mut guides = self.base.clone();
        unpack_guides(x, &mut guides);
        match build_binonic_ring_from_binonic_guides(&guides, self.input) {
            Ok(binonic) => !any_self_intersects(&binonic),
            Err(_) => false,
        }
    }
}

fn pack_guides(guides: &[KpGuidePatchBiNonic]) -> Vec<f64> {
    let n = guides.len();
    let mut x = Vec::with_capacity(n * 64 * 3);
    for guide in guides {
        for slot in 0..64 {
            let (i, j) = cp_index(slot);
            x.push(guide.grid[i][j].x);
            x.push(guide.grid[i][j].y);
            x.push(guide.grid[i][j].z);
        }
    }
    x
}

fn unpack_guides(x: &[f64], guides: &mut [KpGuidePatchBiNonic]) {
    let mut idx = 0;
    for guide in guides.iter_mut() {
        for slot in 0..64 {
            let (i, j) = cp_index(slot);
            guide.grid[i][j].x = x[idx];
            guide.grid[i][j].y = x[idx + 1];
            guide.grid[i][j].z = x[idx + 2];
            idx += 3;
        }
    }
}

/// Maps FD slot index to (row, col) in the 8×8 interior block.
fn cp_index(slot: usize) -> (usize, usize) {
    let r = 1 + slot / 8;
    let c = 1 + slot % 8;
    (r, c)
}

fn any_self_intersects(binonic: &[RectangularGregoryBiNonic]) -> bool {
    for patch in binonic {
        let gp = GregoryPatch::RectangularBiNonic(Box::new(patch.clone()));
        if gp.self_intersects_with_tolerance(1e-6) {
            return true;
        }
    }
    false
}

/// Residual *vector*: per-sample coplanarity (not summed). Sample
/// count matches `measure_ring_residual`: 50 samples × N spokes.
/// Returns 0 for ill-conditioned normals so the vector length stays
/// predictable.
fn residual_vector(patches: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let n = patches.len();
    let mut out = Vec::with_capacity(n * 50);
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let cross = xtan_a.cross(&xtan_b);
            let numer = (cross.x * t_spoke.x + cross.y * t_spoke.y + cross.z * t_spoke.z).abs();
            let na = (xtan_a.x * xtan_a.x + xtan_a.y * xtan_a.y + xtan_a.z * xtan_a.z).sqrt();
            let nb = (xtan_b.x * xtan_b.x + xtan_b.y * xtan_b.y + xtan_b.z * xtan_b.z).sqrt();
            let ns = (t_spoke.x * t_spoke.x + t_spoke.y * t_spoke.y + t_spoke.z * t_spoke.z).sqrt();
            let denom = na * nb * ns;
            if denom > 1e-300 {
                out.push(numer / denom);
            } else {
                out.push(0.0);
            }
        }
    }
    out
}

fn measure_ring_residual(patches: &[RectangularGregoryBiNonic]) -> f64 {
    residual_vector(patches).iter().map(|v| v * v).sum()
}
