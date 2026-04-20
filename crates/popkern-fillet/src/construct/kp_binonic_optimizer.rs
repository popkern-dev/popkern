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
//! Third FD descent phase layered on top of phase 2's converged
//! bi-quintic guides:
//!
//! 1. Run phase 1 (bicubic C-twist) and phase 2 (bi-quintic
//!    interior) via `optimize_biquintic_and_build_binonic`.
//! 2. If phase-2 residual is already below the target, return.
//! 3. Degree-elevate each bi-quintic guide to bi-nonic (pointwise-
//!    identical, so no regression baseline).
//! 4. FD gradient-descent on the **64 interior CPs** of each bi-
//!    nonic guide (`grid[i][j]` for `i, j ∈ 1..=8`). 64 × 3 × N =
//!    192N DoF — 4× the phase-2 cost but only invoked when phase 2
//!    fails to converge.
//! 5. Return the final bi-nonic ring.

use super::kp_biquintic_optimizer::optimize_biquintic_and_build_binonic;
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
/// optimizer, then runs phase 3 (bi-nonic interior CPs) if residual
/// is still above target.
pub(super) fn optimize_binonic_and_build_binonic(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Option<KpBiNonicOptimizerResult> {
    const TARGET_RESIDUAL: f64 = 1e-4;
    const MAX_ITER: usize = 100;
    const ALPHA_INIT: f64 = 1e-2;
    const ALPHA_MIN: f64 = 1e-8;
    const FD_EPS: f64 = 1e-4;

    // Phase 1 + 2: delegate to existing bi-quintic optimizer.
    let bq_result = optimize_biquintic_and_build_binonic(input, options)?;
    if bq_result.final_residual < TARGET_RESIDUAL {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        });
    }

    // Elevate phase-2 converged bi-quintic guides to bi-nonic
    // (pointwise-identical at start, so no regression baseline).
    let n = bq_result.bq_guides.len();
    if n == 0 {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        });
    }
    let mut bn_guides: Vec<KpGuidePatchBiNonic> = bq_result
        .bq_guides
        .iter()
        .map(KpGuidePatchBiNonic::from_biquintic_elevated)
        .collect();

    let initial = build_binonic_ring_from_binonic_guides(&bn_guides, input).ok()?;
    if any_self_intersects(&initial) {
        // Initial bi-nonic build self-intersects — fall back to the
        // phase-2 result (which was not folded).
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        });
    }
    let mut best_residual = measure_ring_residual(&initial);
    let mut best_binonic = initial;
    let mut alpha = ALPHA_INIT;

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!("[kp_binonic N={n}] phase3 init: residual={best_residual:.3e}");
    }

    let mut iter_count = 0usize;
    let mut accepted_count = 0usize;
    for _iter in 0..MAX_ITER {
        iter_count += 1;
        if best_residual < TARGET_RESIDUAL {
            break;
        }

        // FD gradient: perturb each interior guide CP along each
        // axis, rebuild, measure. N × 64 × 3 = 192N partials.
        let mut grad: Vec<Vec<[f64; 3]>> = (0..n).map(|_| vec![[0.0; 3]; 64]).collect();
        #[allow(clippy::needless_range_loop)]
        for k in 0..n {
            for slot in 0..64 {
                let (i, j) = cp_index(slot);
                for (axis, out) in grad[k][slot].iter_mut().enumerate() {
                    let mut bn_perturbed = bn_guides.clone();
                    let p = &mut bn_perturbed[k].grid[i][j];
                    match axis {
                        0 => p.x += FD_EPS,
                        1 => p.y += FD_EPS,
                        _ => p.z += FD_EPS,
                    }
                    let Ok(binonic) = build_binonic_ring_from_binonic_guides(&bn_perturbed, input)
                    else {
                        continue;
                    };
                    let r = measure_ring_residual(&binonic);
                    *out = (r - best_residual) / FD_EPS;
                }
            }
        }

        let grad_norm: f64 = grad
            .iter()
            .flatten()
            .flatten()
            .map(|g| g * g)
            .sum::<f64>()
            .sqrt();
        if grad_norm < 1e-10 {
            break;
        }

        // Normalize gradient to unit direction; `alpha` controls
        // absolute step size in guide-CP units. Without normalization
        // step magnitude scales with gradient norm, which can blow
        // the candidate out of the fold-safe neighbourhood in a
        // single iteration.
        let inv_norm = 1.0 / grad_norm;
        let mut cand = bn_guides.clone();
        #[allow(clippy::needless_range_loop)]
        for k in 0..n {
            for slot in 0..64 {
                let (i, j) = cp_index(slot);
                let p = &mut cand[k].grid[i][j];
                p.x -= alpha * grad[k][slot][0] * inv_norm;
                p.y -= alpha * grad[k][slot][1] * inv_norm;
                p.z -= alpha * grad[k][slot][2] * inv_norm;
            }
        }
        let Ok(cand_binonic) = build_binonic_ring_from_binonic_guides(&cand, input) else {
            alpha *= 0.5;
            if alpha < ALPHA_MIN {
                break;
            }
            continue;
        };
        if any_self_intersects(&cand_binonic) {
            alpha *= 0.5;
            if alpha < ALPHA_MIN {
                break;
            }
            continue;
        }
        let cand_residual = measure_ring_residual(&cand_binonic);
        if cand_residual < best_residual {
            best_residual = cand_residual;
            best_binonic = cand_binonic;
            bn_guides = cand;
            alpha = (alpha * 1.25).min(ALPHA_INIT);
            accepted_count += 1;
        } else {
            alpha *= 0.5;
            if alpha < ALPHA_MIN {
                break;
            }
        }
    }

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!(
            "[kp_binonic N={n}] phase3 done: residual={best_residual:.3e} iters={iter_count} accepted={accepted_count}"
        );
    }

    // If phase 3 couldn't improve on phase 2 (fold-locked or
    // gradient vanished without progress), return whichever is
    // better.
    if best_residual < bq_result.final_residual {
        Some(KpBiNonicOptimizerResult {
            binonic: best_binonic,
            final_residual: best_residual,
        })
    } else {
        Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
        })
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

fn measure_ring_residual(patches: &[RectangularGregoryBiNonic]) -> f64 {
    let n = patches.len();
    let mut acc = 0.0_f64;
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
                let coplan = numer / denom;
                acc += coplan * coplan;
            }
        }
    }
    acc
}
