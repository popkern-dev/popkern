//! Phase H.7 gap #2 — K-P 2004 §4 iterative guide-twist
//! optimizer for the K-P bidegree-9 construction.
//!
//! # Problem
//!
//! The one-shot K-P build (H.1–H.5) produces a ring of bicubic
//! guides whose central-vertex twists come from the Shirman-Séquin
//! Lagrangian QP. Under extreme radius asymmetry (e.g. cube_extreme
//! `r = [0.02, 0.10, 0.20]`), those twists are sub-optimal — the
//! downstream bidegree-9 patches show worst-case G1 coplanarity
//! ≈ 0.4 at `t ≈ 0.02`, concentrated at `C`. Dense-sampling
//! diagnostic confirmed all 3 spokes peak at the same location.
//!
//! # Iterative fix (K-P 2004 §4)
//!
//! Treat the N per-guide central-vertex twists as free parameters
//! (3N DoF in 3-space). The objective is the ring-wide worst-case
//! G1 coplanarity on shared spokes, with 50 samples per spoke on
//! `t ∈ [0.05, 0.95]`. Loop:
//!
//! 1. Build the binonic ring from the current twists.
//! 2. Measure residual (sum of coplanarities across all samples).
//! 3. Numerical FD gradient (3N perturbations).
//! 4. Gradient-descent step with backtracking line search.
//! 5. Repeat until residual converges or max iterations.
//!
//! The iterate starts at the Shirman-Séquin baseline, so the
//! optimization can only improve — a revert path keeps the best
//! seen iterate in case a descent step overshoots.
//!
//! # Scope
//!
//! MVP optimizes only the **central-vertex twists**. The outer-
//! corner twists are set by the K-P outer-corner QP and kept
//! fixed; the diagnostic localised the failure at `C`, not at
//! outer vertices. Extending to outer twists is a follow-up.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GuideCornerData, RectangularGregoryBiNonic,
    build_binonic_ring_from_guides, build_guide_ring_from_centrals,
    extract_central_vertex_hermite_with_twist,
};

/// Outcome of the iterative optimization. Includes the final
/// binonic ring and the residual achieved. The residual is
/// reported so the caller can decide whether to commit or revert
/// based on downstream quality gates.
pub(super) struct KpTwistOptimizerResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    /// Sum-of-squared-coplanarities at the returned iterate.
    /// Currently unused by the caller (the promotion ladder's
    /// quality guard re-measures in its own metric) but kept
    /// for diagnostics / future use.
    #[allow(dead_code)]
    pub final_residual: f64,
}

/// Iterative central-vertex-twist optimizer. Returns `None` if
/// the initial Shirman-Séquin solve or the initial binonic build
/// fails (same failure modes as `build_binonic_ring`).
///
/// Convergence heuristics:
///
/// - Max `MAX_ITER` iterations. Typical convergence: 3–5.
/// - Early exit when the residual drops below `TARGET_RESIDUAL`.
/// - Backtracking step size: start at `ALPHA_INIT`, scale `×0.5`
///   on a failed step, `×1.25` on a successful one. Bail when
///   `α < ALPHA_MIN` (no further progress possible).
/// - Best-iterate tracking: always return the lowest-residual
///   iterate seen, even if later iterations regress.
pub(super) fn optimize_c_twists_and_build_binonic(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Option<KpTwistOptimizerResult> {
    const MAX_ITER: usize = 10;
    const TARGET_RESIDUAL: f64 = 1e-4;
    const ALPHA_INIT: f64 = 1.0;
    const ALPHA_MIN: f64 = 1e-4;
    const FD_EPS: f64 = 1e-4;

    let initial_central = extract_central_vertex_hermite_with_twist(input, options).ok()?;
    let n = initial_central.len();
    if n == 0 {
        return None;
    }

    let mut central = initial_central.clone();
    let initial_binonic = build_initial_binonic(input, &central)?;
    let mut best_residual = measure_ring_residual(&initial_binonic);
    let mut best_binonic = initial_binonic;

    let mut alpha = ALPHA_INIT;

    for _iter in 0..MAX_ITER {
        if best_residual < TARGET_RESIDUAL {
            break;
        }

        // Numerical FD gradient w.r.t. each central twist component.
        let mut gradient = vec![Vec3::new(0.0, 0.0, 0.0); n];
        for k in 0..n {
            for axis in 0..3 {
                let mut perturbed = central.clone();
                match axis {
                    0 => perturbed[k].twist.x += FD_EPS,
                    1 => perturbed[k].twist.y += FD_EPS,
                    _ => perturbed[k].twist.z += FD_EPS,
                }
                let Some(binonic) = build_initial_binonic(input, &perturbed) else {
                    continue;
                };
                let r = measure_ring_residual(&binonic);
                let partial = (r - best_residual) / FD_EPS;
                match axis {
                    0 => gradient[k].x = partial,
                    1 => gradient[k].y = partial,
                    _ => gradient[k].z = partial,
                }
            }
        }

        // Gradient-norm scaling prevents a flat region's tiny
        // gradient from shrinking the step to nothing.
        let grad_norm_sq: f64 = gradient
            .iter()
            .map(|g| g.x * g.x + g.y * g.y + g.z * g.z)
            .sum();
        if grad_norm_sq < 1e-20 {
            break;
        }

        // Candidate step.
        let mut candidate = central.clone();
        for k in 0..n {
            candidate[k].twist.x -= alpha * gradient[k].x;
            candidate[k].twist.y -= alpha * gradient[k].y;
            candidate[k].twist.z -= alpha * gradient[k].z;
        }
        let Some(candidate_binonic) = build_initial_binonic(input, &candidate) else {
            alpha *= 0.5;
            if alpha < ALPHA_MIN {
                break;
            }
            continue;
        };
        let candidate_residual = measure_ring_residual(&candidate_binonic);

        if candidate_residual < best_residual {
            best_residual = candidate_residual;
            best_binonic = candidate_binonic;
            central = candidate;
            alpha = (alpha * 1.25).min(ALPHA_INIT);
        } else {
            alpha *= 0.5;
            if alpha < ALPHA_MIN {
                break;
            }
        }
    }

    Some(KpTwistOptimizerResult {
        binonic: best_binonic,
        final_residual: best_residual,
    })
}

/// Build a binonic ring given the current central-vertex Hermite
/// data. Returns `None` if the downstream build fails (should not
/// happen post-Shirman-Séquin success; the guide substrate is
/// robust under twist perturbation at the magnitudes we use).
fn build_initial_binonic(
    input: &DirectBoundaryInput,
    central: &[GuideCornerData],
) -> Option<Vec<RectangularGregoryBiNonic>> {
    let guides = build_guide_ring_from_centrals(input, central);
    build_binonic_ring_from_guides(&guides, input).ok()
}

/// Ring-wide sum-of-squared-coplanarities across all shared spokes.
/// 50 samples per spoke on `t ∈ [0.05, 0.95]`. Sum (not max)
/// because the optimizer needs a smooth scalar objective; the
/// max-over-samples is used by the commercial gate but is
/// non-differentiable and poor for gradient-based optimization.
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
