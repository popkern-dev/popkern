//! Phase H.8 (Option C step 6) — bi-quintic guide interior-CP
//! optimizer.
//!
//! # Problem
//!
//! The existing `kp_twist_optimizer` converges to the bicubic
//! guide's quality ceiling — on cube_extreme (1:10 asymmetry) it
//! lands near 2.5e-1 worst-case G1 coplanarity, well above the
//! 5e-2 commercial gate. Confirmed diagnostic: the bicubic's 16
//! CPs per sub-patch cap the row-1 frame's cross-derivative
//! flexibility at degree-3 in v; asymmetric cross-curvature needs
//! degree-4+ representation along each boundary.
//!
//! # Approach
//!
//! 1. Run the existing bicubic twist optimizer to convergence.
//! 2. Degree-elevate each bicubic guide to bi-quintic (initial
//!    state is pointwise-identical — no regression baseline).
//! 3. FD gradient-descent on the **4 center interior CPs** of
//!    each sub-patch's bi-quintic guide (grid[2][2], [2][3],
//!    [3][2], [3][3]). 4 × 3 × N = 12N DoF — manageable for FD.
//! 4. Return the final binonic ring.
//!
//! The 4-CP inner block is a compromise between expressiveness
//! (all 16 non-boundary interior CPs) and convergence cost
//! (16 × 3 × N = 48N would be 4× more FD evals per iteration).
//! If it doesn't close cube_extreme, widening the DoF set is
//! the next natural step.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, GuideCornerData,
    KpGuidePatchBiQuintic, RectangularGregoryBiNonic, build_binonic_ring_from_biquintic_guides,
    build_binonic_ring_from_guides, build_guide_ring_from_centrals,
    extract_central_vertex_hermite_with_twist,
};

/// Interior-CP indices optimized: all 16 bi-quintic interior CPs
/// `grid[i][j]` for `i, j ∈ 1..=4`. These are the CPs that were
/// previously determined by degree-elevation and are now free
/// for the optimizer. Boundary-adjacent subsets change cross-
/// derivatives at shared spokes directly; deep-interior subsets
/// (e.g. [2][2]) mostly reshape the patch interior without
/// changing G1. Including all 16 gives the optimizer every
/// available DoF and lets the gradient pick the ones that
/// matter.
const OPT_INDICES: [(usize, usize); 16] = [
    (1, 1),
    (1, 2),
    (1, 3),
    (1, 4),
    (2, 1),
    (2, 2),
    (2, 3),
    (2, 4),
    (3, 1),
    (3, 2),
    (3, 3),
    (3, 4),
    (4, 1),
    (4, 2),
    (4, 3),
    (4, 4),
];

pub(super) struct KpBiQuinticOptimizerResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    #[allow(dead_code)]
    pub final_residual: f64,
    /// Phase-2 converged bi-quintic guides, so a caller can
    /// continue gradient descent at a higher degree without
    /// re-paying phases 1 and 2.
    pub bq_guides: Vec<KpGuidePatchBiQuintic>,
}

/// Two-phase optimizer: first converges bicubic C-twists (cheap,
/// ~9N DoF), then converges bi-quintic interior CPs (12N DoF).
pub(super) fn optimize_biquintic_and_build_binonic(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Option<KpBiQuinticOptimizerResult> {
    const MAX_ITER_TWIST: usize = 10;
    const MAX_ITER_BQ: usize = 20;
    const ALPHA_INIT: f64 = 1.0;
    const ALPHA_MIN: f64 = 1e-5;
    const FD_EPS: f64 = 1e-4;
    const TARGET_RESIDUAL: f64 = 1e-4;

    // Phase 1: bicubic C-twist optimization (mirrors kp_twist_optimizer).
    let mut central = extract_central_vertex_hermite_with_twist(input, options).ok()?;
    let n = central.len();
    if n == 0 {
        return None;
    }
    let mut best_bicubic_residual = bicubic_residual(input, &central).unwrap_or(f64::INFINITY);
    let mut alpha = ALPHA_INIT;
    for _iter in 0..MAX_ITER_TWIST {
        if best_bicubic_residual < TARGET_RESIDUAL {
            break;
        }
        let mut grad = vec![Point3::new(0.0, 0.0, 0.0); n];
        for k in 0..n {
            for axis in 0..3 {
                let mut perturbed = central.clone();
                match axis {
                    0 => perturbed[k].twist.x += FD_EPS,
                    1 => perturbed[k].twist.y += FD_EPS,
                    _ => perturbed[k].twist.z += FD_EPS,
                }
                if let Some(r) = bicubic_residual(input, &perturbed) {
                    let partial = (r - best_bicubic_residual) / FD_EPS;
                    match axis {
                        0 => grad[k].x = partial,
                        1 => grad[k].y = partial,
                        _ => grad[k].z = partial,
                    }
                }
            }
        }
        let mut cand = central.clone();
        for k in 0..n {
            cand[k].twist.x -= alpha * grad[k].x;
            cand[k].twist.y -= alpha * grad[k].y;
            cand[k].twist.z -= alpha * grad[k].z;
        }
        // Accept any improving step; SI filter is applied on the
        // final best-non-folded snapshot below.
        if let Some(r) = bicubic_residual(input, &cand)
            && r < best_bicubic_residual
        {
            best_bicubic_residual = r;
            central = cand;
            alpha = (alpha * 1.25).min(ALPHA_INIT);
            continue;
        }
        alpha *= 0.5;
        if alpha < ALPHA_MIN {
            break;
        }
    }

    // Phase 1 → phase 2 handoff: verify the converged bicubic is
    // non-folded. If not, fall back to the pre-optimisation twist.
    {
        let guides = build_guide_ring_from_centrals(input, &central);
        if let Ok(binonic) = build_binonic_ring_from_guides(&guides, input)
            && any_self_intersects(&binonic)
        {
            // Revert to initial Shirman-Séquin twist.
            central = extract_central_vertex_hermite_with_twist(input, options).ok()?;
        }
    }

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!(
            "[kp_biquintic N={n}] phase1 (bicubic twist) done: residual={best_bicubic_residual:.3e}"
        );
    }

    // Phase 2: build bicubic guides from converged twists, elevate
    // to bi-quintic, optimize interior CPs.
    let bicubic_guides = build_guide_ring_from_centrals(input, &central);
    let mut bq_guides: Vec<KpGuidePatchBiQuintic> = bicubic_guides
        .iter()
        .map(KpGuidePatchBiQuintic::from_bicubic_elevated)
        .collect();

    let initial = build_binonic_ring_from_biquintic_guides(&bq_guides, input).ok()?;
    if any_self_intersects(&initial) {
        // Initial bi-quintic build self-intersects — can't optimize
        // from a folded starting point. Return it anyway; the
        // downstream quality guard will revert.
        let r = measure_ring_residual(&initial);
        return Some(KpBiQuinticOptimizerResult {
            binonic: initial,
            final_residual: r,
            bq_guides,
        });
    }
    let mut best_residual = measure_ring_residual(&initial);
    let mut best_binonic = initial;
    let mut alpha = ALPHA_INIT;

    for _iter in 0..MAX_ITER_BQ {
        if best_residual < TARGET_RESIDUAL {
            break;
        }

        // FD gradient: N sub-patches × 16 CPs × 3 axes = 48N partials.
        let mut grad: Vec<[[f64; 3]; 16]> = vec![[[0.0; 3]; 16]; n];
        #[allow(clippy::needless_range_loop)]
        for k in 0..n {
            for (slot, &(i, j)) in OPT_INDICES.iter().enumerate() {
                for (axis, out) in grad[k][slot].iter_mut().enumerate() {
                    let mut bq_perturbed = bq_guides.clone();
                    let p = &mut bq_perturbed[k].grid[i][j];
                    match axis {
                        0 => p.x += FD_EPS,
                        1 => p.y += FD_EPS,
                        _ => p.z += FD_EPS,
                    }
                    let Ok(binonic) =
                        build_binonic_ring_from_biquintic_guides(&bq_perturbed, input)
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

        // Candidate step.
        let mut cand = bq_guides.clone();
        for k in 0..n {
            for (slot, &(i, j)) in OPT_INDICES.iter().enumerate() {
                let p = &mut cand[k].grid[i][j];
                p.x -= alpha * grad[k][slot][0];
                p.y -= alpha * grad[k][slot][1];
                p.z -= alpha * grad[k][slot][2];
            }
        }
        let Ok(cand_binonic) = build_binonic_ring_from_biquintic_guides(&cand, input) else {
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
            bq_guides = cand;
            alpha = (alpha * 1.25).min(ALPHA_INIT);
        } else {
            alpha *= 0.5;
            if alpha < ALPHA_MIN {
                break;
            }
        }
    }

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!(
            "[kp_biquintic N={n}] phase2 (bi-quintic interior) done: residual={best_residual:.3e}"
        );
    }

    Some(KpBiQuinticOptimizerResult {
        binonic: best_binonic,
        final_residual: best_residual,
        bq_guides,
    })
}

/// Helper for the bicubic twist phase: builds guides from centrals
/// and measures on the pure bicubic path. Matches the standalone
/// `kp_twist_optimizer` behavior exactly — no self-intersection
/// filter in the FD loop (that check happens in the final
/// accept/reject, not on every FD perturbation). Rejecting
/// perturbations mid-gradient-computation truncates the gradient
/// and stalls convergence around the initial seed.
fn bicubic_residual(input: &DirectBoundaryInput, central: &[GuideCornerData]) -> Option<f64> {
    let guides = build_guide_ring_from_centrals(input, central);
    let binonic = build_binonic_ring_from_guides(&guides, input).ok()?;
    Some(measure_ring_residual(&binonic))
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
