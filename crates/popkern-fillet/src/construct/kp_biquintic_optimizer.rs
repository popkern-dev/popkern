//! Phase H.8 (Option C step 6 + step 7) — bi-quintic guide
//! interior-CP optimizer.
//!
//! # Approach
//!
//! 1. Phase 1: FD gradient descent on the N central-vertex twists
//!    (3N DoF — small, converges fast; keeps existing GD).
//! 2. Phase 2: Degree-elevate each bicubic guide to bi-quintic
//!    (initial state is pointwise-identical). Then run
//!    Levenberg-Marquardt on the 16 interior bi-quintic CPs per
//!    sub-patch (48N DoF). The phase-3 L-M work validated that
//!    Gauss-Newton-damped steps recover from regressed starts where
//!    naive FD-GD stalls; the same reasoning applies here at a
//!    smaller DoF scale.
//! 3. Return the bi-nonic ring (reached via the bi-quintic guide
//!    path) plus the converged bi-quintic guides, so the caller's
//!    phase 3 can continue from the phase-2 optimum without
//!    re-solving.

use super::kp_binonic_optimizer::{analytic_ring_jacobian, flatten_ring};
use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
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
    /// Bi-quintic guides to use as phase-3's starting point when
    /// the caller escalates. Deliberately the *unoptimized*
    /// degree-elevated bicubic (not the phase-2 L-M output): phase
    /// 2's tighter local minimum in the 48-DoF sub-manifold can
    /// trap phase 3 at a worse G1 than starting neutral (observed
    /// on cube_extreme — phase 2 L-M residual 9.1e-2 made phase 3
    /// stall at G1 2.0e-1, while a 1.4e-1 neutral start let phase
    /// 3 reach G1 1.1e-2). The `binonic` + `final_residual` fields
    /// still reflect phase 2's optimized result — used when phase
    /// 3 is skipped.
    pub bq_guides: Vec<KpGuidePatchBiQuintic>,
}

/// Two-phase optimizer: first converges bicubic C-twists (cheap,
/// 3N DoF, FD-GD), then converges bi-quintic interior CPs (48N DoF,
/// Levenberg-Marquardt).
pub(super) fn optimize_biquintic_and_build_binonic(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Option<KpBiQuinticOptimizerResult> {
    const MAX_ITER_TWIST: usize = 10;
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
    // to bi-quintic, then run Levenberg-Marquardt on the 48N interior
    // CPs.
    let bicubic_guides = build_guide_ring_from_centrals(input, &central);
    let bq_guides_init: Vec<KpGuidePatchBiQuintic> = bicubic_guides
        .iter()
        .map(KpGuidePatchBiQuintic::from_bicubic_elevated)
        .collect();

    let initial = build_binonic_ring_from_biquintic_guides(&bq_guides_init, input).ok()?;
    if any_self_intersects(&initial) {
        // Initial bi-quintic build self-intersects — can't optimize
        // from a folded starting point. Return it anyway; the
        // downstream quality guard will revert.
        let r = measure_ring_residual(&initial);
        return Some(KpBiQuinticOptimizerResult {
            binonic: initial,
            final_residual: r,
            bq_guides: bq_guides_init,
        });
    }

    let x0 = pack_bq_guides(&bq_guides_init);
    let eval = BiQuinticEval {
        base: bq_guides_init.clone(),
        input,
        m_matrix: None,
    };
    // Phase-2 L-M can converge aggressively — its output is only
    // used when phase 3 is skipped (caller returns `binonic`
    // directly). When phase 3 engages, it starts from the neutral
    // `bq_guides_neutral` set below, so phase 2's local-minimum
    // choice doesn't leak into phase 3's basin.
    let lm_result = run_lm(
        eval,
        x0,
        &LmConfig {
            fd_eps: FD_EPS,
            lambda_init: 1e-3,
            lambda_max: 1e6,
            max_iter: 15,
            target_residual: TARGET_RESIDUAL,
        },
        "biquintic-phase2",
    );

    let mut best_bq_guides = bq_guides_init;
    unpack_bq_guides(&lm_result.x, &mut best_bq_guides);
    let best_binonic = match build_binonic_ring_from_biquintic_guides(&best_bq_guides, input) {
        Ok(ring) if !any_self_intersects(&ring) => ring,
        _ => {
            // L-M iterate landed in a folded or unbuildable region.
            // Fall back to the initial (degree-elevated) bi-quintic.
            best_bq_guides = bicubic_guides
                .iter()
                .map(KpGuidePatchBiQuintic::from_bicubic_elevated)
                .collect();
            build_binonic_ring_from_biquintic_guides(&best_bq_guides, input).ok()?
        }
    };
    let best_residual = measure_ring_residual(&best_binonic);

    if std::env::var_os("POPKERN_KP_TRACE").is_some() {
        eprintln!(
            "[kp_biquintic N={n}] phase2 (L-M bi-quintic interior) done: residual={best_residual:.3e}"
        );
    }

    // Return the phase-2 L-M optimized result as the shippable
    // bi-nonic ring, but hand phase 3 the *initial* (degree-
    // elevated bicubic) bq_guides so phase 3 can pick its own
    // basin — see the `bq_guides` field doc on
    // `KpBiQuinticOptimizerResult`.
    let bq_guides_neutral = bicubic_guides
        .iter()
        .map(KpGuidePatchBiQuintic::from_bicubic_elevated)
        .collect();

    Some(KpBiQuinticOptimizerResult {
        binonic: best_binonic,
        final_residual: best_residual,
        bq_guides: bq_guides_neutral,
    })
}

/// L-M residual callback for phase 2. Same analytic-Jacobian
/// pattern as `BinonicEval` in `kp_binonic_optimizer.rs`: M =
/// `∂ring_CP/∂guide_CP` is precomputed once (ring build is linear
/// in guide CPs), then `G = ∂f/∂ring_CP` is computed analytically
/// per iteration and `J = G · M` via the shared
/// `analytic_ring_jacobian` helper.
struct BiQuinticEval<'a> {
    base: Vec<KpGuidePatchBiQuintic>,
    input: &'a DirectBoundaryInput,
    /// Cached `D_ring × D_guide` linearization; see `BinonicEval`
    /// for the full argument. Shape here is `(N · 300) × (N · 48)`
    /// because bi-quintic guides have 16 interior CPs × 3 axes
    /// per sub-patch vs bi-nonic's 64 × 3.
    m_matrix: Option<Vec<Vec<f64>>>,
}

impl<'a> BiQuinticEval<'a> {
    fn guide_cp_count(&self) -> usize {
        self.base.len() * OPT_INDICES.len() * 3
    }

    fn build_ring(&self, x: &[f64]) -> Option<Vec<RectangularGregoryBiNonic>> {
        let mut guides = self.base.clone();
        unpack_bq_guides(x, &mut guides);
        build_binonic_ring_from_biquintic_guides(&guides, self.input).ok()
    }

    fn precompute_m_matrix(&mut self) -> Option<()> {
        let n = self.base.len();
        let d_ring = n * 300;
        let d_guide = self.guide_cp_count();

        let x0 = pack_bq_guides(&self.base);
        let baseline_ring = self.build_ring(&x0)?;
        let baseline_flat = flatten_ring(&baseline_ring);

        let eps = 1e-4;
        let mut m = vec![vec![0.0; d_guide]; d_ring];
        let mut x_pert = x0;
        for k in 0..d_guide {
            x_pert[k] += eps;
            let ring_pert = self.build_ring(&x_pert)?;
            let flat_pert = flatten_ring(&ring_pert);
            for i in 0..d_ring {
                m[i][k] = (flat_pert[i] - baseline_flat[i]) / eps;
            }
            x_pert[k] -= eps;
        }
        self.m_matrix = Some(m);
        Some(())
    }
}

impl<'a> LmEval for BiQuinticEval<'a> {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        let binonic = self.build_ring(x)?;
        Some(residual_vector(&binonic))
    }

    fn fold_safe(&mut self, x: &[f64]) -> bool {
        match self.build_ring(x) {
            Some(binonic) => !any_self_intersects(&binonic),
            None => false,
        }
    }

    fn jacobian(&mut self, x: &[f64]) -> Option<Vec<Vec<f64>>> {
        if self.m_matrix.is_none() {
            self.precompute_m_matrix()?;
        }
        let ring = self.build_ring(x)?;
        let d_guide = self.guide_cp_count();
        Some(analytic_ring_jacobian(
            &ring,
            self.m_matrix.as_ref()?,
            d_guide,
        ))
    }
}

fn pack_bq_guides(guides: &[KpGuidePatchBiQuintic]) -> Vec<f64> {
    let n = guides.len();
    let mut x = Vec::with_capacity(n * OPT_INDICES.len() * 3);
    for guide in guides {
        for &(i, j) in OPT_INDICES.iter() {
            x.push(guide.grid[i][j].x);
            x.push(guide.grid[i][j].y);
            x.push(guide.grid[i][j].z);
        }
    }
    x
}

fn unpack_bq_guides(x: &[f64], guides: &mut [KpGuidePatchBiQuintic]) {
    let mut idx = 0;
    for guide in guides.iter_mut() {
        for &(i, j) in OPT_INDICES.iter() {
            guide.grid[i][j].x = x[idx];
            guide.grid[i][j].y = x[idx + 1];
            guide.grid[i][j].z = x[idx + 2];
            idx += 3;
        }
    }
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

/// Residual *vector* — per-sample coplanarity. Matches the shape
/// used by `kp_binonic_optimizer::residual_vector` so both phases
/// share the same metric definition.
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
