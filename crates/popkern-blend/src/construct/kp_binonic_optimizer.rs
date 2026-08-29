//! the bi-9 hardening pass (Option C step 3) — bi-nonic guide interior-CP
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

#[path = "kp_binonic_jacobian.rs"]
mod kp_binonic_jacobian;
use kp_binonic_jacobian::outer_fillet_jacobian;
pub(super) use kp_binonic_jacobian::{
    analytic_ring_jacobian, flatten_ring, outer_fillet_jacobian_sampled, ring_idx,
};

use super::kp_biquintic_optimizer::optimize_biquintic_and_build_binonic;
use super::kp_exact_boundary_binonic_optimizer::{
    optimize_exact_boundary_binonic, ring_crossing_decline_label,
};
use super::kp_lm_solver::{LmConfig, LmEval, run_lm};
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, KpGuidePatchBiNonic,
    RectangularGregoryBiNonic, build_binonic_ring_from_binonic_guides,
};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;

pub(super) struct KpBiNonicOptimizerResult {
    pub binonic: Vec<RectangularGregoryBiNonic>,
    #[allow(dead_code)]
    pub final_residual: f64,
    /// PK-275 rung observability: outcome of the phase-3 L-M when
    /// the outer-fillet block is armed (`None` when unarmed). The
    /// rung previously declined at six sites with no record; the
    /// promotion turns this into a structured
    /// `kp_binonic.promotion.outer_fillet_rung` fact.
    pub outer_rung_outcome: Option<String>,
}

/// Three-phase optimizer: delegates phases 1–2 to the bi-quintic
/// optimizer, then runs phase 3 (bi-nonic interior CPs via L-M)
/// if residual is still above target.
///
/// **`outer_fillets` (the outer-fillet-G1 work addition):** when `Some`, each
/// entry `outer_fillets[k]` is the surface the k-th cap sub-patch's
/// `u=0` boundary sits on (typically a cylinder fillet). The L-M refinement's
/// L-M residual is extended with per-sample cap-vs-fillet
/// G1 residuals along each cap's outer boundary, in addition to
/// the existing inner-spoke G1 residuals. Closes the user-visible
/// "messy top" defect (Task #144) when fed sphere-fit boundaries
/// that lie on the fillet surface.
///
/// When `None`, behaviour is identical to the pre-Phase-6c
/// optimizer (inner-spoke G1 only).
pub(super) fn optimize_binonic_with_outer_fillets(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
    outer_fillets: Option<&[SurfaceGeometry]>,
) -> Option<KpBiNonicOptimizerResult> {
    // Phase-2 residual above this triggers phase-3 L-M. Below it,
    // phase 2 is already in or near the commercial-gate basin (per-
    // sample coplanarity sum ≤ 1e-2 roughly corresponds to worst-
    // case G1 ≤ ~1e-2); skip the expensive L-M to keep common-case
    // test wall-time manageable.
    const PHASE3_ENGAGE_RESIDUAL: f64 = 1e-2;
    const LM_TARGET_RESIDUAL: f64 = 1e-6;

    let bq_result = optimize_biquintic_and_build_binonic(input, options)?;
    // When outer_fillets is set, ALWAYS engage phase 3
    // even if inner-spoke residual is already below threshold —
    // the outer-fillet residual is independent and may still be
    // large.
    let force_engage = outer_fillets.is_some();
    // Rung-outcome recorder: only populated when the outer block is
    // armed, so the unarmed path carries no new state.
    let rung = |outcome: &str| force_engage.then(|| outcome.to_string());
    if bq_result.final_residual < PHASE3_ENGAGE_RESIDUAL && !force_engage {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
            outer_rung_outcome: None,
        });
    }

    let n = bq_result.bq_guides.len();
    if n == 0 {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
            outer_rung_outcome: rung("declined:no_guides"),
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
            outer_rung_outcome: rung("declined:initial_ring_self_intersects"),
        });
    }
    if input.hex_boundaries.is_some()
        && outer_fillets.is_none()
        && let Some(direct) = optimize_exact_boundary_binonic(&initial)
        && direct.final_residual < bq_result.final_residual
    {
        return Some(KpBiNonicOptimizerResult {
            binonic: direct.binonic,
            final_residual: direct.final_residual,
            outer_rung_outcome: None,
        });
    }

    // Pack the 64 interior CPs × 3 axes × N patches into a flat
    // parameter vector. Order: sub-patch k, slot s ∈ 0..64, axis a.
    let x0 = pack_guides(&bn_guides);
    if x0.len() > MAX_LM_REFINE_DOF {
        return Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_result.final_residual,
            outer_rung_outcome: rung("declined:dof_cap"),
        });
    }

    let eval = BinonicEval {
        base: bn_guides.clone(),
        input,
        m_matrix: None,
        outer_fillets: outer_fillets.map(|f| f.to_vec()),
    };
    // Outer-fillet-G1 makes the L-M problem harder
    // (more constraints), so bump max_iter substantially when the
    // outer block is active. Also tighter target since we want
    // commercial-grade G1.
    let (max_iter, target_residual) = if outer_fillets.is_some() {
        (60, 1e-9)
    } else {
        (15, LM_TARGET_RESIDUAL)
    };
    let lm_result = run_lm(
        eval,
        x0,
        &LmConfig {
            fd_eps: 1e-4,
            lambda_init: 1e-3,
            lambda_max: 1e6,
            max_iter,
            target_residual,
            require_feasible_start: true,
        },
        "binonic-phase3",
    );

    // Unpack optimal parameter back into guides; rebuild ring.
    let lm_stats = format!("iters={};accepted={}", lm_result.iters, lm_result.accepted);
    let mut best_guides = bn_guides;
    unpack_guides(&lm_result.x, &mut best_guides);
    let best_binonic = match build_binonic_ring_from_binonic_guides(&best_guides, input) {
        // PK-283 ring-level crossing honesty: the per-patch arm cannot
        // see two ring patches crossing each other. Reject the LM
        // product and keep the phase-2 ring (whose own honesty the
        // promotion-level gate enforces).
        Ok(ring) if !any_self_intersects(&ring) => match ring_crossing_decline_label(&ring) {
            None => ring,
            Some(label) => {
                return Some(KpBiNonicOptimizerResult {
                    binonic: bq_result.binonic,
                    final_residual: bq_result.final_residual,
                    outer_rung_outcome: rung(&format!(
                        "lm_rejected:ring_crossing;{label};{lm_stats}"
                    )),
                });
            }
        },
        _ => {
            return Some(KpBiNonicOptimizerResult {
                binonic: bq_result.binonic,
                final_residual: bq_result.final_residual,
                outer_rung_outcome: rung(&format!("lm_rejected:rebuild_or_fold;{lm_stats}")),
            });
        }
    };

    // Guard against regression: if L-M ended worse than phase 2
    // on the metric we're optimising for, return phase-2 result.
    // When outer_fillets is set, the optimisation
    // metric is inner-spoke + outer-fillet residual squared sum.
    // Check on the same combined metric.
    let post_inner = measure_ring_residual(&best_binonic);
    let bq_post_inner = bq_result.final_residual;
    let (post, bq_post) = if let Some(fillets) = outer_fillets {
        let post_outer: f64 = outer_fillet_residual_vector(&best_binonic, fillets)
            .iter()
            .map(|v| v * v)
            .sum();
        let bq_post_outer: f64 = outer_fillet_residual_vector(&bq_result.binonic, fillets)
            .iter()
            .map(|v| v * v)
            .sum();
        (post_inner + post_outer, bq_post_inner + bq_post_outer)
    } else {
        (post_inner, bq_post_inner)
    };
    if post < bq_post {
        Some(KpBiNonicOptimizerResult {
            binonic: best_binonic,
            final_residual: post,
            outer_rung_outcome: rung(&format!(
                "served_lm;{lm_stats};post={post:.3e};presolve={bq_post:.3e}"
            )),
        })
    } else {
        Some(KpBiNonicOptimizerResult {
            binonic: bq_result.binonic,
            final_residual: bq_post,
            outer_rung_outcome: rung(&format!(
                "kept_presolve;{lm_stats};post={post:.3e};presolve={bq_post:.3e}"
            )),
        })
    }
}

const MAX_LM_REFINE_DOF: usize = 1_000;

/// L-M residual callback for the refinement stage. Keeps a base copy of the
/// guide ring and a reference to the input; each `eval` clones
/// the base, overwrites interior CPs from the parameter vector,
/// rebuilds the binonic ring, and returns a per-sample coplanarity
/// vector.
struct BinonicEval<'a> {
    base: Vec<KpGuidePatchBiNonic>,
    input: &'a DirectBoundaryInput,
    /// Cached `∂ring_CP / ∂guide_CP` linearization, shape `D_ring × D_guide`. The ring builder is linear in guide CPs (row-1/row-2 derivations + interior Bernstein-Gram LSQ + fixed boundary data), so this matrix is constant across optimization iterates — populated lazily on first `jacobian()` call and reused.
    m_matrix: Option<Vec<Vec<f64>>>,
    /// the outer-fillet-G1 work addition: when `Some`, each entry is the
    /// surface that sub-patch k's `u=0` boundary sits on. The
    /// residual vector is extended with per-sample cap-vs-fillet
    /// G1 residuals along each cap's outer boundary.
    outer_fillets: Option<Vec<SurfaceGeometry>>,
}

impl<'a> BinonicEval<'a> {
    fn guide_cp_count(&self) -> usize {
        self.base.len() * 192
    }

    fn build_ring(&self, x: &[f64]) -> Option<Vec<RectangularGregoryBiNonic>> {
        let mut guides = self.base.clone();
        unpack_guides(x, &mut guides);
        build_binonic_ring_from_binonic_guides(&guides, self.input).ok()
    }

    /// One-time FD pre-computation of `∂ring_CP/∂guide_CP`. Since
    /// the ring build is linear in guide CPs, this matrix is
    /// constant — one FD pass over the guide DoF gives the exact
    /// linear operator that we'll reuse forever.
    fn precompute_m_matrix(&mut self) -> Option<()> {
        let n = self.base.len();
        let d_ring = n * 300;
        let d_guide = n * 192;

        let x0 = pack_guides(&self.base);
        let baseline_ring = self.build_ring(&x0)?;
        let baseline_flat = flatten_ring(&baseline_ring);

        // FD step scaled for numerical stability. The ring builder's
        // linear operations mean any nonzero eps produces the same
        // gradient up to roundoff; 1e-4 trades truncation against
        // subtraction cancellation.
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

impl<'a> LmEval for BinonicEval<'a> {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        let binonic = self.build_ring(x)?;
        let mut r = residual_vector(&binonic);
        if let Some(fillets) = &self.outer_fillets {
            r.extend(outer_fillet_residual_vector(&binonic, fillets));
        }
        Some(r)
    }

    fn fold_safe(&mut self, x: &[f64]) -> bool {
        match self.build_ring(x) {
            Some(binonic) => !any_self_intersects(&binonic),
            None => false,
        }
    }

    fn jacobian(&mut self, x: &[f64]) -> Option<Vec<Vec<f64>>> {
        // PK-275: both residual blocks have analytic Jacobians now.
        // The inner-spoke block chains ∂f/∂(cap tangents) through the
        // constant ring linearization M; the outer-fillet block adds
        // the ∂f/∂(u=0 boundary) chain including footpoint motion on
        // the fillet. Before this, `outer_fillets.is_some()` returned
        // `None` here, forcing the LM driver onto full-DoF FD — one
        // complete ring rebuild per DoF per iteration (576 rebuilds ×
        // iteration at N = 3).
        if self.m_matrix.is_none() {
            self.precompute_m_matrix()?;
        }
        let ring = self.build_ring(x)?;
        let d_guide = self.guide_cp_count();
        let m = self.m_matrix.as_ref()?;
        let mut jac = analytic_ring_jacobian(&ring, m, d_guide);
        if let Some(fillets) = &self.outer_fillets {
            jac.extend(outer_fillet_jacobian(&ring, fillets, m, d_guide));
        }
        Some(jac)
    }
}

/// Per-sample cap-vs-fillet G1 residual along each
/// sub-patch's `u=0` boundary. Sample count: 50 per sub-patch.
///
/// G1 condition: at every boundary point, the cap and fillet
/// share a tangent plane (their normals are parallel).
///
/// Residual = `|n_cap × n_fillet| / (|n_cap| · |n_fillet|)`,
/// the sine of the angle between the two normals. Range [0, 1];
/// 0 = perfect tangent-plane match, 1 = orthogonal normals.
fn outer_fillet_residual_vector(
    patches: &[RectangularGregoryBiNonic],
    outer_fillets: &[SurfaceGeometry],
) -> Vec<f64> {
    let n = patches.len();
    let mut out = Vec::with_capacity(n * 50);
    for (k, cap) in patches.iter().enumerate() {
        let Some(fillet) = outer_fillets.get(k) else {
            // Pad with zeros so the residual-vector length stays
            // predictable across `n` and across iterations.
            out.extend(std::iter::repeat_n(0.0, 50));
            continue;
        };
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            // Cap's tangent plane at u=0, v=t:
            //   spanned by ∂cap/∂u and ∂cap/∂v.
            let su_cap = cap.partial_u(0.0, t);
            let sv_cap = cap.partial_v(0.0, t);
            let n_cap = su_cap.cross(&sv_cap);
            // Cap surface point at this UV — also the 3-D location
            // we project onto the fillet to get the fillet normal.
            let p_3d = cap.point_at(0.0, t);
            let Some((u_f, v_f)) = project_point_to_uv(fillet, &p_3d) else {
                out.push(0.0);
                continue;
            };
            let su_f = fillet.partial_u(u_f, v_f);
            let sv_f = fillet.partial_v(u_f, v_f);
            let n_f = su_f.cross(&sv_f);
            let cross = n_cap.cross(&n_f);
            let nc_len = (n_cap.x.powi(2) + n_cap.y.powi(2) + n_cap.z.powi(2)).sqrt();
            let nf_len = (n_f.x.powi(2) + n_f.y.powi(2) + n_f.z.powi(2)).sqrt();
            let cross_len = (cross.x.powi(2) + cross.y.powi(2) + cross.z.powi(2)).sqrt();
            let denom = nc_len * nf_len;
            if denom > 1e-300 {
                out.push(cross_len / denom);
            } else {
                out.push(0.0);
            }
        }
    }
    out
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

#[cfg(test)]
#[path = "kp_binonic_jacobian_fd_tests.rs"]
mod kp_binonic_jacobian_fd_tests;
