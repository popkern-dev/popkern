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
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, KpGuidePatchBiNonic,
    RectangularGregoryBiNonic, bernstein_binonic_all, bernstein_binonic_deriv_all,
    build_binonic_ring_from_binonic_guides,
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
        m_matrix: None,
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
    /// Cached `∂ring_CP / ∂guide_CP` linearization, shape `D_ring × D_guide`. The ring builder is linear in guide CPs (row-1/row-2 derivations + interior Bernstein-Gram LSQ + fixed boundary data), so this matrix is constant across optimization iterates — populated lazily on first `jacobian()` call and reused.
    m_matrix: Option<Vec<Vec<f64>>>,
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

/// Flat index into the ring-CP vector: patch `p`, cell `(i, j)`,
/// axis `a ∈ {0=x, 1=y, 2=z}`. `pub(super)` because the
/// analytic-Jacobian helper is shared with `kp_biquintic_optimizer`.
#[inline]
pub(super) fn ring_idx(p: usize, i: usize, j: usize, a: usize) -> usize {
    p * 300 + (i * 10 + j) * 3 + a
}

/// Flatten a bi-nonic ring's CPs to a row-major `D_ring` vector
/// `[patch_0_grid (row-major), patch_1_grid, ...]` with each CP
/// stored as `(x, y, z)`. Shared with `kp_biquintic_optimizer`
/// for the phase-2 analytic-Jacobian precompute.
pub(super) fn flatten_ring(ring: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let mut out = Vec::with_capacity(ring.len() * 300);
    for patch in ring {
        for i in 0..10 {
            for j in 0..10 {
                out.push(patch.grid[i][j].x);
                out.push(patch.grid[i][j].y);
                out.push(patch.grid[i][j].z);
            }
        }
    }
    out
}

/// Analytic Jacobian `J[sample, k] = ∂f[sample]/∂guide_CP_k`.
///
/// Both phase 2 (bi-quintic guides) and phase 3 (bi-nonic guides)
/// produce bi-nonic rings; this function takes the ring and the
/// pre-computed `M = ∂ring_CP/∂guide_CP` matrix and produces the
/// residual-vector Jacobian via the chain rule. It does NOT care
/// about the guide degree — all that matters is that `m.len() ==
/// ring.len() * 300` and each row of `m` has exactly `d_guide`
/// columns.
///
/// Cost: `O(samples · D_ring_relevant · D_guide)`. For the bi-nonic
/// output that's `50N · 30 · D_guide` per iteration.
pub(super) fn analytic_ring_jacobian(
    ring: &[RectangularGregoryBiNonic],
    m: &[Vec<f64>],
    d_guide: usize,
) -> Vec<Vec<f64>> {
    let n = ring.len();
    let total_samples = n * 50;
    let mut jac = vec![vec![0.0; d_guide]; total_samples];

    for spoke in 0..n {
        let a_idx = spoke;
        let b_idx = (spoke + 1) % n;
        let a = &ring[a_idx];
        let b = &ring[b_idx];

        for s in 0..50 {
            let t = 0.05 + 0.90 * (s as f64) / 49.0;
            let t_b = 1.0 - t;

            // Bernstein bases at this sample — constants the
            // chain-rule inner loop needs.
            let b_u = bernstein_binonic_all(t);
            let b_prime = bernstein_binonic_deriv_all(t);
            let b_v = bernstein_binonic_all(t_b);

            // Reference vectors at current iterate.
            let u = a.partial_u(1.0, t);
            let v = b.partial_v(t_b, 0.0);
            let w = a.partial_v(1.0, t);

            let cross_uv = u.cross(&v);
            let triple = cross_uv.x * w.x + cross_uv.y * w.y + cross_uv.z * w.z;
            let sign_t = if triple >= 0.0 { 1.0 } else { -1.0 };

            let nu_sq = u.x * u.x + u.y * u.y + u.z * u.z;
            let nv_sq = v.x * v.x + v.y * v.y + v.z * v.z;
            let nw_sq = w.x * w.x + w.y * w.y + w.z * w.z;
            let denom = nu_sq.sqrt() * nv_sq.sqrt() * nw_sq.sqrt();

            let row = spoke * 50 + s;
            if denom < 1e-300 {
                continue; // matches residual_vector's ill-conditioned zero
            }

            let f_val = triple.abs() / denom;

            // ∂f/∂u = sign_t · (v × w) / denom − f · u / |u|²  (similarly for v, w)
            let v_cross_w = v.cross(&w);
            let w_cross_u = w.cross(&u);
            let df_du = Vec3::new(
                sign_t * v_cross_w.x / denom - f_val * u.x / nu_sq,
                sign_t * v_cross_w.y / denom - f_val * u.y / nu_sq,
                sign_t * v_cross_w.z / denom - f_val * u.z / nu_sq,
            );
            let df_dv = Vec3::new(
                sign_t * w_cross_u.x / denom - f_val * v.x / nv_sq,
                sign_t * w_cross_u.y / denom - f_val * v.y / nv_sq,
                sign_t * w_cross_u.z / denom - f_val * v.z / nv_sq,
            );
            let df_dw = Vec3::new(
                sign_t * cross_uv.x / denom - f_val * w.x / nw_sq,
                sign_t * cross_uv.y / denom - f_val * w.y / nw_sq,
                sign_t * cross_uv.z / denom - f_val * w.z / nw_sq,
            );

            // Row via chain rule.
            // ∂u_α/∂g_k = 9 · Σ_j B_j(t) · (M[a,9,j,α; k] − M[a,8,j,α; k])
            // ∂v_α/∂g_k = 9 · Σ_i B_i(t_b) · (M[b,i,1,α; k] − M[b,i,0,α; k])
            // ∂w_α/∂g_k =      Σ_j B'_j(t) · M[a,9,j,α; k]
            let row_slice = &mut jac[row];
            for k in 0..d_guide {
                let mut du_x = 0.0_f64;
                let mut du_y = 0.0_f64;
                let mut du_z = 0.0_f64;
                let mut dv_x = 0.0_f64;
                let mut dv_y = 0.0_f64;
                let mut dv_z = 0.0_f64;
                let mut dw_x = 0.0_f64;
                let mut dw_y = 0.0_f64;
                let mut dw_z = 0.0_f64;
                for jj in 0..10 {
                    let bj = b_u[jj];
                    let r9x = m[ring_idx(a_idx, 9, jj, 0)][k];
                    let r9y = m[ring_idx(a_idx, 9, jj, 1)][k];
                    let r9z = m[ring_idx(a_idx, 9, jj, 2)][k];
                    let r8x = m[ring_idx(a_idx, 8, jj, 0)][k];
                    let r8y = m[ring_idx(a_idx, 8, jj, 1)][k];
                    let r8z = m[ring_idx(a_idx, 8, jj, 2)][k];
                    du_x += bj * (r9x - r8x);
                    du_y += bj * (r9y - r8y);
                    du_z += bj * (r9z - r8z);

                    let bp = b_prime[jj];
                    dw_x += bp * r9x;
                    dw_y += bp * r9y;
                    dw_z += bp * r9z;
                }
                for ii in 0..10 {
                    let bi = b_v[ii];
                    let rb1x = m[ring_idx(b_idx, ii, 1, 0)][k];
                    let rb1y = m[ring_idx(b_idx, ii, 1, 1)][k];
                    let rb1z = m[ring_idx(b_idx, ii, 1, 2)][k];
                    let rb0x = m[ring_idx(b_idx, ii, 0, 0)][k];
                    let rb0y = m[ring_idx(b_idx, ii, 0, 1)][k];
                    let rb0z = m[ring_idx(b_idx, ii, 0, 2)][k];
                    dv_x += bi * (rb1x - rb0x);
                    dv_y += bi * (rb1y - rb0y);
                    dv_z += bi * (rb1z - rb0z);
                }
                du_x *= 9.0;
                du_y *= 9.0;
                du_z *= 9.0;
                dv_x *= 9.0;
                dv_y *= 9.0;
                dv_z *= 9.0;

                row_slice[k] = df_du.x * du_x
                    + df_du.y * du_y
                    + df_du.z * du_z
                    + df_dv.x * dv_x
                    + df_dv.y * dv_y
                    + df_dv.z * dv_z
                    + df_dw.x * dw_x
                    + df_dw.y * dw_y
                    + df_dw.z * dw_z;
            }
        }
    }

    jac
}

impl<'a> LmEval for BinonicEval<'a> {
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
