//! Phase-10 C — N-spoke ring coupler on top of the single-spoke
//! polynomial-identity machinery in [`super::spoke_poly_g1_solver`].
//!
//! # What this module does
//!
//! Given `N` bidegree-7 two-stripe Gregory patches arranged in a
//! ring (the sub-patches of a three-quad decomposition) and `N`
//! shared spokes, build **one** block-structured linear system
//! whose unknowns are:
//!
//! * every row-1 near-spoke CP of every patch (deduplicated —
//!   a CP shared between adjacent spokes appears once), and
//! * two scalar `β` Bernstein coefficients per spoke.
//!
//! The constraint rows are per-spoke Bernstein-coefficient
//! equations as built in [`super::spoke_poly_g1_solver`], rewritten
//! against the global unknown vector.
//!
//! # Why the coupling matters
//!
//! In the three-quad ring, patch `k`'s row-1 near-spoke CP at
//! grid `(6, 1)` (the central-vertex-adjacent CP) is used by
//! **both** spokes incident on patch `k` (incoming on VZero,
//! outgoing on UOne). Solving each spoke in isolation and
//! writing back would double-write this CP with inconsistent
//! values — destroying the global structure. The ring solver
//! assigns it a single scalar variable index in the unknown
//! vector and lets both spokes' constraints reference it.
//!
//! # The α sign
//!
//! The K-P identity takes `xtan_A + α · xtan_B + β · T ≡ 0`. The
//! sign of `α` depends on whether the two patches' cross-tangents
//! at their spoke boundaries point in the same or opposite
//! geometric direction:
//!
//! * Consistent-orientation (both partials point outward into
//!   their interiors from the same "side" of the spoke):
//!   `α = −1` and the identity is `xtan_A − xtan_B = β · T`.
//! * Opposite-orientation: `α = +1` and the identity is
//!   `xtan_A + xtan_B = β · T`.
//!
//! We auto-detect per-spoke at assembly time by inspecting the
//! initial cross-tangents. The rest of the machinery is linear
//! and degree-escalation-friendly regardless of which sign the
//! geometry requires.

#![allow(clippy::needless_range_loop)]
// Fillet-crate wiring lands later. Dead-code warnings
// here are expected until that integration lands.
#![allow(dead_code)]

use crate::geometry::bernstein_poly::{BernsteinPoly, bernstein_product_scalar_vec3};
use crate::geometry::gregory::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;
use crate::geometry::gregory::spoke_g1_solver::Spoke;
use crate::geometry::gregory::spoke_poly_g1_solver::{
    boundary_bernstein_cps, boundary_tangent_poly, free_row1_indices_for_side,
    row1_indices_for_side, xtan_bernstein_cps,
};
use crate::geometry::{Point3, Vec3};
use std::collections::BTreeMap;

/// Outcome of a ring solve. Residual and β are reported per spoke
/// so callers can gate promotion and drive degree-escalation when
/// a particular spoke's system is inconsistent.
#[derive(Debug, Clone)]
pub struct RingPolyG1Result {
    /// Post-solve `‖LHS‖₂` maxed over each spoke's 8 Bernstein
    /// coefficient equations. A healthy solve produces values at
    /// machine precision.
    pub per_spoke_residual_max: Vec<f64>,
    /// Per-spoke `(β₀, β₁)` after the solve.
    pub per_spoke_beta: Vec<(f64, f64)>,
    /// Per-spoke `α` sign chosen from geometry. `−1` or `+1`.
    pub per_spoke_alpha: Vec<f64>,
}

/// Solve the polynomial-identity G1 system over a ring of N
/// two-stripe patches and N spokes. Mutates `patches` in place
/// on success. Returns `None` only on a singular KKT matrix —
/// the caller should then fall back to the sampled-damped
/// iteration path (or escalate `β` degree in future solver work).
pub fn solve_ring_poly_g1(
    patches: &mut [RectangularGregoryTwoStripe],
    spokes: &[Spoke],
) -> Option<RingPolyG1Result> {
    if patches.is_empty() || spokes.is_empty() {
        return Some(RingPolyG1Result {
            per_spoke_residual_max: Vec::new(),
            per_spoke_beta: Vec::new(),
            per_spoke_alpha: Vec::new(),
        });
    }

    // 1. Build variable-index map over the GENUINELY FREE row-1
    //    near-spoke CPs (`j ∈ 1..7`). The `j = 0` / `j = 7` row-1
    //    positions lie on the patch's OTHER (hex-edge) boundaries
    //    and MUST NOT move — moving them creates a G0 position
    //    gap with the adjacent sub-patch across that hex edge.
    //    Shared CPs (the corner-adjacent position where two
    //    spokes incident on the same patch overlap) get one
    //    variable index, not two.
    let mut var_offset: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut next_var = 0_usize;
    for spoke in spokes {
        for end in [spoke.a, spoke.b] {
            if end.patch_idx >= patches.len() {
                return None;
            }
            for ij in free_row1_indices_for_side(end.side) {
                let key = (end.patch_idx, ij.0, ij.1);
                if let std::collections::btree_map::Entry::Vacant(e) = var_offset.entry(key) {
                    e.insert(next_var);
                    next_var += 3;
                }
            }
        }
    }
    let beta_var_start = next_var;
    let n_unknowns = beta_var_start + 2 * spokes.len();

    // 2. Collect the seed vector from current patch state.
    let mut seed = vec![0.0_f64; n_unknowns];
    for (&(patch_idx, i, j), &off) in &var_offset {
        let p = patches[patch_idx].grid[i][j];
        seed[off] = p.x;
        seed[off + 1] = p.y;
        seed[off + 2] = p.z;
    }
    // β seeds at zero (identity flat).

    // 3. Assemble 24·N constraint rows.
    let n_spokes = spokes.len();
    let mut a_mat: Vec<Vec<f64>> = Vec::with_capacity(24 * n_spokes);
    let mut rhs: Vec<f64> = Vec::with_capacity(24 * n_spokes);
    let mut per_spoke_alpha = Vec::with_capacity(n_spokes);

    // Reusable per-spoke linear β basis: β = B_0^1 and β = B_1^1.
    let beta0_poly = BernsteinPoly::<f64>::new(vec![1.0, 0.0]);
    let beta1_poly = BernsteinPoly::<f64>::new(vec![0.0, 1.0]);

    for (s_idx, spoke) in spokes.iter().enumerate() {
        // xtan initial values for α sign detection.
        let xtan_a_init = xtan_bernstein_cps(&patches[spoke.a.patch_idx], spoke.a.side);
        let xtan_b_init_raw = xtan_bernstein_cps(&patches[spoke.b.patch_idx], spoke.b.side);
        let xtan_b_init: [Vec3; 8] = if spoke.b.reversed {
            let mut out = [Vec3::zero(); 8];
            for k in 0..8 {
                out[k] = xtan_b_init_raw[7 - k];
            }
            out
        } else {
            xtan_b_init_raw
        };
        let dot_mid = xtan_a_init[3].x * xtan_b_init[3].x
            + xtan_a_init[3].y * xtan_b_init[3].y
            + xtan_a_init[3].z * xtan_b_init[3].z;
        let alpha = if dot_mid >= 0.0 { -1.0 } else { 1.0 };
        per_spoke_alpha.push(alpha);

        // Boundary tangent T(t) for this spoke.
        let boundary_cps = boundary_bernstein_cps(&patches[spoke.a.patch_idx], spoke.a.side);
        let t_poly = boundary_tangent_poly(&boundary_cps);
        let beta0_t = bernstein_product_scalar_vec3(&beta0_poly, &t_poly);
        let beta1_t = bernstein_product_scalar_vec3(&beta1_poly, &t_poly);

        let a_row1_indices_all = row1_indices_for_side(spoke.a.side);
        let b_row1_indices_raw_all = row1_indices_for_side(spoke.b.side);

        let beta0_global_idx = beta_var_start + 2 * s_idx;
        let beta1_global_idx = beta_var_start + 2 * s_idx + 1;

        // the abandoned polynomial-identity direct-formula attempt+ — all 8 Bernstein coefficient equations.
        // At k ∈ {1..6}: row-1 CPs on BOTH ends are free, so
        // they appear in the constraint rows with `+7` and `α·7`.
        // At k ∈ {0, 7}: the row-1 CPs are pinned (they sit on
        // adjacent-spoke boundaries / hex-external edges). Their
        // values contribute to the RHS; the only unknowns in
        // these rows are `β_0`, `β_1`. The constraint then says
        // `β · T[k] = (known constant)`, which β absorbs if the
        // seed geometry is G1-at-endpoint consistent (Shirman-
        // Séquin-valid at t=0, hex-boundary-G1-valid at t=1).
        //
        // If the endpoint G1 condition happens to be INCONSISTENT
        // with any linear β (3 scalar equations per pinned
        // endpoint, 2 β DOFs), the KKT Tikhonov produces a
        // least-squares approximation, and the caller's
        // `per_spoke_residual_max` gate detects it and falls
        // back to the sampled solver.
        for k in 0..8 {
            // Grid indices on A, B for Bernstein coefficient k.
            let a_ij = a_row1_indices_all[k];
            let k_b = if spoke.b.reversed { 7 - k } else { k };
            let b_ij = b_row1_indices_raw_all[k_b];

            // Look up solver offsets — present iff the CP is free.
            let a_off = var_offset.get(&(spoke.a.patch_idx, a_ij.0, a_ij.1));
            let b_off = var_offset.get(&(spoke.b.patch_idx, b_ij.0, b_ij.1));

            // Pinned CP values (used when a_off / b_off is None,
            // or to compute the seed contribution for free CPs).
            let a_pinned_val = patches[spoke.a.patch_idx].grid[a_ij.0][a_ij.1];
            let b_pinned_val = patches[spoke.b.patch_idx].grid[b_ij.0][b_ij.1];

            let w0 = beta0_t.coeff(k);
            let w1 = beta1_t.coeff(k);
            let boundary_k = boundary_cps[k];
            let rhs_factor = (1.0 + alpha) * 7.0;

            for c in 0..3 {
                let mut row = vec![0.0_f64; n_unknowns];
                let (w0_c, w1_c, c_c, a_p_c, b_p_c) = match c {
                    0 => (w0.x, w1.x, boundary_k.x, a_pinned_val.x, b_pinned_val.x),
                    1 => (w0.y, w1.y, boundary_k.y, a_pinned_val.y, b_pinned_val.y),
                    _ => (w0.z, w1.z, boundary_k.z, a_pinned_val.z, b_pinned_val.z),
                };
                let mut rhs_c = rhs_factor * c_c;

                if let Some(off) = a_off {
                    row[off + c] = 7.0;
                } else {
                    // A's row-1 CP at this k is pinned; move its
                    // contribution (`7·A_pinned`) to the RHS.
                    rhs_c -= 7.0 * a_p_c;
                }
                if let Some(off) = b_off {
                    row[off + c] = alpha * 7.0;
                } else {
                    rhs_c -= alpha * 7.0 * b_p_c;
                }
                row[beta0_global_idx] = w0_c;
                row[beta1_global_idx] = w1_c;
                a_mat.push(row);
                rhs.push(rhs_c);
            }
        }
    }

    // 4. Min-norm KKT: x = seed + Aᵀ · λ, (A · Aᵀ + λI) · λ = rhs − A · seed.
    let m = a_mat.len();
    let mut a_seed = vec![0.0_f64; m];
    for r in 0..m {
        for k in 0..n_unknowns {
            a_seed[r] += a_mat[r][k] * seed[k];
        }
    }
    let rhs_adj: Vec<f64> = (0..m).map(|r| rhs[r] - a_seed[r]).collect();

    let mut aat = vec![vec![0.0_f64; m]; m];
    let mut trace = 0.0_f64;
    for r in 0..m {
        for c in r..m {
            let mut s = 0.0_f64;
            for k in 0..n_unknowns {
                s += a_mat[r][k] * a_mat[c][k];
            }
            aat[r][c] = s;
            aat[c][r] = s;
        }
        trace += aat[r][r];
    }
    let lambda_reg = 1e-12 * (trace / (m as f64)).max(1.0);
    for r in 0..m {
        aat[r][r] += lambda_reg;
    }
    let lambdas = solve_linear_system(aat, rhs_adj)?;

    let mut x = seed.clone();
    for r in 0..m {
        for k in 0..n_unknowns {
            x[k] += lambdas[r] * a_mat[r][k];
        }
    }

    // 5. Write back CP values.
    for (&(patch_idx, i, j), &off) in &var_offset {
        patches[patch_idx].grid[i][j] = Point3::new(x[off], x[off + 1], x[off + 2]);
    }

    // 6. Post-solve residual per spoke. 8 Bernstein equations × 3
    //    scalar coords = 24 rows per spoke (k ∈ 0..8).
    let rows_per_spoke = 24;
    let mut per_spoke_residual_max = Vec::with_capacity(n_spokes);
    let mut per_spoke_beta = Vec::with_capacity(n_spokes);
    for s_idx in 0..n_spokes {
        let mut r_max = 0.0_f64;
        for k_local in 0..8 {
            let row_base = rows_per_spoke * s_idx + 3 * k_local;
            let mut lhs = [0.0_f64; 3];
            for c in 0..3 {
                let row = row_base + c;
                let mut acc = -rhs[row];
                for (j, &v) in a_mat[row].iter().enumerate() {
                    acc += v * x[j];
                }
                lhs[c] = acc;
            }
            let mag = (lhs[0].powi(2) + lhs[1].powi(2) + lhs[2].powi(2)).sqrt();
            if mag > r_max {
                r_max = mag;
            }
        }
        per_spoke_residual_max.push(r_max);
        let b0_idx = beta_var_start + 2 * s_idx;
        per_spoke_beta.push((x[b0_idx], x[b0_idx + 1]));
    }

    Some(RingPolyG1Result {
        per_spoke_residual_max,
        per_spoke_beta,
        per_spoke_alpha,
    })
}

fn solve_linear_system(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = m.len();
    for i in 0..n {
        let mut pivot = i;
        let mut best = m[i][i].abs();
        for (r, row) in m.iter().enumerate().skip(i + 1) {
            let v = row[i].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if pivot != i {
            m.swap(i, pivot);
            b.swap(i, pivot);
        }
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            let pivot_row = m[i].clone();
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                m[r][c] -= factor * pc;
            }
            b[r] -= factor * b[i];
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for (c, xc) in x.iter().enumerate().skip(i + 1) {
            s -= m[i][c] * xc;
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

#[cfg(test)]
#[path = "spoke_poly_g1_ring_tests.rs"]
mod tests;
