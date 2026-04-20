//! N-spoke ring-level row-1 G1 solver.
//!
//! Extends [`super::direct_row1_g1::solve_single_spoke_row1_g1`]
//! from one spoke to the full `N`-spoke ring of a three-quad
//! decomposition. The key addition is **global CP coupling**:
//! the corner-C-adjacent CP `grid_k[6][1]` on each patch `k`
//! appears in TWO spokes' row-1 unknown sets — spoke `k`'s A
//! side at `j=1` and spoke `(k-1)`'s B side at `k'=1`. Solving
//! the spokes independently would over-write that CP with
//! inconsistent values; the ring solver registers it as a
//! single variable index so both spokes' constraints converge
//! on one value.
//!
//! Pattern mirrors the Phase C polynomial-identity ring solver
//! in [`super::spoke_poly_g1_ring`]: BTreeMap-based variable
//! indexing for determinism, stacked constraint matrix,
//! min-norm KKT solve.

#![allow(clippy::needless_range_loop)]

use crate::geometry::bernstein_poly::{BernsteinPoly, bernstein_product_scalar_vec3};
use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::gregory::direct_blend_coeffs::SpokeBlendCoeffs;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::gregory::direct_linalg::{binomial, solve_dense};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::gregory::direct_types::CornerKind;
use crate::geometry::{Point3, Vec3};
use std::collections::BTreeMap;

/// Output of [`solve_ring_row1_g1`].
#[derive(Debug, Clone)]
pub struct RingRow1Result {
    /// The updated patches with row-1 CPs refined to satisfy the
    /// polynomial-identity G1 condition.
    pub updated_patches: Vec<RectangularGregoryTwoStripe>,
    /// Max `‖LHS coefficient‖₂` across every spoke's 9 Bernstein
    /// equations. Machine precision after a successful solve.
    pub max_residual: f64,
    /// Per-spoke residual maxes, for UI / telemetry.
    pub per_spoke_residual: Vec<f64>,
}

/// Solve the row-1 G1 identity across the full N-spoke ring.
///
/// Requires `input.seed_patches` to be `Some(_)` — we need the
/// full 8×8 grids as the seed for the min-norm solve and as the
/// source of pinned boundary CPs at `j=0` and `j=7` on each
/// side of each spoke.
///
/// `blend_coeffs.len()` must equal `input.spokes.len()`.
pub fn solve_ring_row1_g1(
    input: &DirectBoundaryInput,
    blend_coeffs: &[SpokeBlendCoeffs],
    options: &DirectConstructionOptions,
) -> Result<RingRow1Result, DirectConstructionError> {
    let Some(seed) = input.seed_patches.as_ref() else {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: "ring row-1 solve requires seed_patches (full biseptic-promoted grids)".into(),
        });
    };
    let n = input.n();
    if seed.len() != n {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: format!("seed_patches length {} ≠ ring size {n}", seed.len()),
        });
    }
    if blend_coeffs.len() != n {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: format!("blend_coeffs length {} ≠ ring size {n}", blend_coeffs.len()),
        });
    }

    // 1. Enumerate unknowns via BTreeMap (patch_idx, i, j) →
    //    starting variable offset. Each CP takes 3 scalar
    //    variables (x, y, z).
    let mut var_offset: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut next_var: usize = 0;
    for k in 0..n {
        // Spoke k's A side: patch k, grid[6][j] for j=1..=6.
        for j in 1..=6 {
            var_offset.entry((k, 6, j)).or_insert_with(|| {
                let v = next_var;
                next_var += 3;
                v
            });
        }
        // Spoke k's B side: patch (k+1), grid[7-k'][1] for k'=1..=6.
        // i.e. grid[1][1]..grid[6][1] on patch (k+1).
        let patch_b = (k + 1) % n;
        for kp in 1..=6 {
            let i = 7 - kp;
            var_offset.entry((patch_b, i, 1)).or_insert_with(|| {
                let v = next_var;
                next_var += 3;
                v
            });
        }
    }
    let n_unknowns = next_var;

    // 2. Seed vector from current patch state.
    let mut seed_vec = vec![0.0_f64; n_unknowns];
    for (&(patch_idx, i, j), &off) in &var_offset {
        let p = seed[patch_idx].grid[i][j];
        seed_vec[off] = p.x;
        seed_vec[off + 1] = p.y;
        seed_vec[off + 2] = p.z;
    }

    // 3. Assemble the global constraint matrix. Each spoke
    //    contributes 9 Bernstein vector equations × 3 coords =
    //    27 scalar rows.
    let mut a_mat: Vec<Vec<f64>> = Vec::with_capacity(27 * n);
    let mut rhs: Vec<f64> = Vec::with_capacity(27 * n);
    let beta0_basis = BernsteinPoly::<f64>::new(vec![1.0, 0.0]);
    let beta1_basis = BernsteinPoly::<f64>::new(vec![0.0, 1.0]);
    let _ = (&beta0_basis, &beta1_basis); // unused for now; reserved for potential higher-degree β expansion.

    for (s_idx, spoke) in input.spokes.iter().enumerate() {
        let q = &spoke.boundary_cps;
        let bc = blend_coeffs[s_idx];

        // Build spoke-t-ordered row-1 CPs for A and B sides.
        let patch_a = s_idx;
        let patch_b = (s_idx + 1) % n;
        let mut patch_a_row1 = [Point3::new(0.0, 0.0, 0.0); 8];
        let mut patch_b_row1 = [Point3::new(0.0, 0.0, 0.0); 8];
        for jj in 0..8 {
            patch_a_row1[jj] = seed[patch_a].grid[6][jj];
            // Patch B's VZero reversed: spoke-t = jj/7 → u = 1 - jj/7.
            patch_b_row1[jj] = seed[patch_b].grid[7 - jj][1];
        }

        // Build the per-spoke LHS constant polynomial (xtan_A +
        // α·xtan_B + β·T) at degree 8, with j=1..=6 / k'=1..=6
        // contributions stripped out (carried via per-unknown
        // weights below).
        let mut xtan_a_const = [Vec3::zero(); 8];
        for j in 0..8 {
            let q_j = q[j];
            let row1_j = patch_a_row1[j];
            xtan_a_const[j] = Vec3::new(
                7.0 * (q_j.x - row1_j.x),
                7.0 * (q_j.y - row1_j.y),
                7.0 * (q_j.z - row1_j.z),
            );
        }
        for j in 1..=6 {
            let q_j = q[j];
            xtan_a_const[j] = Vec3::new(7.0 * q_j.x, 7.0 * q_j.y, 7.0 * q_j.z);
        }

        let mut xtan_b_const = [Vec3::zero(); 8];
        for kk in 0..8 {
            let q_k = q[kk];
            let row1_k = patch_b_row1[kk];
            xtan_b_const[kk] = Vec3::new(
                7.0 * (row1_k.x - q_k.x),
                7.0 * (row1_k.y - q_k.y),
                7.0 * (row1_k.z - q_k.z),
            );
        }
        for kk in 1..=6 {
            let q_k = q[kk];
            xtan_b_const[kk] = Vec3::new(-7.0 * q_k.x, -7.0 * q_k.y, -7.0 * q_k.z);
        }

        let mut t_coeffs = Vec::with_capacity(7);
        for jj in 0..7 {
            let a = q[jj];
            let b = q[jj + 1];
            t_coeffs.push(Vec3::new(
                7.0 * (b.x - a.x),
                7.0 * (b.y - a.y),
                7.0 * (b.z - a.z),
            ));
        }
        let t_poly = BernsteinPoly::<Vec3>::new(t_coeffs);

        let alpha_poly = BernsteinPoly::<f64>::new(vec![bc.alpha_0, bc.alpha_1]);
        let beta_poly = BernsteinPoly::<f64>::new(vec![bc.beta_0, bc.beta_1]);

        let beta_t = bernstein_product_scalar_vec3(&beta_poly, &t_poly);

        let xtan_a_const_poly = BernsteinPoly::<Vec3>::new(xtan_a_const.to_vec()).degree_elevate(1);
        let xtan_b_const_poly = BernsteinPoly::<Vec3>::new(xtan_b_const.to_vec());
        let alpha_xtan_b_const = bernstein_product_scalar_vec3(&alpha_poly, &xtan_b_const_poly);
        let beta_t_elev = beta_t.degree_elevate(1);

        let lhs_const = xtan_a_const_poly
            .add_elevated(&alpha_xtan_b_const)
            .add_elevated(&beta_t_elev);

        // Per-unknown weights (9 per unknown, scalar each) at
        // degree 8. Side A: -7·(8-j)/8 at k=j, -7·(j+1)/8 at k=j+1.
        // Side B: +7·(C(7,k')/C(8,k'))·α_0 at k=k', +7·(C(7,k')/C(8,k'+1))·α_1 at k=k'+1.
        let mut weights_a: [[f64; 9]; 6] = [[0.0; 9]; 6];
        for (idx, j) in (1..=6).enumerate() {
            weights_a[idx][j] = -7.0 * (8 - j) as f64 / 8.0;
            weights_a[idx][j + 1] = -7.0 * (j + 1) as f64 / 8.0;
        }
        let mut weights_b: [[f64; 9]; 6] = [[0.0; 9]; 6];
        for (idx, kp) in (1..=6).enumerate() {
            let c_7_kp = binomial(7, kp);
            let c_8_kp = binomial(8, kp);
            let c_8_kpp1 = binomial(8, kp + 1);
            weights_b[idx][kp] = 7.0 * (c_7_kp / c_8_kp) * bc.alpha_0;
            weights_b[idx][kp + 1] = 7.0 * (c_7_kp / c_8_kpp1) * bc.alpha_1;
        }

        // Map the 12 per-spoke unknowns to global variable offsets.
        let mut a_global_offsets = [0usize; 6];
        for (idx, j) in (1..=6).enumerate() {
            a_global_offsets[idx] = var_offset[&(patch_a, 6, j)];
        }
        let mut b_global_offsets = [0usize; 6];
        for (idx, kp) in (1..=6).enumerate() {
            let i = 7 - kp;
            b_global_offsets[idx] = var_offset[&(patch_b, i, 1)];
        }

        // Emit 27 scalar rows.
        for kk in 0..9 {
            let coeff = lhs_const.coeff(kk);
            for c in 0..3 {
                let mut row = vec![0.0_f64; n_unknowns];
                for idx_a in 0..6 {
                    row[a_global_offsets[idx_a] + c] = weights_a[idx_a][kk];
                }
                for idx_b in 0..6 {
                    row[b_global_offsets[idx_b] + c] = weights_b[idx_b][kk];
                }
                let r = match c {
                    0 => -coeff.x,
                    1 => -coeff.y,
                    _ => -coeff.z,
                };
                a_mat.push(row);
                rhs.push(r);
            }
        }
    }

    // 4. Min-norm KKT solve: x = seed + Aᵀ · λ,
    //    (A · Aᵀ + μI) · λ = rhs − A · seed.
    let m = a_mat.len();
    let mut a_seed = vec![0.0_f64; m];
    for r in 0..m {
        let mut s = 0.0;
        for k in 0..n_unknowns {
            s += a_mat[r][k] * seed_vec[k];
        }
        a_seed[r] = s;
    }
    let rhs_adj: Vec<f64> = (0..m).map(|r| rhs[r] - a_seed[r]).collect();

    let mut aat = vec![vec![0.0_f64; m]; m];
    let mut trace = 0.0_f64;
    for r in 0..m {
        for cc in r..m {
            let mut s = 0.0_f64;
            for k in 0..n_unknowns {
                s += a_mat[r][k] * a_mat[cc][k];
            }
            aat[r][cc] = s;
            aat[cc][r] = s;
        }
        trace += aat[r][r];
    }
    let mu = 1e-12 * (trace / (m as f64)).max(1.0);
    for r in 0..m {
        aat[r][r] += mu;
    }
    let Some(lambdas) = solve_dense(aat, rhs_adj) else {
        return Err(DirectConstructionError::CornerTangentPlaneSingular {
            corner_idx: 0,
            kind: CornerKind::Central,
            smallest_sv: 0.0,
            tolerance: options.corner_singularity_eps,
            edge_hint: "ring row-1 KKT singular".into(),
        });
    };

    let mut x = seed_vec;
    for r in 0..m {
        for k in 0..n_unknowns {
            x[k] += lambdas[r] * a_mat[r][k];
        }
    }

    // 5. Apply updated CPs to patches.
    let mut updated_patches = seed.clone();
    for (&(patch_idx, i, j), &off) in &var_offset {
        updated_patches[patch_idx].grid[i][j] = Point3::new(x[off], x[off + 1], x[off + 2]);
    }

    // 6. Per-spoke residual for telemetry.
    let mut per_spoke_residual = Vec::with_capacity(n);
    let mut max_residual = 0.0_f64;
    for s_idx in 0..n {
        let mut spoke_max = 0.0_f64;
        for kk in 0..9 {
            let row_base = 27 * s_idx + 3 * kk;
            let mut lhs = [0.0_f64; 3];
            for c in 0..3 {
                let row = row_base + c;
                let mut v = -rhs[row];
                for (j, coef) in a_mat[row].iter().enumerate() {
                    v += coef * x[j];
                }
                lhs[c] = v;
            }
            let mag = (lhs[0].powi(2) + lhs[1].powi(2) + lhs[2].powi(2)).sqrt();
            if mag > spoke_max {
                spoke_max = mag;
            }
        }
        per_spoke_residual.push(spoke_max);
        if spoke_max > max_residual {
            max_residual = spoke_max;
        }
    }

    Ok(RingRow1Result {
        updated_patches,
        max_residual,
        per_spoke_residual,
    })
}

#[cfg(test)]
#[path = "direct_row1_ring_tests.rs"]
mod tests;
