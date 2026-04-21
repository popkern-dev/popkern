//! Phase I.9 — K-P 2015 F_5 quality functional for bi-6 patches.
//!
//! Source: [KP15 §3.6 step 4], transcribed in
//! `docs/journal/polar-kp-math.md` §1.4 and §2.5:
//!
//! ```text
//!   F_5(f) = ∫₀¹ ∫₀¹  Σ_{i+j=5} (5! / (i!j!)) · (∂ⁱ_s ∂ʲ_t f)²  ds dt
//!          = ∫∫  (∂⁵_s f)² + 5(∂⁴_s ∂_t f)² + 10(∂³_s ∂²_t f)²
//!               + 10(∂²_s ∂³_t f)² + 5(∂_s ∂⁴_t f)² + (∂⁵_t f)²
//!             ds dt
//! ```
//!
//! Per KP15 step 4, the cap is finalised by minimising F_5 over
//! the free parameters `P` (per §2.5 step 4). F_5 is quadratic in
//! the BB-coefficients, so the minimisation is an SPD linear
//! solve.
//!
//! # Scope of I.9
//!
//! This module provides:
//! - [`f5_quadratic_form`] — the 49×49 PSD matrix `Q` such that
//!   `F_5(f) = flat(f)ᵀ · Q · flat(f)` for the scalar-valued case,
//!   or `Σ_{d∈{x,y,z}} flat(f_d)ᵀ · Q · flat(f_d)` for vector f.
//! - [`f5_evaluate`] — scalar F_5 value for a given bi-6 patch.
//!
//! I.10 will use `Q` to build the global quadratic form over the
//! full free-parameter vector P (by propagating `P → [patch_k]`
//! through rules (S1)–(S4)) and solve.
//!
//! # Flattening convention
//!
//! `[[Point3; 7]; 7]` flattens row-major: `flat[i * 7 + j] =
//! grid[i][j]` for `i, j ∈ 0..7`.
//!
//! # Quadrature
//!
//! F_5 involves 5th-order mixed partials; each squared has max
//! degree 12 in one axis and 2 in the other. 7-point Gauss-
//! Legendre on `[0, 1]` is exact for polynomials up to degree 13,
//! so 7×7 product quadrature integrates F_5 exactly.

// Cross-index loops `phi[a * 7 + b]` and `grid[a][b]` don't
// rewrite cleanly to iterator form.
#![allow(clippy::needless_range_loop)]

use super::bernstein_bi6::{
    bernstein_bi6_all, bernstein_bi6_deriv_all, bernstein_bi6_fifth_deriv_all,
    bernstein_bi6_fourth_deriv_all, bernstein_bi6_second_deriv_all, bernstein_bi6_third_deriv_all,
};
use crate::geometry::Point3;

/// Flattened bi-6 patch size: 7 × 7 = 49 BB-CPs.
pub const BI6_CP_COUNT: usize = 49;

/// 7-point Gauss-Legendre nodes on `[0, 1]`.
const GL7_NODES: [f64; 7] = [
    0.02544604382862074, // (1 + x₁)/2 where x₁ ≈ -0.9491079123427585
    0.12923440720030277,
    0.29707742431130146,
    0.5,
    0.7029225756886986,
    0.8707655927996972,
    0.9745539561713793,
];

/// 7-point Gauss-Legendre weights on `[0, 1]` (half the
/// `[-1, 1]` weights).
const GL7_WEIGHTS: [f64; 7] = [
    0.06474248308443484,
    0.13985269574463833,
    0.19091502525255947,
    0.2089795918367347,
    0.19091502525255947,
    0.13985269574463833,
    0.06474248308443484,
];

/// Scalar weights `5! / (i!·j!)` in front of `(∂ⁱ_s ∂ʲ_t f)²` for
/// the six terms in `F_5`, ordered as
/// `(i, j) ∈ {(5,0), (4,1), (3,2), (2,3), (1,4), (0,5)}`.
const F5_COEFFS: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];

/// Build the 49 × 49 `F_5` quadratic-form matrix.
/// Symmetric, positive semi-definite; invariant across all bi-6
/// patches (depends only on Bernstein basis structure).
///
/// # Usage
/// ```ignore
/// let q = f5_quadratic_form();
/// let flat = flatten_patch(&grid);  // 49-long per coord
/// let f5 = flat.iter().zip(/* Q × flat */).map(|(f, qf)| f * qf).sum();
/// ```
///
/// Prefer [`f5_evaluate`] for the scalar case.
pub fn f5_quadratic_form() -> [[f64; BI6_CP_COUNT]; BI6_CP_COUNT] {
    let mut q = [[0.0_f64; BI6_CP_COUNT]; BI6_CP_COUNT];

    // Precompute all six 1D derivative basis vectors at each GL
    // node (for both u and v axes); 7 × 6 × 7 = 294 f64 each.
    // Indices: [gl_node_idx][deriv_order] → [f64; 7]
    let u_basis: Vec<[[f64; 7]; 6]> = GL7_NODES.iter().map(|&u| basis_all_derivs(u)).collect();
    let v_basis: Vec<[[f64; 7]; 6]> = GL7_NODES.iter().map(|&v| basis_all_derivs(v)).collect();

    // For each (i_partial, j_partial) with i+j = 5, accumulate
    // Σ_{quad} w_u · w_v · c_{ij} · outer(phi, phi) where
    // phi_{ab} = B_a^{(i_partial)}(u) · B_b^{(j_partial)}(v).
    for (term_idx, &(ip, jp)) in F5_PARTIALS.iter().enumerate() {
        let c = F5_COEFFS[term_idx];
        for (pu, &wu) in GL7_WEIGHTS.iter().enumerate() {
            let bu_ip = &u_basis[pu][ip];
            for (pv, &wv) in GL7_WEIGHTS.iter().enumerate() {
                let bv_jp = &v_basis[pv][jp];
                let weight = wu * wv * c;

                // Build the basis-derivative vector phi ∈ ℝ⁴⁹
                // and accumulate weight · phi · phiᵀ into q.
                let mut phi = [0.0_f64; BI6_CP_COUNT];
                for a in 0..7 {
                    for b in 0..7 {
                        phi[a * 7 + b] = bu_ip[a] * bv_jp[b];
                    }
                }
                // Symmetric rank-1 update: q += weight · phi · phiᵀ.
                for r in 0..BI6_CP_COUNT {
                    let pr = phi[r];
                    if pr == 0.0 {
                        continue;
                    }
                    let w_pr = weight * pr;
                    for (c_idx, &pc) in phi.iter().enumerate().take(BI6_CP_COUNT) {
                        q[r][c_idx] += w_pr * pc;
                    }
                }
            }
        }
    }
    q
}

/// The six `(i_partial, j_partial)` pairs with `i + j = 5`, in
/// `F5_COEFFS` order.
const F5_PARTIALS: [(usize, usize); 6] = [(5, 0), (4, 1), (3, 2), (2, 3), (1, 4), (0, 5)];

/// Given a parameter `t ∈ [0, 1]`, return all six derivative
/// orders 0..=5 of the 7 bi-6 basis values. Output
/// `[deriv_order][basis_idx]`.
fn basis_all_derivs(t: f64) -> [[f64; 7]; 6] {
    [
        bernstein_bi6_all(t),
        bernstein_bi6_deriv_all(t),
        bernstein_bi6_second_deriv_all(t),
        bernstein_bi6_third_deriv_all(t),
        bernstein_bi6_fourth_deriv_all(t),
        bernstein_bi6_fifth_deriv_all(t),
    ]
}

/// Evaluate `F_5` for a bi-6 patch. Sums across x/y/z coordinates.
///
/// Uses 7×7 Gauss-Legendre quadrature directly (no matrix build),
/// so this is cheap — O(49 · 36) per quadrature point × 49
/// quadrature points ≈ 86k multiplies.
pub fn f5_evaluate(grid: &[[Point3; 7]; 7]) -> f64 {
    let u_basis: Vec<[[f64; 7]; 6]> = GL7_NODES.iter().map(|&u| basis_all_derivs(u)).collect();
    let v_basis: Vec<[[f64; 7]; 6]> = GL7_NODES.iter().map(|&v| basis_all_derivs(v)).collect();

    let mut total = 0.0_f64;
    for (term_idx, &(ip, jp)) in F5_PARTIALS.iter().enumerate() {
        let c = F5_COEFFS[term_idx];
        for (pu, &wu) in GL7_WEIGHTS.iter().enumerate() {
            let bu_ip = &u_basis[pu][ip];
            for (pv, &wv) in GL7_WEIGHTS.iter().enumerate() {
                let bv_jp = &v_basis[pv][jp];
                // ∂^(ip+jp) f at (u_p, v_p) = Σ_{a, b} B_a^{(ip)}(u) · B_b^{(jp)}(v) · grid[a][b]
                // Three coord accumulators.
                let mut dx = 0.0_f64;
                let mut dy = 0.0_f64;
                let mut dz = 0.0_f64;
                for a in 0..7 {
                    if bu_ip[a] == 0.0 {
                        continue;
                    }
                    let wa = bu_ip[a];
                    for b in 0..7 {
                        let w_ab = wa * bv_jp[b];
                        let p = grid[a][b];
                        dx += w_ab * p.x;
                        dy += w_ab * p.y;
                        dz += w_ab * p.z;
                    }
                }
                total += wu * wv * c * (dx * dx + dy * dy + dz * dz);
            }
        }
    }
    total
}

#[cfg(test)]
#[path = "kp15_bi6_f5_tests.rs"]
mod tests;
