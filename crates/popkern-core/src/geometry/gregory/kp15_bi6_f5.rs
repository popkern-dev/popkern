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

/// Build the 49 × 49 `F_5` quadratic-form matrix. Thin wrapper
/// over [`f_kappa_quadratic_form`]`(5)` kept for API compatibility.
pub fn f5_quadratic_form() -> [[f64; BI6_CP_COUNT]; BI6_CP_COUNT] {
    f_kappa_quadratic_form(5)
}

/// Build the 49 × 49 `F_κ` quadratic-form matrix for κ ∈ {1, 3,
/// 4, 5}. F_κ(f) := Σ_{i+j=κ} (κ! / (i!·j!)) · ∫∫ (∂ⁱ_s ∂ʲ_t f)²
/// du dv (KP15 eq 14 / math spec §1.4).
///
/// F_5 is the paper's primary choice (KP15 §3.6 step 4). F_4 is
/// used in KP13 n=3; F_3 is KP13's general choice for p; F_1 is
/// Dirichlet energy (KP13 p refinement). See `f5-tuning-plan.md`
/// option C — the F_κ fallback ladder for fold-avoidance on
/// sharp apex geometries.
///
/// # Panics
/// Panics if `kappa` isn't in `{1, 3, 4, 5}`. F_2 is intentionally
/// unsupported (paper notes F*_2 produces visible highlight-line
/// oscillations; see math spec §1.4).
pub fn f_kappa_quadratic_form(kappa: usize) -> [[f64; BI6_CP_COUNT]; BI6_CP_COUNT] {
    let (partials, coeffs) = f_kappa_partials_and_coeffs(kappa);
    assemble_f_kappa_q(&partials, &coeffs)
}

/// Return the `(i, j)` partial indices and `κ!/(i!·j!)` coefficients
/// for F_κ, κ ∈ {1, 3, 4, 5}. Panics otherwise.
fn f_kappa_partials_and_coeffs(kappa: usize) -> (Vec<(usize, usize)>, Vec<f64>) {
    match kappa {
        1 => (vec![(1, 0), (0, 1)], vec![1.0, 1.0]),
        3 => (
            vec![(3, 0), (2, 1), (1, 2), (0, 3)],
            vec![1.0, 3.0, 3.0, 1.0],
        ),
        4 => (
            vec![(4, 0), (3, 1), (2, 2), (1, 3), (0, 4)],
            vec![1.0, 4.0, 6.0, 4.0, 1.0],
        ),
        5 => (F5_PARTIALS.to_vec(), F5_COEFFS.to_vec()),
        _ => panic!("f_kappa: unsupported κ = {kappa} (expected 1, 3, 4, or 5)"),
    }
}

/// Build an F_κ quadratic form with **caller-supplied weight
/// overrides** — the option-E (per-valence tuning) path.
///
/// Starts from the paper's `κ!/(i!·j!)` defaults, then multiplies
/// each term by `weight_multipliers[t]`. Length must equal κ + 1
/// (the number of (i, j) pairs with `i + j = κ`). Panics if the
/// length mismatches or κ ∉ {1, 3, 4, 5}.
///
/// # Example
/// For κ = 5 the partials are ordered [(5,0), (4,1), (3,2), (2,3),
/// (1,4), (0,5)]; a profile like `[0.5, 1.0, 1.0, 1.0, 1.0, 0.5]`
/// halves the pure radial/angular penalty while keeping the mixed
/// terms at paper defaults.
pub fn f_kappa_quadratic_form_weighted(
    kappa: usize,
    weight_multipliers: &[f64],
) -> [[f64; BI6_CP_COUNT]; BI6_CP_COUNT] {
    let (partials, mut coeffs) = f_kappa_partials_and_coeffs(kappa);
    assert_eq!(
        weight_multipliers.len(),
        coeffs.len(),
        "f_kappa_quadratic_form_weighted: expected {} weight multipliers for κ={kappa}, got {}",
        coeffs.len(),
        weight_multipliers.len()
    );
    for (c, &m) in coeffs.iter_mut().zip(weight_multipliers.iter()) {
        *c *= m;
    }
    assemble_f_kappa_q(&partials, &coeffs)
}

/// Internal: assemble the 49×49 F_κ matrix from an explicit set
/// of partial indices + weights. Shared by `f_kappa_quadratic_form`
/// and `f_kappa_quadratic_form_weighted`.
fn assemble_f_kappa_q(
    partials: &[(usize, usize)],
    coeffs: &[f64],
) -> [[f64; BI6_CP_COUNT]; BI6_CP_COUNT] {
    let mut q = [[0.0_f64; BI6_CP_COUNT]; BI6_CP_COUNT];
    let u_basis: Vec<[[f64; 7]; 6]> = GL7_NODES.iter().map(|&u| basis_all_derivs(u)).collect();
    let v_basis: Vec<[[f64; 7]; 6]> = GL7_NODES.iter().map(|&v| basis_all_derivs(v)).collect();
    for (term_idx, &(ip, jp)) in partials.iter().enumerate() {
        let c = coeffs[term_idx];
        for (pu, &wu) in GL7_WEIGHTS.iter().enumerate() {
            let bu_ip = &u_basis[pu][ip];
            for (pv, &wv) in GL7_WEIGHTS.iter().enumerate() {
                let bv_jp = &v_basis[pv][jp];
                let weight = wu * wv * c;
                let mut phi = [0.0_f64; BI6_CP_COUNT];
                for a in 0..7 {
                    for b in 0..7 {
                        phi[a * 7 + b] = bu_ip[a] * bv_jp[b];
                    }
                }
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
