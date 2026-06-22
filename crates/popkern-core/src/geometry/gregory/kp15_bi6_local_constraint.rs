//! Local cross-sector smoothness constraints for the K-P bi-6
//! polar cap — [KP15] §3.2, eq (2) and (3).
//!
//! These are the polynomial identities the cap-construction
//! algorithm must satisfy along every cross-sector edge:
//!
//! ```text
//!   eq (2) (G¹):  ∂_v f̃(u, 0) + ∂_v f(u, 0)
//!                  − 2c(1 − u)² · ∂_u f(u, 0) ≡ 0 on u ∈ [0, 1]
//!
//!   eq (3) (G²):  second-derivative analogue (paper §3.1 line 84)
//! ```
//!
//! Here `f̃` is sector `k`'s bi-6 patch, `f` is sector `k − 1`'s
//! bi-6 patch, both indexed so that their shared boundary lies
//! at `v = 0` of each (i.e. `ṕ_{i, 0} = p̀_{i, 0}` for `i ∈ 0..=6`).
//! The two sectors' `v`-axes point *away* from each other across
//! the shared edge, which is why eq (2) has the `+` sign in front
//! of `∂_v f` instead of the `-` you'd expect from a "smooth
//! patch" identity.
//!
//! # Why this module exists
//!
//! KP15 line 152: "their 7 + 8 BB-coefficients must be set to
//! zero. We solve the corresponding local system of 15 linear
//! equations with the symbolic computation system Maple." The
//! Maple script the paper references is `n6`, which the kernel
//! cannot access (see
//! `docs/decisions/2026-06-22-kp15-bi6-deep-interior-without-maple.md`).
//!
//! Phase B of that ADR's plan begins here: mechanise the constraint
//! enumeration directly from BB-form polynomial arithmetic. The
//! exit criterion for Phase B is that, evaluated on the cap that
//! `assemble_bi6_polar_cap` currently produces, the eq (2)
//! residual at every cross-sector edge vanishes to machine zero
//! across all 7 BB-coefficients (after rule (5) on `p̀_{30}`
//! collapses the degree-7 polynomial to degree 6 — see line 149
//! of `docs/papers/15bi6g2.txt`). This is a positive control on
//! rules (S1)-(S4) of [KP15 §3.2], independent of the F_5 path.
//!
//! [KP15]: docs/papers/15bi6g2.txt

// Phase B: the residual extractor and its strip-derivative
// helpers are not yet wired into a production caller — only the
// Phase B / C closure-verification tests consume them. Phases
// C-E will route the cap assembly's deep-interior derivation
// through this module, removing the allow(dead_code) below.
#![allow(dead_code)]

use super::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::bernstein_poly::{
    BernsteinPoly, PolyCoeff, bernstein_product_scalar, bernstein_product_scalar_vec3,
};
use crate::geometry::{Point3, Vec3};

/// Build `∂_v f(u, 0)` as a Bernstein polynomial of degree 6 in
/// `u`, from a sector's bi-6 BB grid. Uses the standard Bezier
/// derivative rule: the BB-coefs are `6 · (p_{i, 1} − p_{i, 0})`.
fn strip_dot_v(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..7)
        .map(|i| (grid[i][1] - grid[i][0]).scale(6.0))
        .collect();
    BernsteinPoly::new(coeffs)
}

/// Build `∂_u f(u, 0)` as a Bernstein polynomial of degree 5 in
/// `u`, from a sector's bi-6 BB grid. BB-coefs are
/// `6 · (p_{i + 1, 0} − p_{i, 0})`.
fn strip_dot_u(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..6)
        .map(|i| (grid[i + 1][0] - grid[i][0]).scale(6.0))
        .collect();
    BernsteinPoly::new(coeffs)
}

/// [KP15] eq (2) residual at the cross-sector edge between
/// `grid_prev` (= sector `k − 1`, the "previous" `f`) and
/// `grid_curr` (= sector `k`, the "current" `f̃`).
///
/// Returns the 8 Bernstein coefficients of the identity polynomial
///
/// ```text
///   ∂_v f̃(u, 0) + ∂_v f(u, 0) − 2c(1 − u)² · ∂_u f(u, 0)
/// ```
///
/// in degree-7 BB form. (The `∂_v` parts are degree 6 and get
/// implicitly degree-elevated by `add_elevated` to match the
/// degree-7 product `2c(1 − u)² · ∂_u f`.)
///
/// If `assemble_bi6_polar_cap` populated the grids per the paper:
/// 1. Rule (5) on `p̀_{30}` makes the BB coefficient at `u^7`
///    vanish (paper line 149: "the degree of `e` reduces to 5
///    and the left-side in (2) becomes a polynomial of degree 6").
/// 2. Rules (S1) through (S4) jointly determine `ṕ_{i, 1}` for
///    `i ∈ 0..=6`, plus the cross-sector identifications, so the
///    other 7 BB-coefs also vanish.
///
/// Any non-zero residual flags an off-paper choice somewhere in
/// the cap pipeline.
pub fn eq2_g1_residual_bb_degree7(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    params: &Bi6Reparams,
) -> [Vec3; 8] {
    let dv_curr = strip_dot_v(grid_curr);
    let dv_prev = strip_dot_v(grid_prev);
    let du_prev = strip_dot_u(grid_prev);

    // b(u) = 2c · (1 − u)² in Bernstein basis of degree 2.
    // (1 − u)² = 1·B_0^2 + 0·B_1^2 + 0·B_2^2, then scaled by 2c.
    let b = BernsteinPoly::new(vec![2.0 * params.c, 0.0, 0.0]);
    // b(u) · ∂_u f(u, 0): degree 2 + 5 = 7.
    let b_times_du = bernstein_product_scalar_vec3(&b, &du_prev);

    // Sum = ∂_v f̃ + ∂_v f − b · ∂_u f. Lift the two degree-6
    // strips to degree 7 by `add_elevated` so all three live on
    // the common output degree 7.
    let negated: Vec<Vec3> = b_times_du.coeffs().iter().map(|c| c.scale(-1.0)).collect();
    let residual = dv_curr
        .add_elevated(&dv_prev)
        .add_elevated(&BernsteinPoly::new(negated));
    debug_assert_eq!(residual.degree(), 7);

    let mut out = [Vec3::zero(); 8];
    for (k, &c) in residual.coeffs().iter().enumerate() {
        out[k] = c;
    }
    out
}

/// `∂²_v f(u, 0)` as a Bernstein polynomial of degree 6 in `u`.
/// Standard second-derivative BB rule:
/// `c_i = 30 · (p_{i, 2} − 2 · p_{i, 1} + p_{i, 0})`.
/// The `30 = 6 · 5` factor is `n · (n − 1)` for degree `n = 6`.
fn strip_dot_dot_v(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..7)
        .map(|i| {
            let p2 = grid[i][2];
            let p1 = grid[i][1];
            let p0 = grid[i][0];
            let v = (p2 - p1) - (p1 - p0);
            v.scale(30.0)
        })
        .collect();
    BernsteinPoly::new(coeffs)
}

/// `∂²_u f(u, 0)` as a Bernstein polynomial of degree 4 in `u`.
/// `c_i = 30 · (p_{i + 2, 0} − 2 · p_{i + 1, 0} + p_{i, 0})`.
fn strip_dot_dot_u(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..5)
        .map(|i| {
            let p2 = grid[i + 2][0];
            let p1 = grid[i + 1][0];
            let p0 = grid[i][0];
            let v = (p2 - p1) - (p1 - p0);
            v.scale(30.0)
        })
        .collect();
    BernsteinPoly::new(coeffs)
}

/// `∂²_{uv} f(u, 0)` as a Bernstein polynomial of degree 5 in `u`.
/// Mixed second partial at the boundary `v = 0`.
/// `c_i = 36 · (p_{i + 1, 1} − p_{i + 1, 0} − p_{i, 1} + p_{i, 0})`.
/// The factor `36 = 6 · 6` comes from differentiating once in each
/// direction on a bi-6 patch.
fn strip_dot_dot_uv(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..6)
        .map(|i| {
            let p11 = grid[i + 1][1];
            let p10 = grid[i + 1][0];
            let p01 = grid[i][1];
            let p00 = grid[i][0];
            let v = (p11 - p10) - (p01 - p00);
            v.scale(36.0)
        })
        .collect();
    BernsteinPoly::new(coeffs)
}

/// [KP15] eq (3) residual at the cross-sector edge, multiplied by
/// `(1 − u + τu)` to clear the rational denominator in `d(u)` (the
/// paper's `δ(1 − u)² / (1 − u + τu)` term — see eq (4) at paper
/// line 102 and the "After multiplication by 1 − u + τu" remark
/// at line 150).
///
/// With `a := −1, b := 2c(1 − u)², d := δ(1 − u)² / (1 − u + τu),
/// δ := 24c(τ − 1), e := b · (b' − d / 2)` the (φ-multiplied)
/// residual polynomial in `u` is:
///
/// ```text
///   φ · ∂²_v f̃ − φ · ∂²_v f + 4c(1 − u)² φ · ∂²_{uv} f
///                            − 4c²(1 − u)⁴ φ · ∂²_u f
///                            − b · [φ b' − (δ/2)(1 − u)²] · ∂_u f
///                            − δ(1 − u)² · ∂_v f
///                            ≡ 0
/// ```
///
/// where `φ := 1 − u + τu`. Returns the BB coefficients of the
/// polynomial as a `Vec<Vec3>`. The polynomial degree (per term)
/// is at most 9 — paper line 150 claims degree 7 after rule (5)
/// reduces the boundary curve. We return all coefficients up to
/// the actual polynomial degree so the eventual Phase C solver
/// can either trust the paper's claim (asserting BB indices ≥ 8
/// vanish under rule (5)) or take the residual as-is.
pub fn eq3_g2_residual_phi_multiplied(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    params: &Bi6Reparams,
) -> Vec<Vec3> {
    let c = params.c;
    let tau = params.tau;
    let delta = 24.0 * c * (tau - 1.0);

    // Scalar building blocks in BB form on [0, 1].
    let phi = BernsteinPoly::new(vec![1.0, tau]);
    let one_minus_u_sq = BernsteinPoly::new(vec![1.0, 0.0, 0.0]);
    let one_minus_u_quartic = BernsteinPoly::new(vec![1.0, 0.0, 0.0, 0.0, 0.0]);
    let b = BernsteinPoly::new(vec![2.0 * c, 0.0, 0.0]);
    // `b'(u) = -4c(1 - u)` as BB degree 1: at u=0, value -4c; at u=1, 0.
    let b_prime = BernsteinPoly::new(vec![-4.0 * c, 0.0]);

    // φ · b' (degree 2).
    let phi_b_prime = bernstein_product_scalar(&phi, &b_prime);
    // (δ/2) · (1 − u)² (degree 2).
    let half_delta_omu_sq = BernsteinPoly::new(vec![0.5 * delta, 0.0, 0.0]);
    // bracket = φ · b' − (δ/2) · (1 − u)² (degree 2).
    let neg_half_delta: Vec<f64> = half_delta_omu_sq.coeffs().iter().map(|&c| -c).collect();
    let bracket = phi_b_prime.add_elevated(&BernsteinPoly::new(neg_half_delta));

    // Per-sector derivative strips.
    let dvv_curr = strip_dot_dot_v(grid_curr);
    let dvv_prev = strip_dot_dot_v(grid_prev);
    let duu_prev = strip_dot_dot_u(grid_prev);
    let duv_prev = strip_dot_dot_uv(grid_prev);
    let du_prev = strip_dot_u(grid_prev);
    let dv_prev = strip_dot_v(grid_prev);

    // Term 1: φ · ∂²_v f̃.
    let t1 = bernstein_product_scalar_vec3(&phi, &dvv_curr);
    // Term 2: −φ · ∂²_v f.
    let phi_dvv_prev = bernstein_product_scalar_vec3(&phi, &dvv_prev);
    let t2 = BernsteinPoly::new(
        phi_dvv_prev
            .coeffs()
            .iter()
            .map(|c| c.scale(-1.0))
            .collect(),
    );

    // Coefficient polynomial 4c · (1 − u)² · φ (degree 3).
    let four_c_omu_sq = BernsteinPoly::new(vec![4.0 * c, 0.0, 0.0]);
    let four_c_omu_sq_phi = bernstein_product_scalar(&four_c_omu_sq, &phi);
    // Term 3: 4c (1 − u)² φ · ∂²_{uv} f.
    let t3 = bernstein_product_scalar_vec3(&four_c_omu_sq_phi, &duv_prev);

    // Coefficient polynomial −4c² · (1 − u)⁴ · φ (degree 5).
    let neg_4c2_omu_4 = BernsteinPoly::new(
        one_minus_u_quartic
            .coeffs()
            .iter()
            .map(|&v| -4.0 * c * c * v)
            .collect(),
    );
    let neg_4c2_omu_4_phi = bernstein_product_scalar(&neg_4c2_omu_4, &phi);
    // Term 4: −4c² (1 − u)⁴ φ · ∂²_u f.
    let t4 = bernstein_product_scalar_vec3(&neg_4c2_omu_4_phi, &duu_prev);

    // Coefficient polynomial −b · bracket (degree 4).
    let b_times_bracket = bernstein_product_scalar(&b, &bracket);
    let neg_b_times_bracket = BernsteinPoly::new(
        b_times_bracket
            .coeffs()
            .iter()
            .map(|&v| -v)
            .collect::<Vec<f64>>(),
    );
    // Term 5: −b [φ b' − (δ/2)(1 − u)²] · ∂_u f.
    let t5 = bernstein_product_scalar_vec3(&neg_b_times_bracket, &du_prev);

    // Coefficient polynomial −δ · (1 − u)² (degree 2).
    let neg_delta_omu_sq = BernsteinPoly::new(
        one_minus_u_sq
            .coeffs()
            .iter()
            .map(|&v| -delta * v)
            .collect(),
    );
    // Term 6: −δ (1 − u)² · ∂_v f.
    let t6 = bernstein_product_scalar_vec3(&neg_delta_omu_sq, &dv_prev);

    // Sum all six terms, each potentially of a different degree.
    // `add_elevated` lifts the lower-degree polynomial to match.
    let sum = t1
        .add_elevated(&t2)
        .add_elevated(&t3)
        .add_elevated(&t4)
        .add_elevated(&t5)
        .add_elevated(&t6);
    sum.coeffs().to_vec()
}

/// Index a flat BB-coef position `(side, i, j)`. `side = 0` is the
/// previous sector `f`; `side = 1` is the current sector `f̃`.
/// Used as the column index of the Jacobian rows below.
fn flat_grid_index(side: usize, i: usize, j: usize) -> usize {
    debug_assert!(side < 2 && i < 7 && j < 7);
    side * 49 + i * 7 + j
}

/// Build a zero bi-6 grid pair (`[grid_prev, grid_curr]`).
fn empty_grids() -> [[[Point3; 7]; 7]; 2] {
    [[[Point3::new(0.0, 0.0, 0.0); 7]; 7]; 2]
}

/// Extract a dense `8 × 98` Jacobian of the eq (2) BB-residual w.r.t.
/// every BB-coef position across the two sectors. Each column at
/// index `flat_grid_index(side, i, j)` is the residual produced by
/// perturbing `grid[side][i][j]` by the unit vector `(1, 0, 0)`
/// — equivalently the partial derivative `∂R / ∂p_{i, j}^{(side)}`
/// along the x-axis (the y- and z-axis columns equal it by linearity
/// of the residual in each coordinate).
///
/// Returns a `[[f64; 98]; 8]` row-major dense matrix. The entries
/// are exact (no finite-difference noise) because the residual is
/// linear in the BB-coefs.
///
/// This is the building block for Phase C's `M(n, c, τ) · x = b_input`
/// linear system. Each cross-sector edge contributes one such block
/// of 8 rows; the full M stacks `n` of them with sector-relative
/// indexing.
pub fn eq2_jacobian_dense(params: &Bi6Reparams) -> [[f64; 98]; 8] {
    let mut jac = [[0.0_f64; 98]; 8];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb = eq2_g1_residual_bb_degree7(&grids[0], &grids[1], params);
                let col = flat_grid_index(side, i, j);
                for row in 0..8 {
                    jac[row][col] = bb[row].x;
                }
            }
        }
    }
    jac
}

/// Mirror of [`eq2_jacobian_dense`] for the φ-multiplied eq (3).
/// Output row count equals the residual polynomial's BB-vector
/// length (degree + 1 ≤ 10). Caller is responsible for matching
/// rows to the paper's "8 constraints at degree 7" claim — see the
/// notes on `eq3_g2_residual_phi_multiplied`.
pub fn eq3_jacobian_dense(params: &Bi6Reparams) -> Vec<[f64; 98]> {
    // First evaluation establishes the residual polynomial's degree
    // (which depends on whether any terms collapse — at our pre-
    // rule-5, pre-(S1-4) level, all terms are present so the degree
    // is the maximum 9).
    let zero_grids = empty_grids();
    let baseline = eq3_g2_residual_phi_multiplied(&zero_grids[0], &zero_grids[1], params);
    let n_rows = baseline.len();
    let mut jac: Vec<[f64; 98]> = vec![[0.0_f64; 98]; n_rows];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb = eq3_g2_residual_phi_multiplied(&grids[0], &grids[1], params);
                let col = flat_grid_index(side, i, j);
                for (row, v) in bb.iter().enumerate() {
                    jac[row][col] = v.x;
                }
            }
        }
    }
    jac
}

#[cfg(test)]
#[path = "kp15_bi6_local_constraint_tests.rs"]
mod tests;
