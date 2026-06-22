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
use crate::geometry::bernstein_poly::{BernsteinPoly, PolyCoeff, bernstein_product_scalar_vec3};
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

#[cfg(test)]
#[path = "kp15_bi6_local_constraint_tests.rs"]
mod tests;
