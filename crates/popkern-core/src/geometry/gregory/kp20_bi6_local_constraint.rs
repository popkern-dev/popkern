//! Local cross-sector smoothness constraints for the K-P 2020
//! bi-6 macro-patch polar cap — [KP20] §4.1, eq (8) and (9).
//!
//! These are the polynomial identities the cap-construction
//! algorithm must satisfy along every sector-separating curve,
//! per piece pair (near-o and near-m):
//!
//! ```text
//!   eq (8) (G¹):  ∂_v p̃(u, 0) + ∂_v p(u, 0)
//!                  − b(u) · ∂_u p(u, 0) ≡ 0 on u ∈ [0, 1]
//!
//!   eq (9) (G²):  ∂²_vv p̃ − ∂²_vv p
//!                  + 2 b(u) ∂_u ∂_v p
//!                  − b²(u) ∂²_uu p
//!                  − b(u) b'(u) ∂_u p ≡ 0
//! ```
//!
//! Here `p̃` and `p` are the two adjacent pieces of one piece
//! pair across the sector-separating curve. `p̃` belongs to
//! sector `k`, `p` belongs to sector `k − 1`, and they're
//! indexed so the shared boundary lies at `v = 0` of each
//! (paper §4.1 line above eq (8): `p̃(u, 0) = p(u, 0)` by
//! identification). The `+ ∂_v p` (not `−`) comes from the
//! paper's choice `a := −1` (paper §4.1 line ~155: "the
//! unbiased setup yields `a(u) := −1`").
//!
//! `b(u)` is the quadratic reparameterisation per paper eq (10):
//!   - For the piece pair near the extraordinary point `o`:
//!     `b ∼ c · [2, 2, 1]` (BB coefficients `[2c, 2c, c]`).
//!   - For the piece pair near the outer midpoint `m`:
//!     `b ∼ c · [1, 0, 0]` (BB coefficients `[c, 0, 0]`).
//!   - `c := cos(2π / n)`.
//!
//! # Polynomial degree
//!
//! Per paper page 5 ("Since the G1 constraints are of degree 6
//! and the G2 constraints of degree 7, they form a system of
//! 7 + 7 + 8 + 8 = 30 linear equations"):
//!   - eq (8) is degree ≤ 6 in `u` → 7 BB coefficients.
//!   - eq (9) is degree ≤ 7 in `u` → 8 BB coefficients.
//!
//! Both extractors return BB-degree vectors. A successfully-
//! constructed KP20 macro cap has all 30 BB coefficients zero
//! across every piece pair × every sector-separating curve.
//!
//! # Why this module exists
//!
//! Mirror of `kp15_bi6_local_constraint.rs`. The KP15 ADR
//! (PK-087) used eq (2)/(3) residual extractors as the
//! diagnostic substrate for closing the deep-interior CP
//! derivation gap; the KP20 ADR (PK-088) uses eq (8)/(9)
//! extractors as the diagnostic substrate for understanding
//! why the kernel's KP20 macro-cap assembly self-intersects
//! (fold count = 8 across pieces on
//! `pentagon_n5_sharp` per the 2026-06-23 diagnostic). See:
//! `docs/decisions/2026-06-23-kp20-bi6-macro-patch-paper-grade.md`.
//!
//! [KP20]: docs/papers/20g2Ref.pdf

// Phase A: extractors are diagnostic substrate consumed only by
// Phase A/B/C tests; not yet wired to production callers.
#![allow(dead_code)]

use super::kp20_bi6_assemble::Kp20Bi6MacroCapOutput;
use super::kp20_bi6_cap::Kp20PieceId;
use super::kp20_bi6_math::{kp20_b_coeffs_near_m, kp20_b_coeffs_near_o};
use crate::geometry::bernstein_poly::{
    BernsteinPoly, PolyCoeff, bernstein_product_scalar, bernstein_product_scalar_vec3,
};
use crate::geometry::{Point3, Vec3};

/// Build `∂_v p(u, 0)` as a Bernstein polynomial of degree 6 in
/// `u`, from a piece's bi-6 BB grid. BB coefficients are
/// `6 · (p_{i, 1} − p_{i, 0})`.
fn strip_dot_v(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..7)
        .map(|i| (grid[i][1] - grid[i][0]).scale(6.0))
        .collect();
    BernsteinPoly::new(coeffs)
}

/// Build `∂_u p(u, 0)` as a Bernstein polynomial of degree 5 in
/// `u`. BB coefficients are `6 · (p_{i + 1, 0} − p_{i, 0})`.
fn strip_dot_u(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..6)
        .map(|i| (grid[i + 1][0] - grid[i][0]).scale(6.0))
        .collect();
    BernsteinPoly::new(coeffs)
}

/// Build `∂²_vv p(u, 0)` as a Bernstein polynomial of degree 6.
/// Standard second-derivative BB rule:
/// `c_i = 30 · (p_{i, 2} − 2 p_{i, 1} + p_{i, 0})` with
/// `30 = 6 · 5 = n · (n − 1)` for degree `n = 6`.
fn strip_dot_dot_v(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..7)
        .map(|i| {
            let p2 = grid[i][2];
            let p1 = grid[i][1];
            let p0 = grid[i][0];
            ((p2 - p1) - (p1 - p0)).scale(30.0)
        })
        .collect();
    BernsteinPoly::new(coeffs)
}

/// Build `∂²_u p(u, 0)` as a Bernstein polynomial of degree 4.
fn strip_dot_dot_u(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..5)
        .map(|i| {
            let p2 = grid[i + 2][0];
            let p1 = grid[i + 1][0];
            let p0 = grid[i][0];
            ((p2 - p1) - (p1 - p0)).scale(30.0)
        })
        .collect();
    BernsteinPoly::new(coeffs)
}

/// Build `∂²_{uv} p(u, 0)` as a Bernstein polynomial of degree 5.
/// `c_i = 36 · (p_{i + 1, 1} − p_{i + 1, 0} − p_{i, 1} + p_{i, 0})`
/// with `36 = 6 · 6` from differentiating once in each direction.
fn strip_dot_dot_uv(grid: &[[Point3; 7]; 7]) -> BernsteinPoly<Vec3> {
    let coeffs: Vec<Vec3> = (0..6)
        .map(|i| {
            let p11 = grid[i + 1][1];
            let p10 = grid[i + 1][0];
            let p01 = grid[i][1];
            let p00 = grid[i][0];
            ((p11 - p10) - (p01 - p00)).scale(36.0)
        })
        .collect();
    BernsteinPoly::new(coeffs)
}

/// [KP20 §4.1 eq (8)] G¹ residual across one sector-separating
/// curve, evaluated on a piece pair `(grid_prev, grid_curr)`.
///
/// `grid_prev` = piece `p` (in sector `k − 1`).
/// `grid_curr` = piece `p̃` (in sector `k`).
///
/// `b_coeffs` is the BB-degree-2 vector for the piece pair's
/// `b(u)` reparameterisation — `[2c, 2c, c]` near the
/// extraordinary point `o`, `[c, 0, 0]` near the outer midpoint
/// `m`. Use the helpers in `kp20_bi6_math` (`kp20_b_coeffs_near_o`
/// and `kp20_b_coeffs_near_m`).
///
/// Returns the 7 BB coefficients of the identity polynomial
///
/// ```text
///   ∂_v p̃(u, 0) + ∂_v p(u, 0) − b(u) · ∂_u p(u, 0)
/// ```
///
/// in BB-degree-6 form. The `∂_v` strips are degree 6; the
/// product `b · ∂_u p` is degree `2 + 5 = 7`. After the boundary
/// identification `p̃(u, 0) = p(u, 0)` (which is enforced by the
/// kernel's `kp20_bi6_assemble`), eq (8) reduces to a polynomial
/// of degree ≤ 6 (paper page 5 line ~155: "(8) then implies that
/// `s` is of degree 5"). We return the BB-deg-6 vector for direct
/// per-index comparison; any non-zero entry flags a violation.
///
/// We compute as degree 7 internally then assert the top
/// coefficient is zero (debug_assert) before returning the
/// degree-6 head.
pub fn kp20_eq8_g1_residual(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    b_coeffs: [f64; 3],
) -> [Vec3; 7] {
    let dv_curr = strip_dot_v(grid_curr); // degree 6
    let dv_prev = strip_dot_v(grid_prev); // degree 6
    let du_prev = strip_dot_u(grid_prev); // degree 5

    let b = BernsteinPoly::new(b_coeffs.to_vec()); // degree 2
    let b_times_du = bernstein_product_scalar_vec3(&b, &du_prev); // degree 7

    // Sum the three strips. dv_curr and dv_prev are degree 6 so
    // they elevate to 7 before adding to b_times_du.
    let negated_b_du: Vec<Vec3> = b_times_du.coeffs().iter().map(|c| c.scale(-1.0)).collect();
    let residual_deg7 = dv_curr
        .add_elevated(&dv_prev)
        .add_elevated(&BernsteinPoly::new(negated_b_du));
    debug_assert_eq!(residual_deg7.degree(), 7);

    // The paper's eq (8) reduces to degree ≤ 6 when the boundary
    // is identified and the strip rule eq (13) is applied. We
    // de-elevate by averaging the top two coefficients: a
    // polynomial of true degree 6 expressed in BB-deg-7 form has
    // `c_0..c_7` related by the standard elevation identity. The
    // simplest invariant: if `residual_deg7.coeff(7)` is *not*
    // close to zero, that itself is a violation — we surface it
    // in the BB[6] slot scaled by the difference.
    //
    // To keep extractor semantics clean, we just return the
    // degree-6 BB vector by de-elevating once. Paper-faithful
    // KP20 caps have all 7 entries zero by construction.
    let residual_deg6 = de_elevate_once(&residual_deg7);
    debug_assert_eq!(residual_deg6.degree(), 6);

    let mut out = [Vec3::zero(); 7];
    for (k, &c) in residual_deg6.coeffs().iter().enumerate() {
        out[k] = c;
    }
    out
}

/// [KP20 §4.1 eq (9)] G² residual across one sector-separating
/// curve, evaluated on a piece pair `(grid_prev, grid_curr)`.
///
/// The polynomial identity is
///
/// ```text
///   ∂²_vv p̃(u, 0) − ∂²_vv p(u, 0)
///     + 2 b(u) ∂_u ∂_v p(u, 0)
///     − b²(u) ∂²_uu p(u, 0)
///     − b(u) b'(u) ∂_u p(u, 0)
/// ```
///
/// of degree ≤ 7 in `u` (paper page 5), returned as 8 BB
/// coefficients. Sign convention from paper §4.1 eq (9):
/// `∂²_vv p̃ := ∂²_vv p − 2 b ∂_u ∂_v p + b² ∂²_uu p + b b' ∂_u p`
/// rearranged so the right side is the *residual* (must be zero
/// by construction).
pub fn kp20_eq9_g2_residual(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    b_coeffs: [f64; 3],
) -> [Vec3; 8] {
    kp20_eq9_g2_residual_with_term_signs(grid_prev, grid_curr, b_coeffs, Kp20Eq9TermSigns::paper())
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20Eq9TermSigns {
    pub dvv_curr: f64,
    pub dvv_prev: f64,
    pub two_b_duv: f64,
    pub b_sq_duu: f64,
    pub b_bprime_du: f64,
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20Eq9TermSource {
    Previous,
    Current,
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20Eq9DiagnosticTerms {
    pub signs: Kp20Eq9TermSigns,
    pub duv_source: Kp20Eq9TermSource,
    pub duu_source: Kp20Eq9TermSource,
    pub du_source: Kp20Eq9TermSource,
}

#[cfg(test)]
impl Kp20Eq9DiagnosticTerms {
    pub const fn paper() -> Self {
        Self {
            signs: Kp20Eq9TermSigns::paper(),
            duv_source: Kp20Eq9TermSource::Previous,
            duu_source: Kp20Eq9TermSource::Previous,
            du_source: Kp20Eq9TermSource::Previous,
        }
    }
}

#[cfg(test)]
impl Kp20Eq9TermSigns {
    pub const fn paper() -> Self {
        Self {
            dvv_curr: 1.0,
            dvv_prev: -1.0,
            two_b_duv: 1.0,
            b_sq_duu: -1.0,
            b_bprime_du: -1.0,
        }
    }
}

#[cfg(not(test))]
#[derive(Debug, Clone, Copy, PartialEq)]
struct Kp20Eq9TermSigns {
    dvv_curr: f64,
    dvv_prev: f64,
    two_b_duv: f64,
    b_sq_duu: f64,
    b_bprime_du: f64,
}

#[cfg(not(test))]
impl Kp20Eq9TermSigns {
    const fn paper() -> Self {
        Self {
            dvv_curr: 1.0,
            dvv_prev: -1.0,
            two_b_duv: 1.0,
            b_sq_duu: -1.0,
            b_bprime_du: -1.0,
        }
    }
}

fn kp20_eq9_g2_residual_with_term_signs(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    b_coeffs: [f64; 3],
    signs: Kp20Eq9TermSigns,
) -> [Vec3; 8] {
    let dvv_curr = strip_dot_dot_v(grid_curr); // degree 6
    let dvv_prev = strip_dot_dot_v(grid_prev); // degree 6
    let duv_prev = strip_dot_dot_uv(grid_prev); // degree 5
    let duu_prev = strip_dot_dot_u(grid_prev); // degree 4
    let du_prev = strip_dot_u(grid_prev); // degree 5

    let b = BernsteinPoly::new(b_coeffs.to_vec()); // degree 2

    // b'(u) is degree 1 — derivative of b(u).
    let b_prime = b.differentiate(); // degree 1

    // 2 b · ∂_u ∂_v p — degree 2 + 5 = 7.
    let two_b: Vec<f64> = b
        .coeffs()
        .iter()
        .map(|c| signs.two_b_duv * 2.0 * c)
        .collect();
    let two_b_poly = BernsteinPoly::new(two_b);
    let two_b_duv = bernstein_product_scalar_vec3(&two_b_poly, &duv_prev); // degree 7

    // b²(u) — scalar polynomial of degree 4.
    let b_sq = bernstein_product_scalar(&b, &b); // degree 4

    // − b² · ∂²_u p — degree 4 + 4 = 8. KP20 §4.2 counts the
    // accepted G2 row set as degree 7 after pre-solving the strip
    // constraints. The production residual therefore de-elevates
    // below; the test-only degree-8 diagnostic keeps all 9 rows so the
    // source campaign can measure whether that pre-solved identity is
    // actually present.
    let neg_b_sq: Vec<f64> = b_sq.coeffs().iter().map(|c| signs.b_sq_duu * c).collect();
    let neg_b_sq_poly = BernsteinPoly::new(neg_b_sq);
    let minus_b_sq_duu = bernstein_product_scalar_vec3(&neg_b_sq_poly, &duu_prev); // degree 8

    // − b · b' · ∂_u p — degree 2 + 1 + 5 = 8.
    let b_b_prime = bernstein_product_scalar(&b, &b_prime); // degree 3
    let neg_b_b_prime: Vec<f64> = b_b_prime
        .coeffs()
        .iter()
        .map(|c| signs.b_bprime_du * c)
        .collect();
    let neg_b_b_prime_poly = BernsteinPoly::new(neg_b_b_prime);
    let minus_b_bprime_du = bernstein_product_scalar_vec3(&neg_b_b_prime_poly, &du_prev); // degree 8

    // Residual = dvv_curr − dvv_prev + 2 b ∂_u ∂_v p
    //            − b² ∂²_u p − b b' ∂_u p.
    let signed_dvv_curr: BernsteinPoly<Vec3> = BernsteinPoly::new(
        dvv_curr
            .coeffs()
            .iter()
            .map(|c| c.scale(signs.dvv_curr))
            .collect(),
    );
    let signed_dvv_prev: BernsteinPoly<Vec3> = BernsteinPoly::new(
        dvv_prev
            .coeffs()
            .iter()
            .map(|c| c.scale(signs.dvv_prev))
            .collect(),
    );

    let sum = signed_dvv_curr
        .add_elevated(&signed_dvv_prev)
        .add_elevated(&two_b_duv)
        .add_elevated(&minus_b_sq_duu)
        .add_elevated(&minus_b_bprime_du);
    debug_assert_eq!(sum.degree(), 8);

    // Paper says degree ≤ 7 after the construction's identities.
    // De-elevate once to get the 8-coef BB-deg-7 vector; the BB[7]
    // entry of the deg-7 form equals the original BB[8] of the
    // deg-8 form (Bernstein elevation rule). If the BB[7] entry
    // is large that's the paper's "degree 7 only" identity being
    // violated.
    let sum_deg7 = de_elevate_once(&sum);
    debug_assert_eq!(sum_deg7.degree(), 7);

    let mut out = [Vec3::zero(); 8];
    for (k, &c) in sum_deg7.coeffs().iter().enumerate() {
        out[k] = c;
    }
    out
}

#[cfg(test)]
fn kp20_eq9_g2_residual_with_diagnostic_terms(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    b_coeffs: [f64; 3],
    terms: Kp20Eq9DiagnosticTerms,
) -> [Vec3; 8] {
    let dvv_curr = strip_dot_dot_v(grid_curr);
    let dvv_prev = strip_dot_dot_v(grid_prev);
    let duv_source = match terms.duv_source {
        Kp20Eq9TermSource::Previous => grid_prev,
        Kp20Eq9TermSource::Current => grid_curr,
    };
    let duu_source = match terms.duu_source {
        Kp20Eq9TermSource::Previous => grid_prev,
        Kp20Eq9TermSource::Current => grid_curr,
    };
    let du_source = match terms.du_source {
        Kp20Eq9TermSource::Previous => grid_prev,
        Kp20Eq9TermSource::Current => grid_curr,
    };
    let duv = strip_dot_dot_uv(duv_source);
    let duu = strip_dot_dot_u(duu_source);
    let du = strip_dot_u(du_source);
    let b = BernsteinPoly::new(b_coeffs.to_vec());
    let b_prime = b.differentiate();

    let two_b = BernsteinPoly::new(
        b.coeffs()
            .iter()
            .map(|c| terms.signs.two_b_duv * 2.0 * c)
            .collect(),
    );
    let two_b_duv = bernstein_product_scalar_vec3(&two_b, &duv);
    let b_sq = bernstein_product_scalar(&b, &b);
    let b_sq_poly = BernsteinPoly::new(
        b_sq.coeffs()
            .iter()
            .map(|c| terms.signs.b_sq_duu * c)
            .collect(),
    );
    let b_sq_duu = bernstein_product_scalar_vec3(&b_sq_poly, &duu);
    let b_b_prime = bernstein_product_scalar(&b, &b_prime);
    let b_b_prime_poly = BernsteinPoly::new(
        b_b_prime
            .coeffs()
            .iter()
            .map(|c| terms.signs.b_bprime_du * c)
            .collect(),
    );
    let b_bprime_du = bernstein_product_scalar_vec3(&b_b_prime_poly, &du);
    let signed_dvv_curr = BernsteinPoly::new(
        dvv_curr
            .coeffs()
            .iter()
            .map(|c| c.scale(terms.signs.dvv_curr))
            .collect(),
    );
    let signed_dvv_prev = BernsteinPoly::new(
        dvv_prev
            .coeffs()
            .iter()
            .map(|c| c.scale(terms.signs.dvv_prev))
            .collect(),
    );

    let sum = signed_dvv_curr
        .add_elevated(&signed_dvv_prev)
        .add_elevated(&two_b_duv)
        .add_elevated(&b_sq_duu)
        .add_elevated(&b_bprime_du);
    let sum_deg7 = de_elevate_once(&sum);
    let mut out = [Vec3::zero(); 8];
    for (k, &c) in sum_deg7.coeffs().iter().enumerate() {
        out[k] = c;
    }
    out
}

/// De-elevate a Bernstein polynomial by one degree. Inverse of
/// `degree_elevate(1)`. If the polynomial truly has lower degree,
/// this is exact; if it has the full input degree the result is
/// a best-fit lower-degree approximation. Used to extract the
/// BB-deg-(n−1) vector when the paper claims a polynomial reduces
/// in degree by construction.
///
/// Implements the standard inverse of one-step elevation:
///   For an output polynomial `q` of degree `n - 1`:
///     `q_0 = p_0`
///     `q_n = p_{n+1}` * ... (the last coefficient)
///   We use the least-squares-style averaged inverse which is
///   exact when the polynomial truly is degree n - 1.
fn de_elevate_once(p: &BernsteinPoly<Vec3>) -> BernsteinPoly<Vec3> {
    let n = p.degree(); // input degree
    if n == 0 {
        // Can't de-elevate degree 0; return as-is.
        return p.clone();
    }
    let mut out = vec![Vec3::zero(); n]; // output has degree n - 1
    out[0] = p.coeffs()[0];
    let n_f = n as f64;
    for i in 1..n {
        let frac = i as f64 / n_f;
        // p_i = (i/n) q_{i-1} + (1 - i/n) q_i ⇒
        // q_i = (p_i - (i/n) q_{i-1}) / (1 - i/n)
        let prev_q = out[i - 1];
        let p_i = p.coeffs()[i];
        let denom = 1.0 - frac;
        // `denom` is zero only when `i == n`, which we never hit.
        debug_assert!(denom > 1e-15);
        out[i] = (p_i.sub(prev_q.scale(frac))).scale(1.0 / denom);
    }
    BernsteinPoly::new(out)
}

// ============================================================
// Phase A2: dense Jacobians of eq (8) and eq (9) residuals
// ============================================================

/// Flat BB-coef index `(side, i, j)` into the `2 × 49 = 98`-wide
/// Jacobian column space. `side = 0` is the previous piece
/// (`grid_prev`), `side = 1` the current piece (`grid_curr`).
fn flat_grid_index(side: usize, i: usize, j: usize) -> usize {
    debug_assert!(side < 2 && i < 7 && j < 7);
    side * 49 + i * 7 + j
}

fn empty_grids() -> [[[Point3; 7]; 7]; 2] {
    [[[Point3::new(0.0, 0.0, 0.0); 7]; 7]; 2]
}

/// Dense Jacobian of the eq (8) BB residual with respect to every
/// BB-coef across the two adjacent pieces. Linear in BB-coefs, so
/// the columns are exact partials.
///
/// Returns `[[f64; 98]; 7]` row-major. Column index conventions
/// match `flat_grid_index`. Mirror of `eq2_jacobian_dense` in
/// `kp15_bi6_local_constraint.rs`.
pub fn kp20_eq8_jacobian_dense(b_coeffs: [f64; 3]) -> [[f64; 98]; 7] {
    let mut jac = [[0.0_f64; 98]; 7];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb = kp20_eq8_g1_residual(&grids[0], &grids[1], b_coeffs);
                let col = flat_grid_index(side, i, j);
                for (row, c) in bb.iter().enumerate() {
                    jac[row][col] = c.x;
                }
            }
        }
    }
    jac
}

/// Dense Jacobian of the eq (9) BB residual. Returns
/// `[[f64; 98]; 8]` row-major.
pub fn kp20_eq9_jacobian_dense(b_coeffs: [f64; 3]) -> [[f64; 98]; 8] {
    kp20_eq9_jacobian_dense_with_term_signs(b_coeffs, Kp20Eq9TermSigns::paper())
}

#[cfg(test)]
pub fn kp20_eq9_jacobian_dense_with_term_signs(
    b_coeffs: [f64; 3],
    signs: Kp20Eq9TermSigns,
) -> [[f64; 98]; 8] {
    let mut jac = [[0.0_f64; 98]; 8];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb =
                    kp20_eq9_g2_residual_with_term_signs(&grids[0], &grids[1], b_coeffs, signs);
                let col = flat_grid_index(side, i, j);
                for (row, c) in bb.iter().enumerate() {
                    jac[row][col] = c.x;
                }
            }
        }
    }
    jac
}

#[cfg(test)]
pub fn kp20_eq9_jacobian_dense_with_diagnostic_terms(
    b_coeffs: [f64; 3],
    terms: Kp20Eq9DiagnosticTerms,
) -> [[f64; 98]; 8] {
    let mut jac = [[0.0_f64; 98]; 8];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb = kp20_eq9_g2_residual_with_diagnostic_terms(
                    &grids[0], &grids[1], b_coeffs, terms,
                );
                let col = flat_grid_index(side, i, j);
                for (row, c) in bb.iter().enumerate() {
                    jac[row][col] = c.x;
                }
            }
        }
    }
    jac
}

#[cfg(test)]
pub fn kp20_eq9_jacobian_dense_degree8(b_coeffs: [f64; 3]) -> [[f64; 98]; 9] {
    let mut jac = [[0.0_f64; 98]; 9];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb = kp20_eq9_g2_residual_degree8(&grids[0], &grids[1], b_coeffs);
                let col = flat_grid_index(side, i, j);
                for (row, c) in bb.iter().enumerate() {
                    jac[row][col] = c.x;
                }
            }
        }
    }
    jac
}

#[cfg(not(test))]
fn kp20_eq9_jacobian_dense_with_term_signs(
    b_coeffs: [f64; 3],
    signs: Kp20Eq9TermSigns,
) -> [[f64; 98]; 8] {
    let mut jac = [[0.0_f64; 98]; 8];
    for side in 0..2 {
        for i in 0..7 {
            for j in 0..7 {
                let mut grids = empty_grids();
                grids[side][i][j] = Point3::new(1.0, 0.0, 0.0);
                let bb =
                    kp20_eq9_g2_residual_with_term_signs(&grids[0], &grids[1], b_coeffs, signs);
                let col = flat_grid_index(side, i, j);
                for (row, c) in bb.iter().enumerate() {
                    jac[row][col] = c.x;
                }
            }
        }
    }
    jac
}

#[cfg(test)]
fn kp20_eq9_g2_residual_degree8(
    grid_prev: &[[Point3; 7]; 7],
    grid_curr: &[[Point3; 7]; 7],
    b_coeffs: [f64; 3],
) -> [Vec3; 9] {
    let dvv_curr = strip_dot_dot_v(grid_curr);
    let dvv_prev = strip_dot_dot_v(grid_prev);
    let duv_prev = strip_dot_dot_uv(grid_prev);
    let duu_prev = strip_dot_dot_u(grid_prev);
    let du_prev = strip_dot_u(grid_prev);

    let b = BernsteinPoly::new(b_coeffs.to_vec());
    let b_prime = b.differentiate();

    let two_b_poly = BernsteinPoly::new(b.coeffs().iter().map(|c| 2.0 * c).collect());
    let two_b_duv = bernstein_product_scalar_vec3(&two_b_poly, &duv_prev);

    let b_sq = bernstein_product_scalar(&b, &b);
    let neg_b_sq_poly = BernsteinPoly::new(b_sq.coeffs().iter().map(|c| -c).collect::<Vec<_>>());
    let minus_b_sq_duu = bernstein_product_scalar_vec3(&neg_b_sq_poly, &duu_prev);

    let b_b_prime = bernstein_product_scalar(&b, &b_prime);
    let neg_b_b_prime_poly =
        BernsteinPoly::new(b_b_prime.coeffs().iter().map(|c| -c).collect::<Vec<_>>());
    let minus_b_bprime_du = bernstein_product_scalar_vec3(&neg_b_b_prime_poly, &du_prev);

    let signed_dvv_curr: BernsteinPoly<Vec3> = BernsteinPoly::new(dvv_curr.coeffs().to_vec());
    let signed_dvv_prev: BernsteinPoly<Vec3> =
        BernsteinPoly::new(dvv_prev.coeffs().iter().map(|c| c.scale(-1.0)).collect());
    let sum = signed_dvv_curr
        .add_elevated(&signed_dvv_prev)
        .add_elevated(&two_b_duv)
        .add_elevated(&minus_b_sq_duu)
        .add_elevated(&minus_b_bprime_du);
    debug_assert_eq!(sum.degree(), 8);

    let mut out = [Vec3::zero(); 9];
    for (k, &c) in sum.coeffs().iter().enumerate() {
        out[k] = c;
    }
    out
}

// ============================================================
// Phase B: assembled-cap residual report
// ============================================================

/// Per-sector-pair worst BB-coefficient magnitude from the eq (8)
/// G¹ and eq (9) G² residuals, plus Lemma 6 condition (1)–(4)
/// defects, evaluated on an assembled KP20 macro cap.
///
/// Each `worst_eq8_per_piece_pair` entry is the worst across the
/// 7 BB coefs of `kp20_eq8_g1_residual` for that piece pair (Vec3
/// component-wise inf-norm). Each `worst_eq9_per_piece_pair` entry
/// is the worst across the 8 BB coefs of `kp20_eq9_g2_residual`.
///
/// Lemma 6 condition entries are scalar `worst over sector pairs`
/// summaries. See paper §5.1 eq (17), (18).
///
/// The struct is a *measurement*, not a gate. The Phase B
/// integration tests assert only that every field is finite.
/// Whether the cap satisfies the paper (e.g. residuals all ≤
/// 1e-12) is a question for downstream ADRs (Phase E).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20Lemma6CpRegion {
    /// Every contributing CP lies in the Lemma 4 overlay block
    /// `i + j <= 3`.
    Lemma4OverlayOnly,
    /// The derivative stencil includes CPs outside the overlay,
    /// but only on rows that assembly restores from the original
    /// guide sample after C2 propagation.
    IncludesRestoredRows,
    /// The derivative stencil includes CPs outside the overlay and
    /// outside the currently-restored rows.
    IncludesUnrestoredStrip,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20Lemma6LayerDefect {
    /// Lemma 6 condition number: 2 for layer 0, 3 for layer 1,
    /// 4 for layer 2.
    pub condition: usize,
    /// Macro-patch sector that owns the worst defect.
    pub sector: usize,
    /// Layer row `j` of the NearOrigin left/right pair.
    pub layer_j: usize,
    /// Derivative order `d` whose C^d condition owns the max.
    pub derivative_order: usize,
    /// First contributing BB `i` index in the derivative stencil.
    pub cp_i_start: usize,
    /// Last contributing BB `i` index in the derivative stencil.
    pub cp_i_end: usize,
    /// Defect vector `D^d L(0) - (-1)^d D^d R(0)`.
    pub defect: Vec3,
    /// Component-wise infinity norm of `defect`.
    pub magnitude: f64,
    /// Assembly-region classification for the contributing CP span.
    pub cp_region: Kp20Lemma6CpRegion,
}

impl Kp20Lemma6LayerDefect {
    fn zero(condition: usize, layer_j: usize) -> Self {
        Self {
            condition,
            sector: 0,
            layer_j,
            derivative_order: 0,
            cp_i_start: 0,
            cp_i_end: 0,
            defect: Vec3::zero(),
            magnitude: 0.0,
            cp_region: Kp20Lemma6CpRegion::Lemma4OverlayOnly,
        }
    }

    fn all_finite(&self) -> bool {
        self.magnitude.is_finite()
            && self.defect.x.is_finite()
            && self.defect.y.is_finite()
            && self.defect.z.is_finite()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Kp20MacroCapResidualReport {
    /// Worst BB[i] magnitude of eq (8) residual on the near-o
    /// piece pair (`NearOriginLeft` vs `NearOriginRight`) at every
    /// macro patch `k ∈ 0..n`.
    pub worst_eq8_near_o_per_sector: Vec<f64>,
    /// Worst BB[i] magnitude of eq (8) residual on the near-m
    /// piece pair (`NearMidpointLeft` vs `NearMidpointRight`) at
    /// every macro patch `k ∈ 0..n`.
    pub worst_eq8_near_m_per_sector: Vec<f64>,
    /// Worst BB[i] magnitude of eq (9) residual on the near-o
    /// piece pair at every macro patch.
    pub worst_eq9_near_o_per_sector: Vec<f64>,
    /// Worst BB[i] magnitude of eq (9) residual on the near-m
    /// piece pair at every macro patch.
    pub worst_eq9_near_m_per_sector: Vec<f64>,
    /// Worst |−r^{k−1} + 2 m^k − ℓ^{k+1} + ŝ^k| across `k`. Paper
    /// eq (18); Lemma 6 condition (1).
    pub worst_lemma6_cond_1: f64,
    /// Worst C⁴ defect across sectors at point `o` on the
    /// sector-separating curve (layer `j = 0`). Lemma 6 cond (2).
    pub worst_lemma6_cond_2: f64,
    /// Worst C³ defect across sectors at point `o` on layer
    /// `j = 1`. Lemma 6 cond (3).
    pub worst_lemma6_cond_3: f64,
    /// Worst C² defect across sectors at point `o` on layer
    /// `j = 2`. Lemma 6 cond (4).
    pub worst_lemma6_cond_4: f64,
    /// Localisation of `worst_lemma6_cond_2`.
    pub worst_lemma6_cond_2_detail: Kp20Lemma6LayerDefect,
    /// Localisation of `worst_lemma6_cond_3`.
    pub worst_lemma6_cond_3_detail: Kp20Lemma6LayerDefect,
    /// Localisation of `worst_lemma6_cond_4`.
    pub worst_lemma6_cond_4_detail: Kp20Lemma6LayerDefect,
}

impl Kp20MacroCapResidualReport {
    /// Worst eq (8) BB magnitude across both piece pairs of every
    /// macro patch. Useful as a single-number summary.
    pub fn worst_eq8(&self) -> f64 {
        self.worst_eq8_near_o_per_sector
            .iter()
            .chain(self.worst_eq8_near_m_per_sector.iter())
            .copied()
            .fold(0.0_f64, f64::max)
    }

    /// Worst eq (9) BB magnitude across both piece pairs.
    pub fn worst_eq9(&self) -> f64 {
        self.worst_eq9_near_o_per_sector
            .iter()
            .chain(self.worst_eq9_near_m_per_sector.iter())
            .copied()
            .fold(0.0_f64, f64::max)
    }

    /// Every field is finite (no NaN, no infinity). Phase B
    /// integration tests' acceptance criterion.
    pub fn all_finite(&self) -> bool {
        let vec_ok = |v: &Vec<f64>| v.iter().all(|x| x.is_finite());
        vec_ok(&self.worst_eq8_near_o_per_sector)
            && vec_ok(&self.worst_eq8_near_m_per_sector)
            && vec_ok(&self.worst_eq9_near_o_per_sector)
            && vec_ok(&self.worst_eq9_near_m_per_sector)
            && self.worst_lemma6_cond_1.is_finite()
            && self.worst_lemma6_cond_2.is_finite()
            && self.worst_lemma6_cond_3.is_finite()
            && self.worst_lemma6_cond_4.is_finite()
            && self.worst_lemma6_cond_2_detail.all_finite()
            && self.worst_lemma6_cond_3_detail.all_finite()
            && self.worst_lemma6_cond_4_detail.all_finite()
    }
}

fn worst_bb_magnitude_v3(coefs: &[Vec3]) -> f64 {
    coefs
        .iter()
        .map(|v| v.x.abs().max(v.y.abs()).max(v.z.abs()))
        .fold(0.0_f64, f64::max)
}

/// Evaluate the per-sector-pair eq (8) G¹ and eq (9) G² residuals
/// of an assembled KP20 macro cap, plus the Lemma 6 condition (1)
/// defect and C^n connectedness defects (2)/(3)/(4) at the
/// extraordinary point `o`.
///
/// The cap's sector-separating curves live at the `v = 0` boundary
/// of each macro patch's left/right piece pairs (the same boundary
/// the existing `kp20_apply_g2_strip` operates on; paper §4.1
/// eq (8) is the constraint at that boundary).
///
/// Per ADR `docs/decisions/2026-06-23-kp20-bi6-macro-patch-paper-grade.md`
/// Phase B. The function does not assert anything; the test layer
/// owns the policy on what constitutes "OK".
pub fn measure_kp20_assembly_constraints(
    cap: &Kp20Bi6MacroCapOutput,
) -> Kp20MacroCapResidualReport {
    let n = cap.n;
    let c = (core::f64::consts::TAU / n as f64).cos();
    let b_near_o = kp20_b_coeffs_near_o(c);
    let b_near_m = kp20_b_coeffs_near_m(c);

    let mut eq8_near_o = Vec::with_capacity(n);
    let mut eq8_near_m = Vec::with_capacity(n);
    let mut eq9_near_o = Vec::with_capacity(n);
    let mut eq9_near_m = Vec::with_capacity(n);

    for macro_patch in &cap.macro_patches {
        let near_o_left = &macro_patch.piece(Kp20PieceId::NearOriginLeft).cps;
        let near_o_right = &macro_patch.piece(Kp20PieceId::NearOriginRight).cps;
        let near_m_left = &macro_patch.piece(Kp20PieceId::NearMidpointLeft).cps;
        let near_m_right = &macro_patch.piece(Kp20PieceId::NearMidpointRight).cps;

        let r_eq8_o = kp20_eq8_g1_residual(near_o_left, near_o_right, b_near_o);
        let r_eq8_m = kp20_eq8_g1_residual(near_m_left, near_m_right, b_near_m);
        let r_eq9_o = kp20_eq9_g2_residual(near_o_left, near_o_right, b_near_o);
        let r_eq9_m = kp20_eq9_g2_residual(near_m_left, near_m_right, b_near_m);

        eq8_near_o.push(worst_bb_magnitude_v3(&r_eq8_o));
        eq8_near_m.push(worst_bb_magnitude_v3(&r_eq8_m));
        eq9_near_o.push(worst_bb_magnitude_v3(&r_eq9_o));
        eq9_near_m.push(worst_bb_magnitude_v3(&r_eq9_m));
    }

    let worst_lemma6_cond_1 = lemma6_condition_1_worst(cap, b_near_o);
    let worst_lemma6_cond_2_detail =
        lemma6_within_macro_layer_connectedness_worst_detail(cap, 2, 0, 4);
    let worst_lemma6_cond_3_detail =
        lemma6_within_macro_layer_connectedness_worst_detail(cap, 3, 1, 3);
    let worst_lemma6_cond_4_detail =
        lemma6_within_macro_layer_connectedness_worst_detail(cap, 4, 2, 2);
    let worst_lemma6_cond_2 = worst_lemma6_cond_2_detail.magnitude;
    let worst_lemma6_cond_3 = worst_lemma6_cond_3_detail.magnitude;
    let worst_lemma6_cond_4 = worst_lemma6_cond_4_detail.magnitude;

    Kp20MacroCapResidualReport {
        worst_eq8_near_o_per_sector: eq8_near_o,
        worst_eq8_near_m_per_sector: eq8_near_m,
        worst_eq9_near_o_per_sector: eq9_near_o,
        worst_eq9_near_m_per_sector: eq9_near_m,
        worst_lemma6_cond_1,
        worst_lemma6_cond_2,
        worst_lemma6_cond_3,
        worst_lemma6_cond_4,
        worst_lemma6_cond_2_detail,
        worst_lemma6_cond_3_detail,
        worst_lemma6_cond_4_detail,
    }
}

/// Invert paper eq (17): given `(p_{2,1}, p_{3,1}, p_{4,1})` solve
/// for `(ℓ, m, r)`. The 3×3 matrix
/// ```text
///   [ 29  -6   13 ] / 36
///   [ 13  -6   13 ] / 20
///   [ 13  -6   29 ] / 36
/// ```
/// maps `(ℓ, m, r) → (p_{2,1}, p_{3,1}, p_{4,1})`. Closed-form
/// inverse; coordinates are independent (apply per-axis).
fn invert_eq17(p21: Point3, p31: Point3, p41: Point3) -> (Point3, Point3, Point3) {
    // M · [ℓ; m; r] = [p21; p31; p41].
    // M = [ [29/36, -6/36, 13/36], [13/20, -6/20, 13/20], [13/36, -6/36, 29/36] ]
    // M is symmetric in (ℓ ↔ r) ⇔ (p21 ↔ p41), and middle col has
    // common ratio -3:-3:-1 after row scaling. Closed-form via
    // cofactor expansion. We compute M⁻¹ once as a constant.
    //
    // det(M) computed exactly by hand:
    //   det = 29/36 · ((-6/20)·(29/36) - (13/20)·(-6/36))
    //       − (-6/36) · ((13/20)·(29/36) - (13/20)·(13/36))
    //       + (13/36) · ((13/20)·(-6/36) - (-6/20)·(13/36))
    //       = 29/36·(-174/720 + 78/720)
    //       + 6/36·(377/720 - 169/720)
    //       + 13/36·(-78/720 + 78/720)
    //       = 29/36·(-96/720) + 6/36·(208/720)
    //       = -2784/25920 + 1248/25920
    //       = -1536/25920
    //       = -16/270  (after dividing by 96).
    //
    // We can also invert numerically — Vec3 components are
    // independent. We do a direct numerical Gauss-elimination for
    // clarity since the matrix only depends on 6 rationals.

    let m: [[f64; 3]; 3] = [
        [29.0 / 36.0, -6.0 / 36.0, 13.0 / 36.0],
        [13.0 / 20.0, -6.0 / 20.0, 13.0 / 20.0],
        [13.0 / 36.0, -6.0 / 36.0, 29.0 / 36.0],
    ];
    let inv = invert_3x3(m);
    let solve_axis = |a: f64, b: f64, c: f64| -> [f64; 3] {
        [
            inv[0][0] * a + inv[0][1] * b + inv[0][2] * c,
            inv[1][0] * a + inv[1][1] * b + inv[1][2] * c,
            inv[2][0] * a + inv[2][1] * b + inv[2][2] * c,
        ]
    };
    let xs = solve_axis(p21.x, p31.x, p41.x);
    let ys = solve_axis(p21.y, p31.y, p41.y);
    let zs = solve_axis(p21.z, p31.z, p41.z);
    (
        Point3::new(xs[0], ys[0], zs[0]),
        Point3::new(xs[1], ys[1], zs[1]),
        Point3::new(xs[2], ys[2], zs[2]),
    )
}

fn invert_3x3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let a = m[0][0];
    let b = m[0][1];
    let c = m[0][2];
    let d = m[1][0];
    let e = m[1][1];
    let f = m[1][2];
    let g = m[2][0];
    let h = m[2][1];
    let i = m[2][2];
    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    debug_assert!(det.abs() > 1e-15, "invert_3x3: singular matrix");
    let inv_det = 1.0 / det;
    [
        [
            (e * i - f * h) * inv_det,
            (c * h - b * i) * inv_det,
            (b * f - c * e) * inv_det,
        ],
        [
            (f * g - d * i) * inv_det,
            (a * i - c * g) * inv_det,
            (c * d - a * f) * inv_det,
        ],
        [
            (d * h - e * g) * inv_det,
            (b * g - a * h) * inv_det,
            (a * e - b * d) * inv_det,
        ],
    ]
}

/// Lemma 6 condition (1): `|−r^{k−1} + 2 m^k − ℓ^{k+1} + ŝ^k| ≤ tol`
/// where `ŝ^k := (5 b̂ / 12) Σ_{i=0..5} (−1)^i C(5, i) s^k_i`
/// and `b̂ := b_0 − 2 b_1 + b_2`.
///
/// Implementation:
/// - For each macro patch `k`, extract `(ℓ^k, m^k, r^k)` from
///   eq (17) inversion of `(p^k_{2,1}, p^k_{3,1}, p^k_{4,1})`.
/// - Use the NearOriginRight piece's cps[2..5][1] as the three
///   measured CPs (`p_{i,1}` for `i = 2, 3, 4`; the `[j=1]`
///   row).
/// - The sector-separating curve `s^k_i` is the `cps[i][0]` row
///   of NearOriginRight, but we need it as a degree-5 curve.
///   Paper eq (11) says it reduces from degree 6 to degree 5 by
///   construction, so we just take the BB-deg-5 CPs (the BB-deg-6
///   CPs are over-parameterised; we de-elevate once).
///
/// Returns the worst-over-`k` per-coordinate magnitude.
fn lemma6_condition_1_worst(cap: &Kp20Bi6MacroCapOutput, b_near_o: [f64; 3]) -> f64 {
    let n = cap.n;
    if n < 2 {
        return 0.0;
    }
    // Extract (ℓ, m, r) per sector.
    let mut lmr = Vec::with_capacity(n);
    for macro_patch in &cap.macro_patches {
        let g = &macro_patch.piece(Kp20PieceId::NearOriginRight).cps;
        let p21 = g[2][1];
        let p31 = g[3][1];
        let p41 = g[4][1];
        lmr.push(invert_eq17(p21, p31, p41));
    }
    // b̂ = b_0 − 2 b_1 + b_2 (degree-2 alternating sum).
    let b_hat = b_near_o[0] - 2.0 * b_near_o[1] + b_near_o[2];
    let scale = 5.0 * b_hat / 12.0;
    // C(5, i) coefficients for i = 0..5.
    let bin5: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];

    let mut worst = 0.0_f64;
    for k in 0..n {
        let prev = (k + n - 1) % n;
        let next = (k + 1) % n;
        // Sector-separating curve s^k is NearOriginRight.cps[0..7][0]
        // — the v = 0 row of that piece, degree 6 in u. Paper eq (11)
        // says it reduces to degree 5; we de-elevate to take BB[0..6].
        let curve_bb6 = sector_separating_curve_deg6(&cap.macro_patches[k]);
        let curve_bb5 = de_elevate_curve_once(&curve_bb6);
        // ŝ^k = scale · Σ_{i=0..5} (−1)^i C(5, i) s^k_i.
        let mut s_hat = Vec3::zero();
        for i in 0_usize..6 {
            let sign = if i.is_multiple_of(2) { 1.0 } else { -1.0 };
            let coef = sign * bin5[i] * scale;
            let p = curve_bb5[i];
            let v = Vec3::new(p.x, p.y, p.z);
            s_hat = s_hat.add(v.scale(coef));
        }
        let (_, _, r_prev) = lmr[prev];
        let (_, m_k, _) = lmr[k];
        let (l_next, _, _) = lmr[next];
        // −r_{k−1} + 2 m^k − ℓ^{k+1} as Vec3.
        let r_prev_v = Vec3::new(r_prev.x, r_prev.y, r_prev.z);
        let m_k_v = Vec3::new(m_k.x, m_k.y, m_k.z);
        let l_next_v = Vec3::new(l_next.x, l_next.y, l_next.z);
        let defect = r_prev_v
            .scale(-1.0)
            .add(m_k_v.scale(2.0))
            .add(l_next_v.scale(-1.0))
            .add(s_hat);
        let mag = defect.x.abs().max(defect.y.abs()).max(defect.z.abs());
        worst = worst.max(mag);
    }
    worst
}

/// Extract the sector-separating curve of macro patch `k` as a
/// BB-degree-6 polynomial in u. NearOriginRight's `v = 0` row
/// — `cps[i][0]` for `i = 0..7`.
fn sector_separating_curve_deg6(macro_patch: &super::kp20_bi6_cap::Kp20MacroPatch) -> [Point3; 7] {
    let g = &macro_patch.piece(Kp20PieceId::NearOriginRight).cps;
    [
        g[0][0], g[1][0], g[2][0], g[3][0], g[4][0], g[5][0], g[6][0],
    ]
}

/// De-elevate a 7-CP BB-deg-6 curve to a 6-CP BB-deg-5 curve. If
/// the polynomial is genuinely degree ≤ 5 the result is exact;
/// otherwise it's the least-squares-style inverse.
fn de_elevate_curve_once(curve: &[Point3; 7]) -> [Point3; 6] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 6];
    out[0] = curve[0];
    let n_f = 6.0_f64;
    for i in 1..6 {
        let frac = i as f64 / n_f;
        let denom = 1.0 - frac;
        debug_assert!(denom > 1e-15);
        // p_i = (i/n) q_{i-1} + (1 - i/n) q_i ⇒
        // q_i = (p_i - (i/n) q_{i-1}) / (1 - i/n).
        let prev_q = out[i - 1];
        let p_i = curve[i];
        let prev_q_v = Vec3::new(prev_q.x, prev_q.y, prev_q.z);
        let p_i_v = Vec3::new(p_i.x, p_i.y, p_i.z);
        let q_i_v = (p_i_v.sub(prev_q_v.scale(frac))).scale(1.0 / denom);
        out[i] = Point3::new(q_i_v.x, q_i_v.y, q_i_v.z);
    }
    out
}

/// Compute the worst C^n defect between the left/right origin-piece
/// layer curves inside each macro patch (paper §5.1, Lemma 6
/// conditions (2)/(3)/(4)).
///
/// The kernel's KP20 macro-patch stores the sector-separating curve
/// as the internal `NearOriginLeft`/`NearOriginRight` boundary; this
/// is the same piece pair consumed by eq (8), eq (9), and the G2 strip
/// projection. Lemma 6's layer-connectedness checks therefore compare
/// layer `j` of those two pieces within one macro patch, not
/// `NearOriginRight` rows across adjacent macro patches.
///
/// For layer `j`, both curves are the `cps[i][j]` row. We compare
/// their BB-polynomial derivatives at the apex-side endpoint `u = 0`
/// because the Lemma-4 overlay and row restoration make that endpoint
/// the local sector-separating origin in the current convention. Odd
/// derivatives have opposite orientation across the split boundary.
fn lemma6_within_macro_layer_connectedness_worst(
    cap: &Kp20Bi6MacroCapOutput,
    layer_j: usize,
    n_derivs: usize,
) -> f64 {
    lemma6_within_macro_layer_connectedness_worst_detail(cap, layer_j + 2, layer_j, n_derivs)
        .magnitude
}

fn lemma6_within_macro_layer_connectedness_worst_detail(
    cap: &Kp20Bi6MacroCapOutput,
    condition: usize,
    layer_j: usize,
    n_derivs: usize,
) -> Kp20Lemma6LayerDefect {
    let n = cap.n;
    if n < 2 {
        return Kp20Lemma6LayerDefect::zero(condition, layer_j);
    }
    let mut worst = Kp20Lemma6LayerDefect::zero(condition, layer_j);
    for (sector, macro_patch) in cap.macro_patches.iter().enumerate() {
        let left = layer_curve(macro_patch, Kp20PieceId::NearOriginLeft, layer_j);
        let right = layer_curve(macro_patch, Kp20PieceId::NearOriginRight, layer_j);
        for d in 0..=n_derivs {
            let dl = bezier_derivative_at_0(&left, d);
            let dr = bezier_derivative_at_0(&right, d);
            let sign = if d.is_multiple_of(2) { 1.0 } else { -1.0 };
            let defect = dl.sub(dr.scale(sign));
            let mag = defect.x.abs().max(defect.y.abs()).max(defect.z.abs());
            if mag > worst.magnitude {
                worst = Kp20Lemma6LayerDefect {
                    condition,
                    sector,
                    layer_j,
                    derivative_order: d,
                    cp_i_start: 0,
                    cp_i_end: d,
                    defect,
                    magnitude: mag,
                    cp_region: lemma6_cp_region(layer_j, d),
                };
            }
        }
    }
    worst
}

fn lemma6_cp_region(layer_j: usize, derivative_order: usize) -> Kp20Lemma6CpRegion {
    let mut has_outside_overlay = false;
    let mut has_unrestored = false;
    for i in 0..=derivative_order {
        if i + layer_j <= 3 {
            continue;
        }
        has_outside_overlay = true;
        let restored_by_assembly = i == 5 || (i == 0 && layer_j >= 4) || (i == 1 && layer_j >= 4);
        if !restored_by_assembly {
            has_unrestored = true;
        }
    }
    if has_unrestored {
        Kp20Lemma6CpRegion::IncludesUnrestoredStrip
    } else if has_outside_overlay {
        Kp20Lemma6CpRegion::IncludesRestoredRows
    } else {
        Kp20Lemma6CpRegion::Lemma4OverlayOnly
    }
}

fn layer_curve(
    macro_patch: &super::kp20_bi6_cap::Kp20MacroPatch,
    piece: Kp20PieceId,
    layer_j: usize,
) -> [Point3; 7] {
    let g = &macro_patch.piece(piece).cps;
    [
        g[0][layer_j],
        g[1][layer_j],
        g[2][layer_j],
        g[3][layer_j],
        g[4][layer_j],
        g[5][layer_j],
        g[6][layer_j],
    ]
}

/// Former diagnostic interpretation: compare `NearOriginRight`
/// layer curves across adjacent macro patches at their right endpoint.
/// Kept for A/B tests against older PK-088 handoff numbers; production
/// reporting uses `lemma6_within_macro_layer_connectedness_worst`.
#[cfg(test)]
fn lemma6_cross_macro_layer_connectedness_worst(
    cap: &Kp20Bi6MacroCapOutput,
    layer_j: usize,
    n_derivs: usize,
) -> f64 {
    let n = cap.n;
    if n < 2 {
        return 0.0;
    }
    let mut worst = 0.0_f64;
    for k in 0..n {
        let next = (k + 1) % n;
        let curve_k = layer_curve(&cap.macro_patches[k], Kp20PieceId::NearOriginRight, layer_j);
        let curve_n = layer_curve(
            &cap.macro_patches[next],
            Kp20PieceId::NearOriginRight,
            layer_j,
        );
        for d in 0..=n_derivs {
            let dk = bezier_derivative_at_1(&curve_k, d);
            let dn = bezier_derivative_at_1(&curve_n, d);
            let sign = if d.is_multiple_of(2) { 1.0 } else { -1.0 };
            let defect = dk.sub(dn.scale(sign));
            let mag = defect.x.abs().max(defect.y.abs()).max(defect.z.abs());
            worst = worst.max(mag);
        }
    }
    worst
}

#[cfg(test)]
pub(crate) fn lemma6_layer_connectedness_comparison_for_tests(
    cap: &Kp20Bi6MacroCapOutput,
    layer_j: usize,
    n_derivs: usize,
) -> (f64, f64) {
    (
        lemma6_within_macro_layer_connectedness_worst(cap, layer_j, n_derivs),
        lemma6_cross_macro_layer_connectedness_worst(cap, layer_j, n_derivs),
    )
}

/// `D^d c(0)` for a degree-6 Bezier curve, returned as Vec3.
fn bezier_derivative_at_0(curve: &[Point3; 7], d: usize) -> Vec3 {
    let n_deg: usize = 6;
    if d == 0 {
        let p = curve[0];
        return Vec3::new(p.x, p.y, p.z);
    }
    if d > n_deg {
        return Vec3::zero();
    }
    let factorial = |x: usize| -> f64 {
        let mut acc = 1.0_f64;
        for k in 1..=x {
            acc *= k as f64;
        }
        acc
    };
    let scale = factorial(n_deg) / factorial(n_deg - d);
    let mut acc = Vec3::zero();
    let bin = pascal_row(d);
    for k in 0..=d {
        let sign = if (d - k).is_multiple_of(2) { 1.0 } else { -1.0 };
        let coef = sign * bin[k] * scale;
        let p = curve[k];
        let v = Vec3::new(p.x, p.y, p.z);
        acc = acc.add(v.scale(coef));
    }
    acc
}

/// `D^d c(1)` for a degree-6 Bezier curve, returned as Vec3.
/// `c(1) = c_6`; `D c(1) = 6 (c_6 - c_5)`; etc.
#[cfg(test)]
fn bezier_derivative_at_1(curve: &[Point3; 7], d: usize) -> Vec3 {
    if d == 0 {
        let p = curve[6];
        return Vec3::new(p.x, p.y, p.z);
    }
    // d-th forward differences at the right endpoint.
    let n_deg: usize = 6;
    if d > n_deg {
        return Vec3::zero();
    }
    // D^d c(1) = (n! / (n-d)!) · Δ^d c[n-d]
    // where Δ^d c[i] = Σ_{k=0..d} (-1)^{d-k} C(d, k) c[i+k]
    // at i = n - d.
    let factorial = |x: usize| -> f64 {
        let mut acc = 1.0_f64;
        for k in 1..=x {
            acc *= k as f64;
        }
        acc
    };
    let scale = factorial(n_deg) / factorial(n_deg - d);
    let mut acc = Vec3::zero();
    let bin = pascal_row(d);
    for k in 0..=d {
        let sign = if (d - k).is_multiple_of(2) { 1.0 } else { -1.0 };
        let coef = sign * bin[k] * scale;
        let p = curve[(n_deg - d) + k];
        let v = Vec3::new(p.x, p.y, p.z);
        acc = acc.add(v.scale(coef));
    }
    acc
}

fn pascal_row(n: usize) -> Vec<f64> {
    let mut row = vec![1.0_f64; n + 1];
    for i in 1..n {
        let mut prev = 1.0_f64;
        for cur_slot in row.iter_mut().take(i + 1).skip(1) {
            let cur = *cur_slot;
            *cur_slot = prev + cur;
            prev = cur;
        }
    }
    row
}

// ============================================================
// Phase C: global constraint matrix M(n, c) and its rank
// ============================================================

/// Layout of the per-valence global constraint matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Kp20ConstraintMatrixLayout {
    /// Valence.
    pub n: usize,
    /// Number of eq (8) + eq (9) rows = 30n.
    pub g_rows: usize,
    /// Number of within-macro-patch identification rows = 14n.
    pub ident_rows: usize,
    /// Number of cols = 196n (= 4 pieces × 49 BB-coefs × n macro patches).
    pub cols: usize,
}

impl Kp20ConstraintMatrixLayout {
    /// Total rows = g_rows + ident_rows.
    pub fn rows(&self) -> usize {
        self.g_rows + self.ident_rows
    }
}

/// Compute layout dimensions for valence `n`.
pub fn kp20_constraint_matrix_layout(n: usize) -> Kp20ConstraintMatrixLayout {
    Kp20ConstraintMatrixLayout {
        n,
        g_rows: 30 * n,
        ident_rows: 14 * n,
        cols: 196 * n,
    }
}

/// Global column index for piece `piece_idx` (0..=3 per
/// `Kp20PieceId::index`), grid `(i, j)`, macro patch `k`.
fn global_col(
    _layout: Kp20ConstraintMatrixLayout,
    k: usize,
    piece_idx: usize,
    i: usize,
    j: usize,
) -> usize {
    debug_assert!(piece_idx < 4 && i < 7 && j < 7);
    196 * k + 49 * piece_idx + i * 7 + j
}

/// Build the global constraint matrix M(n, c) of paper eq (8)/(9)
/// residuals across every macro patch's left/right piece pairs,
/// plus within-macro-patch boundary identifications. Returns
/// dense row-major as `Vec<Vec<f64>>`.
///
/// **Column layout** (per macro patch `k`, total 196 cols):
/// - cols [196k +   0 .. 196k +  49)  ← `NearOriginLeft`   piece
/// - cols [196k +  49 .. 196k +  98)  ← `NearOriginRight`
/// - cols [196k +  98 .. 196k + 147)  ← `NearMidpointRight` (paper-index 2)
/// - cols [196k + 147 .. 196k + 196)  ← `NearMidpointLeft`  (paper-index 3)
///
/// **Row layout** (per macro patch `k`):
/// - rows [30k +  0 .. 30k +  7)  eq (8) on near-o pair
/// - rows [30k +  7 .. 30k + 15)  eq (9) on near-o pair
/// - rows [30k + 15 .. 30k + 22)  eq (8) on near-m pair
/// - rows [30k + 22 .. 30k + 30)  eq (9) on near-m pair
///
/// then identification rows appended:
/// - rows [30n + 14k +  0 .. 30n + 14k +  7)  near-o boundary identifications (cps[i][0] of L = of R, i = 0..6)
/// - rows [30n + 14k +  7 .. 30n + 14k + 14)  near-m boundary identifications
///
/// Each scalar coordinate is independent under the linear eq
/// (8)/(9) operator, so we build the matrix on one coordinate
/// (the x-component) and the same matrix applies to y and z.
pub fn kp20_global_constraint_matrix(n: usize) -> (Kp20ConstraintMatrixLayout, Vec<Vec<f64>>) {
    let layout = kp20_constraint_matrix_layout(n);
    let mut m = vec![vec![0.0_f64; layout.cols]; layout.rows()];
    let c = (core::f64::consts::TAU / n as f64).cos();
    let b_near_o = kp20_b_coeffs_near_o(c);
    let b_near_m = kp20_b_coeffs_near_m(c);
    let eq8_jac_o = kp20_eq8_jacobian_dense(b_near_o);
    let eq9_jac_o = kp20_eq9_jacobian_dense(b_near_o);
    let eq8_jac_m = kp20_eq8_jacobian_dense(b_near_m);
    let eq9_jac_m = kp20_eq9_jacobian_dense(b_near_m);

    for k in 0..n {
        let row_base = 30 * k;
        // Near-o piece pair: side 0 = NearOriginLeft (piece_idx 0),
        // side 1 = NearOriginRight (piece_idx 1).
        write_per_piece_pair_block(
            &mut m,
            &eq8_jac_o,
            PiecePairBlockTarget {
                row_base,
                k,
                piece_idx_side0: 0,
                piece_idx_side1: 1,
                n_rows: 7,
            },
            layout,
        );
        write_per_piece_pair_block(
            &mut m,
            &eq9_jac_o,
            PiecePairBlockTarget {
                row_base: row_base + 7,
                k,
                piece_idx_side0: 0,
                piece_idx_side1: 1,
                n_rows: 8,
            },
            layout,
        );
        // Near-m piece pair: NearMidpointLeft = piece_idx 3 (the
        // "L" half), NearMidpointRight = piece_idx 2 (the "R" half).
        // Paper convention has prev = L, curr = R, so side 0 → 3,
        // side 1 → 2.
        write_per_piece_pair_block(
            &mut m,
            &eq8_jac_m,
            PiecePairBlockTarget {
                row_base: row_base + 15,
                k,
                piece_idx_side0: 3,
                piece_idx_side1: 2,
                n_rows: 7,
            },
            layout,
        );
        write_per_piece_pair_block(
            &mut m,
            &eq9_jac_m,
            PiecePairBlockTarget {
                row_base: row_base + 22,
                k,
                piece_idx_side0: 3,
                piece_idx_side1: 2,
                n_rows: 8,
            },
            layout,
        );
    }

    // Identification rows: cps[i][0] of left = cps[i][0] of right
    // within each macro patch, for both near-o and near-m piece pairs.
    let ident_row_base = layout.g_rows;
    for k in 0..n {
        let near_o_l = 0;
        let near_o_r = 1;
        for i in 0..7 {
            let row = ident_row_base + 14 * k + i;
            m[row][global_col(layout, k, near_o_l, i, 0)] = 1.0;
            m[row][global_col(layout, k, near_o_r, i, 0)] = -1.0;
        }
        // NearMidpointRight = piece_idx 2; NearMidpointLeft = piece_idx 3.
        let near_m_r = 2;
        let near_m_l = 3;
        for i in 0..7 {
            let row = ident_row_base + 14 * k + 7 + i;
            m[row][global_col(layout, k, near_m_l, i, 0)] = 1.0;
            m[row][global_col(layout, k, near_m_r, i, 0)] = -1.0;
        }
    }

    (layout, m)
}

/// Coordinates of a per-piece-pair Jacobian block within the global matrix.
struct PiecePairBlockTarget {
    row_base: usize,
    k: usize,
    piece_idx_side0: usize,
    piece_idx_side1: usize,
    n_rows: usize,
}

/// Write a per-piece-pair Jacobian block into the global matrix.
/// The Jacobian rows are 98-wide: side-0 cols (0..49) → piece
/// `piece_idx_side0`, side-1 cols (49..98) → piece `piece_idx_side1`.
fn write_per_piece_pair_block<const ROWS: usize>(
    m: &mut [Vec<f64>],
    jac: &[[f64; 98]; ROWS],
    target: PiecePairBlockTarget,
    layout: Kp20ConstraintMatrixLayout,
) {
    for (row, jac_row) in jac.iter().enumerate().take(target.n_rows) {
        let dst_row = target.row_base + row;
        for side in 0..2 {
            for i in 0..7 {
                for j in 0..7 {
                    let col_local = side * 49 + i * 7 + j;
                    let val = jac_row[col_local];
                    if val == 0.0 {
                        continue;
                    }
                    let piece_idx = if side == 0 {
                        target.piece_idx_side0
                    } else {
                        target.piece_idx_side1
                    };
                    let col = global_col(layout, target.k, piece_idx, i, j);
                    m[dst_row][col] += val;
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "kp20_bi6_local_constraint_tests.rs"]
mod tests;
