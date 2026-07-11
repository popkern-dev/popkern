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

use super::direct_linalg::solve_dense;
use super::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::bernstein_poly::{
    BernsteinPoly, PolyCoeff, bernstein_product_scalar, bernstein_product_scalar_vec3,
};
use crate::geometry::{Point3, Vec3};

/// One polar-cap sector's bi-6 patch as its 7×7 Bernstein–Bézier
/// control grid, indexed `[i][j]` per KP15's `p_{i,j}` convention.
pub type Bi6SectorGrid = [[Point3; 7]; 7];

/// Per-sector cap grids plus the achieved per-coordinate F_κ
/// functional values `(fx, fy, fz)` — the common return shape of
/// the F_κ population entry points.
pub type CapGridsWithFKappa = Result<(Vec<Bi6SectorGrid>, [f64; 3]), String>;

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

/// Global BB-coef column index for sector `k`, grid position
/// `(i, j)`. Sector wraps modulo `n`.
pub fn global_grid_col(n: usize, k: usize, i: usize, j: usize) -> usize {
    debug_assert!(k < n && i < 7 && j < 7);
    k * 49 + i * 7 + j
}

/// Append [KP15 §3.2 eq (5)] per-sector constraints to a global
/// matrix `m` of `(rows × cols)`. Rule (5) reads:
///
/// ```text
///   p̀_{30} := (3/4)(p̀_{20} + p̀_{40})
///            − (3/10)(p̀_{10} + p̀_{50})
///            + (1/20)(p̀_{00} + p̀_{60})
/// ```
///
/// rearranged to a homogeneous form. One row per sector is
/// appended.
fn append_rule5_rows(n: usize, cols: usize, m: &mut Vec<f64>, new_rows: usize) -> usize {
    // Reserve n new rows.
    let extra = vec![0.0_f64; n * cols];
    m.extend(extra);
    // Coefficients for `p̀_{0,0}, p̀_{1,0}, ..., p̀_{6,0}`.
    let weights: [f64; 7] = [
        -1.0 / 20.0,
        3.0 / 10.0,
        -3.0 / 4.0,
        1.0,
        -3.0 / 4.0,
        3.0 / 10.0,
        -1.0 / 20.0,
    ];
    for k in 0..n {
        let row_idx = new_rows + k;
        for (i, &w) in weights.iter().enumerate() {
            let col = global_grid_col(n, k, i, 0);
            m[row_idx * cols + col] = w;
        }
    }
    new_rows + n
}

/// Stack the per-edge eq (2) + eq (3) Jacobians into the global
/// constraint matrix `M(n, c, τ)` for the full ring of `n`
/// cross-sector edges.
///
/// Layout:
/// - **Columns** (`49n` total): one per BB-coef of the bi-6 grids,
///   indexed by `global_grid_col(n, k, i, j)`.
/// - **Rows** (`18n` total): for each edge `(k − 1 mod n, k)`,
///   8 rows from eq (2) followed by 10 rows from eq (3).
///
/// At adjacent edge `(k − 1, k)`, side 0 of the per-edge Jacobian
/// maps to global column `(k − 1) · 49 + i · 7 + j`, side 1 maps
/// to `k · 49 + i · 7 + j`. So each sector's 49 BB-coefs appear
/// in two consecutive edges' constraints — exactly the cross-edge
/// coupling the paper's local system encodes.
///
/// Returns the global matrix as `(n_rows, n_cols, entries_row_major)`
/// so the caller can apply standard linear algebra.
pub fn global_constraint_matrix(n: usize, params: &Bi6Reparams) -> (usize, usize, Vec<f64>) {
    let eq2 = eq2_jacobian_dense(params);
    let eq3 = eq3_jacobian_dense(params);
    let rows_per_edge = 8 + eq3.len();
    let n_rows = n * rows_per_edge;
    let n_cols = 49 * n;
    let mut m = vec![0.0_f64; n_rows * n_cols];
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let row_base = k * rows_per_edge;
        // eq (2) rows.
        for (r, eq2_row) in eq2.iter().enumerate() {
            let row_idx = row_base + r;
            for (col_local, &v) in eq2_row.iter().enumerate() {
                if v == 0.0 {
                    continue;
                }
                let side = col_local / 49;
                let local_pos = col_local % 49;
                let i = local_pos / 7;
                let j = local_pos % 7;
                let sector = if side == 0 { km1 } else { k };
                let global_col = global_grid_col(n, sector, i, j);
                m[row_idx * n_cols + global_col] += v;
            }
        }
        // eq (3) rows.
        for (r, eq3_row) in eq3.iter().enumerate() {
            let row_idx = row_base + 8 + r;
            for (col_local, &v) in eq3_row.iter().enumerate() {
                if v == 0.0 {
                    continue;
                }
                let side = col_local / 49;
                let local_pos = col_local % 49;
                let i = local_pos / 7;
                let j = local_pos % 7;
                let sector = if side == 0 { km1 } else { k };
                let global_col = global_grid_col(n, sector, i, j);
                m[row_idx * n_cols + global_col] += v;
            }
        }
    }
    (n_rows, n_cols, m)
}

/// Rule (5)-augmented constraint matrix: eq (2) + eq (3) per edge
/// plus one per-sector rule (5) row. Rows: `18n + n = 19n`.
pub fn global_constraint_matrix_with_rule5(
    n: usize,
    params: &Bi6Reparams,
) -> (usize, usize, Vec<f64>) {
    let (rows, cols, mut m) = global_constraint_matrix(n, params);
    let new_rows = append_rule5_rows(n, cols, &mut m, rows);
    (new_rows, cols, m)
}

/// Append cross-sector identifications `p^k_{i, 0} = p^{k − 1}_{i, 0}`
/// for `i = 0, 1, 2` and `k ∈ 0..n` (mod n). These come from
/// [KP15 eq (6)]'s (S1) re-affirmation that the BB-coefs of the
/// boundary curve near the irregular point are shared between
/// adjacent sectors.
fn append_cross_sector_identifications(
    n: usize,
    cols: usize,
    m: &mut Vec<f64>,
    base: usize,
) -> usize {
    // 3 constraints per edge × n edges. The full circular system
    // has rank 3·(n-1) — one rotation is redundant.
    let new_count = 3 * n;
    m.extend(vec![0.0_f64; new_count * cols]);
    let mut row = base;
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        for i in 0..3 {
            let col_curr = global_grid_col(n, k, i, 0);
            let col_prev = global_grid_col(n, km1, i, 0);
            m[row * cols + col_curr] = 1.0;
            m[row * cols + col_prev] = -1.0;
            row += 1;
        }
    }
    row
}

/// Constraint matrix with cross-sector identifications appended
/// to the per-edge eq (2)/(3) Jacobians. Rows: `18n + 3n = 21n`.
pub fn global_constraint_matrix_with_identifications(
    n: usize,
    params: &Bi6Reparams,
) -> (usize, usize, Vec<f64>) {
    let (rows, cols, mut m) = global_constraint_matrix(n, params);
    let new_rows = append_cross_sector_identifications(n, cols, &mut m, rows);
    (new_rows, cols, m)
}

/// Constraint matrix with the previous side of each local KP15 edge
/// transposed into the global grid before adding the S1 shared-edge
/// identifications. This matches the local two-sided notation where the
/// two adjacent patches see the same sector-separating curve with
/// opposite tensor axes.
pub fn global_constraint_matrix_with_transposed_previous_side(
    n: usize,
    params: &Bi6Reparams,
) -> (usize, usize, Vec<f64>) {
    let eq2 = eq2_jacobian_dense(params);
    let eq3 = eq3_jacobian_dense(params);
    let rows_per_edge = 8 + eq3.len();
    let base_rows = n * rows_per_edge;
    let rows = base_rows + 3 * n;
    let cols = 49 * n;
    let mut m = vec![0.0_f64; rows * cols];
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let row_base = k * rows_per_edge;
        for (r, eq2_row) in eq2.iter().enumerate() {
            append_transposed_previous_side_edge_row(&mut m, cols, row_base + r, km1, k, eq2_row);
        }
        for (r, eq3_row) in eq3.iter().enumerate() {
            append_transposed_previous_side_edge_row(
                &mut m,
                cols,
                row_base + 8 + r,
                km1,
                k,
                eq3_row,
            );
        }
        for i in 0..3 {
            let row = base_rows + 3 * k + i;
            m[row * cols + global_grid_col(n, km1, 0, i)] = 1.0;
            m[row * cols + global_grid_col(n, k, i, 0)] = -1.0;
        }
    }
    (rows, cols, m)
}

/// Constraint matrix with the transposed previous-side edge rows and
/// the full [KP15 §3.2 (S1)] quadratic-near-star source rows.
pub fn global_constraint_matrix_with_transposed_previous_side_and_full_s1(
    n: usize,
    params: &Bi6Reparams,
) -> (usize, usize, Vec<f64>) {
    global_constraint_matrix_with_transposed_previous_side_assignment_prefix(n, params, 6)
}

pub fn global_constraint_matrix_with_transposed_previous_side_assignment_prefix(
    n: usize,
    params: &Bi6Reparams,
    assignment_rows_per_edge: usize,
) -> (usize, usize, Vec<f64>) {
    assert!((6..=18).contains(&assignment_rows_per_edge));
    let eq2 = eq2_jacobian_dense(params);
    let eq3 = eq3_jacobian_dense(params);
    let rows_per_edge = 8 + eq3.len();
    let base_rows = n * rows_per_edge;
    let rows = base_rows + assignment_rows_per_edge * n;
    let cols = 49 * n;
    let mut m = vec![0.0_f64; rows * cols];
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let row_base = k * rows_per_edge;
        for (r, eq2_row) in eq2.iter().enumerate() {
            append_transposed_previous_side_edge_row(&mut m, cols, row_base + r, km1, k, eq2_row);
        }
        for (r, eq3_row) in eq3.iter().enumerate() {
            append_transposed_previous_side_edge_row(
                &mut m,
                cols,
                row_base + 8 + r,
                km1,
                k,
                eq3_row,
            );
        }
        append_kp15_assignment_prefix_rows(
            &mut m,
            cols,
            base_rows + assignment_rows_per_edge * k,
            km1,
            k,
            params,
            assignment_rows_per_edge,
        );
    }
    (rows, cols, m)
}

pub fn global_constraint_matrix_with_transposed_previous_side_full_edge_assignment_prefix(
    n: usize,
    params: &Bi6Reparams,
    assignment_rows_per_edge: usize,
) -> (usize, usize, Vec<f64>) {
    assert!((6..=18).contains(&assignment_rows_per_edge));
    let eq2 = eq2_jacobian_dense(params);
    let eq3 = eq3_jacobian_dense(params);
    let rows_per_edge = 8 + eq3.len();
    let base_rows = n * rows_per_edge;
    let shared_edge_rows = 7 * n;
    let assignment_base = base_rows + shared_edge_rows;
    let rows = assignment_base + assignment_rows_per_edge * n;
    let cols = 49 * n;
    let mut m = vec![0.0_f64; rows * cols];
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let row_base = k * rows_per_edge;
        for (r, eq2_row) in eq2.iter().enumerate() {
            append_transposed_previous_side_edge_row(&mut m, cols, row_base + r, km1, k, eq2_row);
        }
        for (r, eq3_row) in eq3.iter().enumerate() {
            append_transposed_previous_side_edge_row(
                &mut m,
                cols,
                row_base + 8 + r,
                km1,
                k,
                eq3_row,
            );
        }
        append_full_shared_edge_rows(&mut m, cols, base_rows + 7 * k, km1, k);
        append_kp15_assignment_prefix_rows(
            &mut m,
            cols,
            assignment_base + assignment_rows_per_edge * k,
            km1,
            k,
            params,
            assignment_rows_per_edge,
        );
    }
    (rows, cols, m)
}

/// Constraint matrix with transposed previous-side edge rows and the
/// explicit KP15 §3.2 S1-S4 coefficient assignment rows.
pub fn global_constraint_matrix_with_transposed_previous_side_assignments(
    n: usize,
    params: &Bi6Reparams,
) -> (usize, usize, Vec<f64>) {
    global_constraint_matrix_with_transposed_previous_side_assignment_prefix(n, params, 18)
}

/// Matrix containing only the explicit KP15 §3.2 S1-S4 assignment
/// rows around the ring.
pub fn global_constraint_matrix_with_assignment_rows_only(
    n: usize,
    params: &Bi6Reparams,
) -> (usize, usize, Vec<f64>) {
    let assignment_rows_per_edge = 18;
    let rows = assignment_rows_per_edge * n;
    let cols = 49 * n;
    let mut m = vec![0.0_f64; rows * cols];
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        append_kp15_assignment_rows(&mut m, cols, assignment_rows_per_edge * k, km1, k, params);
    }
    (rows, cols, m)
}

fn append_full_shared_edge_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
) {
    let n = cols / 49;
    for i in 0..=6 {
        let row = &mut m[(base_row + i) * cols..(base_row + i + 1) * cols];
        add_curr_side_coeff(row, n, current_sector, i, 0, 1.0);
        add_prev_side_coeff(row, n, previous_sector, i, 0, -1.0);
    }
}

fn append_kp15_assignment_prefix_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
    params: &Bi6Reparams,
    row_count: usize,
) {
    assert!((6..=18).contains(&row_count));
    append_full_s1_rows(m, cols, base_row, previous_sector, current_sector, params);
    if row_count == 6 {
        return;
    }
    assert!(row_count >= 9);
    append_kp15_s2_rows(
        m,
        cols,
        base_row + 6,
        previous_sector,
        current_sector,
        params,
    );
    if row_count == 9 {
        return;
    }
    assert!(row_count >= 12);
    append_kp15_s3_rows(
        m,
        cols,
        base_row + 9,
        previous_sector,
        current_sector,
        params,
    );
    if row_count == 12 {
        return;
    }
    assert_eq!(row_count, 18);
    append_kp15_s4_rows(
        m,
        cols,
        base_row + 12,
        previous_sector,
        current_sector,
        params,
    );
}

fn add_prev_side_coeff(
    row: &mut [f64],
    n: usize,
    previous_sector: usize,
    local_i: usize,
    local_j: usize,
    value: f64,
) {
    row[global_grid_col(n, previous_sector, local_j, local_i)] += value;
}

fn add_curr_side_coeff(
    row: &mut [f64],
    n: usize,
    current_sector: usize,
    local_i: usize,
    local_j: usize,
    value: f64,
) {
    row[global_grid_col(n, current_sector, local_i, local_j)] += value;
}

fn append_full_s1_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
    params: &Bi6Reparams,
) {
    let n = cols / 49;
    let c = params.c;
    for i in 0..3 {
        let row = &mut m[(base_row + i) * cols..(base_row + i + 1) * cols];
        add_curr_side_coeff(row, n, current_sector, i, 0, 1.0);
        add_prev_side_coeff(row, n, previous_sector, i, 0, -1.0);
    }

    let row = &mut m[(base_row + 3) * cols..(base_row + 4) * cols];
    add_curr_side_coeff(row, n, current_sector, 0, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 0, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 1, 0, -2.0 * c);
    add_prev_side_coeff(row, n, previous_sector, 0, 0, -2.0 * (1.0 - c));

    let row = &mut m[(base_row + 4) * cols..(base_row + 5) * cols];
    add_curr_side_coeff(row, n, current_sector, 1, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 1, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 2, 0, -5.0 * c / 3.0);
    add_prev_side_coeff(row, n, previous_sector, 1, 0, -2.0 * (1.0 - c));
    add_prev_side_coeff(row, n, previous_sector, 0, 0, -c / 3.0);

    let tau = params.tau;
    let row = &mut m[(base_row + 5) * cols..(base_row + 6) * cols];
    add_curr_side_coeff(row, n, current_sector, 0, 2, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 0, 2, -1.0);
    add_prev_side_coeff(row, n, previous_sector, 1, 1, 24.0 * c / 5.0);
    add_prev_side_coeff(row, n, previous_sector, 2, 0, -4.0 * c * c);
    add_prev_side_coeff(
        row,
        n,
        previous_sector,
        0,
        1,
        -4.0 * (6.0 * tau * c / 5.0 - 1.0),
    );
    add_prev_side_coeff(
        row,
        n,
        previous_sector,
        1,
        0,
        -4.0 * c / 5.0 * (11.0 - 6.0 * (1.0 + tau) * c),
    );
    add_prev_side_coeff(
        row,
        n,
        previous_sector,
        0,
        0,
        -4.0 * (1.0 - c - 6.0 * tau * c / 5.0 + (6.0 * tau + 1.0) * c * c / 5.0),
    );
}

fn append_kp15_assignment_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
    params: &Bi6Reparams,
) {
    append_full_s1_rows(m, cols, base_row, previous_sector, current_sector, params);
    let mut row_offset = base_row + 6;
    append_kp15_s2_rows(m, cols, row_offset, previous_sector, current_sector, params);
    row_offset += 3;
    append_kp15_s3_rows(m, cols, row_offset, previous_sector, current_sector, params);
    row_offset += 3;
    append_kp15_s4_rows(m, cols, row_offset, previous_sector, current_sector, params);
}

fn append_kp15_s2_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
    params: &Bi6Reparams,
) {
    let n = cols / 49;
    let c = params.c;
    let tau = params.tau;
    let dot = [2.0 * c / 5.0, 2.0 - c, c, -2.0 * c / 5.0, c / 15.0];
    let ddot = [
        -4.0 * c * tau * tau / 5.0,
        4.0 * (1.0 - c + c * c / 5.0)
            + 24.0 * c * (c - 1.0) * tau / 5.0
            + 4.0 * c * c * tau * tau / 5.0,
        24.0 * c * tau / 5.0 - 4.0,
        2.0 * c * (11.0 - 3.0 * c - 6.0 * c * tau) / 3.0,
        -4.0 * c,
        2.0 * c * c,
        -4.0 * c * c / 5.0,
        2.0 * c * c / 15.0,
    ];
    let tdot = [
        4.0 * c * tau * tau * tau / 25.0,
        4.0 * c * (5.0 - 3.0 * c * tau + 6.0 * (1.0 - c) * tau * tau - c * tau * tau * tau) / 25.0,
        -24.0 * c * tau * tau / 25.0,
        4.0 - 2.0 * c + 6.0 * c * (c - 2.0) * tau / 5.0 + 4.0 * c * c * tau * tau / 5.0,
        12.0 * c * tau / 5.0 - 4.0,
        2.0 * c * (5.0 - 3.0 * c * tau) / 5.0,
        4.0 * c * (3.0 * c * tau - 5.0) / 25.0,
        2.0 * c * (5.0 - 3.0 * c * tau) / 75.0,
        -3.0 * c / (5.0 * tau),
        6.0 * c * (5.0 - c * tau) / (125.0 * tau * tau),
        c * (6.0 * c * tau - 5.0) / (125.0 * tau * tau * tau),
    ];

    let row = &mut m[base_row * cols..(base_row + 1) * cols];
    add_curr_side_coeff(row, n, current_sector, 2, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 2, 1, 1.0);
    add_prev_side_coeff(
        row,
        n,
        previous_sector,
        0,
        0,
        -(2.0 - dot.iter().sum::<f64>()),
    );
    for ((i, j), value) in [(1, 0), (2, 0), (4, 0), (5, 0), (6, 0)]
        .into_iter()
        .zip(dot)
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
    }

    let row = &mut m[(base_row + 1) * cols..(base_row + 2) * cols];
    add_curr_side_coeff(row, n, current_sector, 1, 2, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 1, 2, -1.0);
    add_prev_side_coeff(row, n, previous_sector, 0, 0, ddot.iter().sum::<f64>());
    for ((i, j), value) in [
        (0, 1),
        (1, 0),
        (1, 1),
        (2, 0),
        (2, 1),
        (4, 0),
        (5, 0),
        (6, 0),
    ]
    .into_iter()
    .zip(ddot)
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
    }

    let row = &mut m[(base_row + 2) * cols..(base_row + 3) * cols];
    add_curr_side_coeff(row, n, current_sector, 2, 2, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 2, 2, -1.0);
    add_prev_side_coeff(
        row,
        n,
        previous_sector,
        0,
        0,
        tdot[0..8].iter().sum::<f64>(),
    );
    for ((i, j), value) in [
        (0, 1),
        (1, 0),
        (1, 1),
        (2, 0),
        (2, 1),
        (4, 0),
        (5, 0),
        (6, 0),
    ]
    .into_iter()
    .zip(tdot[0..8].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
    }
    for ((i, j), value) in [(4, 2), (5, 2), (6, 2)]
        .into_iter()
        .zip(tdot[8..11].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
        add_curr_side_coeff(row, n, current_sector, i, j, value);
    }
}

fn append_kp15_s3_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
    params: &Bi6Reparams,
) {
    let n = cols / 49;
    let c = params.c;
    let tau = params.tau;
    let tmu = [
        -11.0 / 20.0,
        (11.0 * c - 10.0) / (8.0 * c),
        5.0 / (4.0 * c),
        1.0 / 8.0 - c / 4.0,
        3.0 * c / 10.0 - 1.0 / 20.0,
        1.0 / 120.0 - c / 20.0,
        3.0 / (8.0 * tau),
        (3.0 * c * tau - 15.0) / (100.0 * tau * tau),
        (5.0 - 6.0 * c * tau) / (200.0 * tau * tau * tau),
        -5.0 / (16.0 * c),
    ];
    let bmu = [
        1.0,
        5.0 * (2.0 - c) / (2.0 * c),
        -5.0 / c,
        5.0 / 2.0,
        -1.0,
        1.0 / 6.0,
        3.0 * (c * tau - 5.0) / (10.0 * tau),
        3.0 * (25.0 - 20.0 * c * tau + c * c * tau * tau) / (125.0 * tau * tau),
        -(25.0 - 45.0 * c * tau + 6.0 * c * c * tau * tau) / (250.0 * tau * tau * tau),
        5.0 / (4.0 * c),
    ];
    let base_coeff = 1.0 - tmu[0..6].iter().sum::<f64>();

    let row = &mut m[base_row * cols..(base_row + 1) * cols];
    add_prev_side_coeff(row, n, previous_sector, 3, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 0, 0, -base_coeff);
    for ((i, j), value) in [(1, 0), (2, 0), (2, 1), (4, 0), (5, 0), (6, 0)]
        .into_iter()
        .zip(tmu[0..6].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
    }
    for ((i, j), value) in [(4, 2), (5, 2), (6, 2), (2, 2)]
        .into_iter()
        .zip(tmu[6..10].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
        add_curr_side_coeff(row, n, current_sector, i, j, value);
    }

    let row = &mut m[(base_row + 1) * cols..(base_row + 2) * cols];
    add_curr_side_coeff(row, n, current_sector, 3, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 0, 0, -base_coeff);
    for ((i, j), value) in [(1, 0), (2, 0), (2, 1), (4, 0), (5, 0), (6, 0)]
        .into_iter()
        .zip(tmu[0..6].iter().copied())
    {
        let source_is_current = (i, j) == (2, 1);
        if source_is_current {
            add_curr_side_coeff(row, n, current_sector, i, j, -value);
        } else {
            add_prev_side_coeff(row, n, previous_sector, i, j, -value);
        }
    }
    for ((i, j), value) in [(4, 2), (5, 2), (6, 2), (2, 2)]
        .into_iter()
        .zip(tmu[6..10].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, value);
        add_curr_side_coeff(row, n, current_sector, i, j, -value);
    }

    let row = &mut m[(base_row + 2) * cols..(base_row + 3) * cols];
    add_curr_side_coeff(row, n, current_sector, 3, 2, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 3, 2, -1.0);
    add_prev_side_coeff(row, n, previous_sector, 0, 0, bmu[0..6].iter().sum::<f64>());
    for ((i, j), value) in [(1, 0), (2, 0), (2, 1), (4, 0), (5, 0), (6, 0)]
        .into_iter()
        .zip(bmu[0..6].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
    }
    for ((i, j), value) in [(4, 2), (5, 2), (6, 2), (2, 2)]
        .into_iter()
        .zip(bmu[6..10].iter().copied())
    {
        add_prev_side_coeff(row, n, previous_sector, i, j, -value);
        add_curr_side_coeff(row, n, current_sector, i, j, value);
    }
}

fn append_kp15_s4_rows(
    m: &mut [f64],
    cols: usize,
    base_row: usize,
    previous_sector: usize,
    current_sector: usize,
    params: &Bi6Reparams,
) {
    let n = cols / 49;
    for (offset, i) in [5, 6].into_iter().enumerate() {
        let row = &mut m[(base_row + 2 * offset) * cols..(base_row + 2 * offset + 1) * cols];
        add_prev_side_coeff(row, n, previous_sector, i, 1, 1.0);
        add_prev_side_coeff(row, n, previous_sector, i, 0, -1.0);
        add_prev_side_coeff(row, n, previous_sector, i, 2, -0.25);
        add_curr_side_coeff(row, n, current_sector, i, 2, 0.25);

        let row = &mut m[(base_row + 2 * offset + 1) * cols..(base_row + 2 * offset + 2) * cols];
        add_curr_side_coeff(row, n, current_sector, i, 1, 1.0);
        add_prev_side_coeff(row, n, previous_sector, i, 0, -1.0);
        add_prev_side_coeff(row, n, previous_sector, i, 2, 0.25);
        add_curr_side_coeff(row, n, current_sector, i, 2, -0.25);
    }

    let c = params.c;
    let tau = params.tau;
    let w1 = -c / (50.0 * tau);
    let w2 = c / 50.0;
    let w3 = c / 15.0;
    let row = &mut m[(base_row + 4) * cols..(base_row + 5) * cols];
    add_prev_side_coeff(row, n, previous_sector, 4, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 4, 0, -1.0);
    add_prev_side_coeff(row, n, previous_sector, 4, 2, -0.25);
    add_curr_side_coeff(row, n, current_sector, 4, 2, 0.25);
    add_prev_side_coeff(row, n, previous_sector, 6, 2, -w1);
    add_curr_side_coeff(row, n, current_sector, 6, 2, w1);
    add_prev_side_coeff(row, n, previous_sector, 5, 2, -w2);
    add_curr_side_coeff(row, n, current_sector, 5, 2, w2);
    add_prev_side_coeff(row, n, previous_sector, 6, 0, -w3);
    add_prev_side_coeff(row, n, previous_sector, 5, 0, w3);

    let row = &mut m[(base_row + 5) * cols..(base_row + 6) * cols];
    add_curr_side_coeff(row, n, current_sector, 4, 1, 1.0);
    add_prev_side_coeff(row, n, previous_sector, 4, 0, -1.0);
    add_prev_side_coeff(row, n, previous_sector, 4, 2, 0.25);
    add_curr_side_coeff(row, n, current_sector, 4, 2, -0.25);
    add_prev_side_coeff(row, n, previous_sector, 6, 2, w1);
    add_curr_side_coeff(row, n, current_sector, 6, 2, -w1);
    add_prev_side_coeff(row, n, previous_sector, 5, 2, w2);
    add_curr_side_coeff(row, n, current_sector, 5, 2, -w2);
    add_prev_side_coeff(row, n, previous_sector, 6, 0, -w3);
    add_prev_side_coeff(row, n, previous_sector, 5, 0, w3);
}

fn append_transposed_previous_side_edge_row(
    m: &mut [f64],
    cols: usize,
    row_idx: usize,
    previous_sector: usize,
    current_sector: usize,
    local_row: &[f64; 98],
) {
    let n = cols / 49;
    for (col_local, &value) in local_row.iter().enumerate() {
        if value == 0.0 {
            continue;
        }
        let side = col_local / 49;
        let local_pos = col_local % 49;
        let i = local_pos / 7;
        let j = local_pos % 7;
        let (sector, mapped_i, mapped_j) = if side == 0 {
            (previous_sector, j, i)
        } else {
            (current_sector, i, j)
        };
        m[row_idx * cols + global_grid_col(n, sector, mapped_i, mapped_j)] += value;
    }
}

/// Reduced row-echelon form of a dense `(rows × cols)` matrix
/// stored row-major. Pivots are recorded as `(row, col)` pairs in
/// increasing column order. Pivot tolerance `tol` is the absolute
/// threshold below which a column is considered to have no pivot
/// available in the remaining rows.
///
/// Used by [`Kp15Bi6GeneratingTable::build`] below to split the
/// global constraint matrix's columns into "pivot" (dependent)
/// and "free" sets — the foundation of Phase D's solver.
fn rref_in_place(rows: usize, cols: usize, m: &mut [f64], tol: f64) -> Vec<(usize, usize)> {
    assert_eq!(m.len(), rows * cols);
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut row = 0usize;
    for col in 0..cols {
        if row >= rows {
            break;
        }
        // Partial pivoting.
        let mut best = row;
        let mut best_abs = m[row * cols + col].abs();
        for r in (row + 1)..rows {
            let v = m[r * cols + col].abs();
            if v > best_abs {
                best = r;
                best_abs = v;
            }
        }
        if best_abs < tol {
            continue;
        }
        if best != row {
            for c in 0..cols {
                m.swap(row * cols + c, best * cols + c);
            }
        }
        // Normalise pivot row.
        let pivot = m[row * cols + col];
        for c in col..cols {
            m[row * cols + c] /= pivot;
        }
        // Eliminate above and below.
        for r in 0..rows {
            if r == row {
                continue;
            }
            let factor = m[r * cols + col];
            if factor == 0.0 {
                continue;
            }
            for c in col..cols {
                m[r * cols + c] -= factor * m[row * cols + c];
            }
        }
        pivots.push((row, col));
        row += 1;
    }
    pivots
}

/// A pre-computed RREF of [`global_constraint_matrix`]. Encodes
/// each pivot (dependent) BB-coef as an affine function of the
/// free BB-coefs: for pivot row `r` with pivot column `p`,
///
/// ```text
///   x_p = − Σ_{c ∈ free_cols} rref[r][c] · x_c.
/// ```
///
/// This is the same `H(n)` matrix the ADR's Phase D names; we
/// store it in RREF form so the structure is explicit and the
/// pivot/free split is auditable per-row.
///
/// **NB on rank deficiency:** for `n = 6, params = Bi6Reparams::
/// for_valence(6)` the pre-rule-(5) measured rank is 86 (out of
/// `49 · 6 = 294` cols). That leaves 208 columns in `free_cols`
/// — more than the paper's `21n + |P| = 126 + 19 = 145` "true
/// free" entries. The extra `208 - 145 = 63` free columns
/// represent the *cross-edge coupling defect* surfaced in Phase
/// C v2's measurement; closing them requires Phase D's
/// rule-(5)-driven substitution and the (S1)–(S4) tightening.
/// For now, this table exposes the *raw* RREF so subsequent
/// phases can layer constraints atop it.
pub struct Kp15Bi6GeneratingTable {
    /// Valence.
    pub n: usize,
    /// Per-valence parameters at which `H(n)` was built.
    pub params: Bi6Reparams,
    /// `(row, col)` pairs of the pivots, in increasing column order.
    pub pivots: Vec<(usize, usize)>,
    /// Free-column indices (= columns not in `pivots`'s second slot).
    pub free_cols: Vec<usize>,
    /// Row count of the RREF (≡ `n × 18`).
    pub rows: usize,
    /// Column count (≡ `n × 49`).
    pub cols: usize,
    /// RREF entries, row-major.
    pub rref: Vec<f64>,
}

impl Kp15Bi6GeneratingTable {
    /// Build the table for valence `n` at parameters `params`.
    pub fn build(n: usize, params: Bi6Reparams) -> Self {
        Self::build_from_matrix(n, params, global_constraint_matrix)
    }

    /// Build with [KP15 eq (5)] applied as `n` extra per-sector
    /// rows on top of the per-edge eq (2)/(3) Jacobians.
    pub fn build_with_rule5(n: usize, params: Bi6Reparams) -> Self {
        Self::build_from_matrix(n, params, global_constraint_matrix_with_rule5)
    }

    /// Build with cross-sector `p^k_{i, 0} = p^{k − 1}_{i, 0}`
    /// identifications (i = 0, 1, 2; k ∈ 0..n) appended as `3n`
    /// extra rows on top of the per-edge eq (2)/(3) Jacobians.
    /// These are the (S1) corner-shared identifications the paper
    /// re-affirms at line 158-161; the rank lift measures whether
    /// they are already implied by eq (2)/(3) (the rule-(5)
    /// pattern from Phase D v2) or genuinely independent.
    pub fn build_with_identifications(n: usize, params: Bi6Reparams) -> Self {
        Self::build_from_matrix(n, params, global_constraint_matrix_with_identifications)
    }

    /// Build with cross-sector identifications, and bias the RREF
    /// so that [KP15 §3.6 step 4]'s "P" set ends up in the **free**
    /// (non-pivot) column set of the resulting table.
    ///
    /// "P" is the paper's notation for the cap's free-parameter
    /// set (sector 0's quadratic-near-star CPs + `p^0_{22}` +
    /// per-sector `p^k_{23}, p^k_{33}` + `p^0_{21}` for `n = 6`).
    /// See [`kp15_p_set_columns`] for the explicit column-index
    /// listing.
    ///
    /// Implication for downstream callers: the F_5 minimisation
    /// step (e.g.
    /// [`Self::populate_cap_grids_from_input_with_f_kappa`])
    /// classifies each free column as VARIABLE iff it isn't already
    /// seeded by the input. After this build, the VARIABLE set is
    /// largely the paper's P (minus the few overlaps with the 13
    /// seeded positions like the star and the row-0/1/2 boundary
    /// CPs), matching the paper's prescribed minimisation domain.
    /// Useful when downstream verification compares against the
    /// paper's published examples or asserts that the cap lives
    /// in the paper's exact solution space.
    ///
    /// Trade-off vs the default [`Self::build_with_identifications`]:
    /// fewer degrees of freedom for F_5 ⇒ higher achieved F_5
    /// (measured ~3.4× on the synthetic hexagon fixture). eq (2)/
    /// (3) residuals stay at machine zero either way; the lever is
    /// purely about *which* solution in the M·x = 0 nullspace
    /// we land on.
    pub fn build_with_kp15_p_set_as_free(n: usize, params: Bi6Reparams) -> Self {
        let (rows, cols, m_unordered) = global_constraint_matrix_with_identifications(n, &params);
        let p_set = kp15_p_set_columns(n);
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &p_set)
    }

    /// Build the paper-source diagnostic table: transpose the previous
    /// side of each local KP15 edge into the global grid, bias the
    /// tensor-border seed columns and central quadratic source columns
    /// into the free set, and leave the central quadratic expansion to
    /// the F_5 solve rather than treating every sector's near-origin
    /// quadratic coefficients as fixed tensor-border data.
    pub fn build_with_transposed_previous_side_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side(n, &params);
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Build with the full KP15 §3.2 (S1) central quadratic rows
    /// instead of only the three shared-edge identifications.
    pub fn build_with_transposed_previous_side_full_s1_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side_and_full_s1(n, &params);
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Diagnostic builder for transposed eq(2)/(3) plus KP15 §3.2
    /// S1-S2 assignment rows.
    pub fn build_with_transposed_previous_side_s1_s2_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side_assignment_prefix(n, &params, 9);
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Diagnostic builder for transposed eq(2)/(3) plus KP15 §3.2
    /// S1-S3 assignment rows.
    pub fn build_with_transposed_previous_side_s1_s3_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side_assignment_prefix(
                n, &params, 12,
            );
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Diagnostic builder for the explicit KP15 §3.2 S1-S4
    /// assignment rows, using the Appendix coefficients for
    /// equations (7), (8), and (9).
    pub fn build_with_transposed_previous_side_assignments_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side_assignments(n, &params);
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Build the Fig. 2b source-closed table: explicit KP15 §3.2
    /// S1-S4 rows plus the paper free set `P`, central quadratic
    /// source columns, and the full two-sided depth-2 tensor-border
    /// L-band biased into the free set.
    ///
    /// When callers fix the full two-sided tensor source and leave only
    /// `P` variable, the remaining row system should have zero nullity.
    pub fn build_with_transposed_previous_side_assignments_and_two_sided_tensor_source(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side_assignments(n, &params);
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.extend(kp15_two_sided_tensor_border_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Diagnostic builder for the explicit KP15 §3.2 S1-S4 rows
    /// plus full shared-curve identifications `p̌_{i0}=ṕ_{i0}`,
    /// `i=0..6`.
    pub fn build_with_transposed_previous_side_full_edge_assignments_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_transposed_previous_side_full_edge_assignment_prefix(
                n, &params, 18,
            );
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    /// Diagnostic builder for the KP15 §3.2 assignment rows without
    /// also appending the raw eq (2)/(3) coefficient rows.
    pub fn build_with_assignments_only_and_central_quadratic_free(
        n: usize,
        params: Bi6Reparams,
    ) -> Self {
        let (rows, cols, m_unordered) =
            global_constraint_matrix_with_assignment_rows_only(n, &params);
        let mut free_bias = kp15_p_set_columns(n);
        free_bias.extend(kp15_seed_columns(n));
        free_bias.extend(kp15_central_quadratic_columns(n));
        free_bias.sort_unstable();
        free_bias.dedup();
        Self::build_from_matrix_with_free_bias(n, params, rows, cols, m_unordered, &free_bias)
    }

    fn build_from_matrix_with_free_bias(
        n: usize,
        params: Bi6Reparams,
        rows: usize,
        cols: usize,
        m_unordered: Vec<f64>,
        free_bias_cols: &[usize],
    ) -> Self {
        let free_bias_lookup: std::collections::HashSet<usize> =
            free_bias_cols.iter().copied().collect();
        let mut perm: Vec<usize> = (0..cols)
            .filter(|c| !free_bias_lookup.contains(c))
            .collect();
        perm.extend(free_bias_cols.iter().copied());
        debug_assert_eq!(perm.len(), cols);
        let mut inverse = vec![0usize; cols];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            inverse[old_idx] = new_idx;
        }
        // Apply the permutation: m_permuted[r, new_c] = m_unordered[r, perm[new_c]].
        let mut m_permuted = vec![0.0_f64; rows * cols];
        for r in 0..rows {
            for (new_c, &old_c) in perm.iter().enumerate() {
                m_permuted[r * cols + new_c] = m_unordered[r * cols + old_c];
            }
        }
        let pivots_in_perm = rref_in_place(rows, cols, &mut m_permuted, 1e-12);
        // Translate pivot/RREF entries back to the unordered column
        // indexing.
        let pivots: Vec<(usize, usize)> =
            pivots_in_perm.iter().map(|&(r, c)| (r, perm[c])).collect();
        let pivot_cols: std::collections::BTreeSet<usize> =
            pivots.iter().map(|&(_, c)| c).collect();
        let free_cols: Vec<usize> = (0..cols).filter(|c| !pivot_cols.contains(c)).collect();
        // Translate the RREF entries: rref[r, c_original] =
        // m_permuted[r, inverse[c_original]].
        let mut rref = vec![0.0_f64; rows * cols];
        for r in 0..rows {
            for c_orig in 0..cols {
                rref[r * cols + c_orig] = m_permuted[r * cols + inverse[c_orig]];
            }
        }
        Self {
            n,
            params,
            pivots,
            free_cols,
            rows,
            cols,
            rref,
        }
    }

    fn build_from_matrix(
        n: usize,
        params: Bi6Reparams,
        builder: fn(usize, &Bi6Reparams) -> (usize, usize, Vec<f64>),
    ) -> Self {
        let (rows, cols, mut m) = builder(n, &params);
        let pivots = rref_in_place(rows, cols, &mut m, 1e-12);
        let pivot_cols: std::collections::BTreeSet<usize> =
            pivots.iter().map(|&(_, c)| c).collect();
        let free_cols: Vec<usize> = (0..cols).filter(|c| !pivot_cols.contains(c)).collect();
        Self {
            n,
            params,
            pivots,
            free_cols,
            rows,
            cols,
            rref: m,
        }
    }

    /// Rank of the global constraint matrix.
    pub fn rank(&self) -> usize {
        self.pivots.len()
    }

    /// Given values for every column of the global grid (one f64 per
    /// `(k, i, j)` position; use `global_grid_col` to compute the
    /// index), compute the residual norm `‖M · x‖`. Useful as a
    /// sanity check that an externally-constructed cap satisfies
    /// the eq (2) + eq (3) BB constraints.
    pub fn constraint_residual_norm(&self, x: &[f64]) -> f64 {
        assert_eq!(x.len(), self.cols);
        // We can compute M·x by re-deriving M, but the RREF stored
        // here suffices: row r of RREF reads
        //     x[col(pivot)] + Σ_free rref[r][col_free] · x[col_free] = 0,
        // i.e. row r's residual is the dot product of the row with x.
        let mut acc = 0.0_f64;
        for r in 0..self.pivots.len() {
            let row_dot: f64 = x
                .iter()
                .enumerate()
                .map(|(c, &xc)| self.rref[r * self.cols + c] * xc)
                .sum();
            acc += row_dot * row_dot;
        }
        acc.sqrt()
    }

    /// Returns true iff global column `col` is a free (non-pivot)
    /// column of this RREF.
    pub fn is_free(&self, col: usize) -> bool {
        self.free_cols.binary_search(&col).is_ok()
    }

    /// Apply the table: given values for every **free** column,
    /// compute the values of all pivot columns to satisfy M·x = 0.
    /// Returns the full vector indexed by global column.
    ///
    /// `free_values` must equal `free_cols.len()` in length and
    /// align with `free_cols`.
    pub fn apply(&self, free_values: &[f64]) -> Vec<f64> {
        assert_eq!(free_values.len(), self.free_cols.len());
        let mut x = vec![0.0_f64; self.cols];
        for (&fc, &v) in self.free_cols.iter().zip(free_values.iter()) {
            x[fc] = v;
        }
        // Pivot rows: x_pivot = − Σ_free rref[row][col_free] · x[col_free].
        for &(row, pcol) in &self.pivots {
            let mut sum = 0.0;
            for &fc in &self.free_cols {
                sum += self.rref[row * self.cols + fc] * x[fc];
            }
            x[pcol] = -sum;
        }
        x
    }

    /// Phase E.2 — populate the `n` per-sector bi-6 grids of a
    /// polar cap from a `Bi6CapInput` using this generating table.
    ///
    /// Strategy. For each of the three coordinate axes
    /// independently:
    /// 1. Build a `49n`-long "seed" vector. At each
    ///    `global_grid_col(n, k, i, j)` position, write whatever
    ///    the `Bi6CapInput` provides for `(k, i, j)`:
    ///    - `(0, 0)` → the Catmull-Clark star (shared across all
    ///      sectors).
    ///    - Other positions in `CirculantSectorRhs`'s fixed 12-
    ///      slot pattern → the corresponding field of
    ///      `input.sectors[k]`.
    ///    - Else → zero (the table will treat it as a free input
    ///      defaulted to zero).
    /// 2. Sub-sample the seed at the table's `free_cols` indices.
    /// 3. Apply the table to compute the pivot-column values
    ///    satisfying `M · x = 0`.
    /// 4. Reshape the resulting `49n`-long vector into `n` per-
    ///    sector grids `[[Point3; 7]; 7]`.
    ///
    /// **Note on the rank gap.** The current table has rank ≈
    /// `17n + 3(n-1)` against `49n` cols. After fixing the 12n + n
    /// = 13n caller-supplied positions and the `13(n-1) + 3n + ...`
    /// pivot positions, ~6n free positions remain unspecified
    /// and default to zero. Geometrically this means the resulting
    /// cap is *partial*: rows j = 3..6 (deep interior) get
    /// well-defined zero seeds, and the pivot positions are
    /// derived to keep `‖M · x‖ ≈ 0`. Use the residual reported
    /// by `constraint_residual_norm` on each coordinate to verify.
    ///
    /// The caller is responsible for the `star` argument — the
    /// `Bi6CapInput.c7` plus `c5s`/`c6s` go through
    /// `catmull_clark_star_point` in the existing assembly path.
    pub fn populate_cap_grids_from_input(
        &self,
        sectors: &[super::kp15_bi6_circulant::CirculantSectorRhs],
        star: Point3,
    ) -> Vec<[[Point3; 7]; 7]> {
        let n = self.n;
        assert_eq!(sectors.len(), n, "sectors length must equal valence");

        // Build per-sector (i, j) → Point3 lookup.
        let mut seed = vec![Point3::new(0.0, 0.0, 0.0); 49 * n];
        for k in 0..n {
            let s = &sectors[k];
            // Star at (0, 0).
            seed[global_grid_col(n, k, 0, 0)] = star;
            // The 12 CirculantSectorRhs positions.
            seed[global_grid_col(n, k, 1, 0)] = s.p_10;
            seed[global_grid_col(n, k, 2, 0)] = s.p_20;
            seed[global_grid_col(n, k, 4, 0)] = s.p_40;
            seed[global_grid_col(n, k, 5, 0)] = s.p_50;
            seed[global_grid_col(n, k, 6, 0)] = s.p_60;
            seed[global_grid_col(n, k, 0, 1)] = s.p_01;
            seed[global_grid_col(n, k, 1, 1)] = s.p_11;
            seed[global_grid_col(n, k, 0, 2)] = s.p_02;
            seed[global_grid_col(n, k, 4, 2)] = s.p_42;
            seed[global_grid_col(n, k, 5, 2)] = s.p_52;
            seed[global_grid_col(n, k, 6, 2)] = s.p_62;
        }

        // Apply table per coordinate axis.
        let coord = |pick: fn(Point3) -> f64| -> Vec<f64> {
            let free_vals: Vec<f64> = self.free_cols.iter().map(|&col| pick(seed[col])).collect();
            self.apply(&free_vals)
        };
        let xs = coord(|p| p.x);
        let ys = coord(|p| p.y);
        let zs = coord(|p| p.z);

        // Reshape into per-sector grids.
        let mut grids = Vec::with_capacity(n);
        for k in 0..n {
            let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
            for (i, row) in grid.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    let col = global_grid_col(n, k, i, j);
                    *cell = Point3::new(xs[col], ys[col], zs[col]);
                }
            }
            grids.push(grid);
        }
        grids
    }

    /// Build the table's implicit linear map `T` from the free-
    /// column value vector to the full `49n`-long column vector,
    /// as a dense `49n × |free_cols|` matrix in row-major order.
    /// Column `j` of T is the action of `apply` on the unit vector
    /// at `j`.
    fn build_linear_map(&self) -> Vec<f64> {
        let n_free = self.free_cols.len();
        let n_total = self.cols;
        let mut t = vec![0.0_f64; n_total * n_free];
        for (j, _fc) in self.free_cols.iter().enumerate() {
            // Apply on the j-th unit vector.
            let mut free = vec![0.0_f64; n_free];
            free[j] = 1.0;
            let x = self.apply(&free);
            for (row, &v) in x.iter().enumerate() {
                t[row * n_free + j] = v;
            }
        }
        t
    }

    /// Phase E.3 — F_5-minimised cap construction via the table.
    ///
    /// Like [`populate_cap_grids_from_input`], but instead of
    /// defaulting all "extra" free columns (those not seeded by
    /// the `Bi6CapInput`) to zero, this method solves three
    /// independent SPD systems (one per coordinate) to choose
    /// those values such that the per-patch [KP15 §3.6 step 4]
    /// quadratic functional `F_κ` is minimised.
    ///
    /// Inputs:
    /// - `sectors`, `star`: as in `populate_cap_grids_from_input`.
    /// - `q_local`: the `49 × 49` per-patch F_κ quadratic form
    ///   (from
    ///   [`super::kp15_bi6_f5::f_kappa_quadratic_form`]). The
    ///   full `49n × 49n` Q is block-diagonal with this `q_local`
    ///   replicated `n` times.
    ///
    /// Algorithm:
    /// 1. Build the table's linear map `T` (49n × |free|).
    /// 2. Form `H = T^T Q T` (the |free| × |free| projected
    ///    quadratic form on free-column values).
    /// 3. Partition free columns into *FIXED* (those seeded by
    ///    the input — the 12 CirculantSectorRhs slots + the star
    ///    at `(0, 0)` per sector) and *VARIABLE* (the rest).
    /// 4. Solve `H_vv · v = − H_vf · f` for the VARIABLE
    ///    values, where `H = [[H_ff, H_fv], [H_vf, H_vv]]` is
    ///    the partitioned form and `f` is the fixed-value
    ///    vector. Three independent solves, one per coordinate
    ///    (paper §3.7: "works for each coordinate separately").
    /// 5. Combine fixed + variable free-values, apply the table,
    ///    reshape into per-sector grids.
    ///
    /// `ridge_lambda` is an optional ridge regularisation parameter.
    /// The SPD solve replaces `H_vv` with `H_vv + λI`, biasing the
    /// solution toward the minimum-norm point in the constraint
    /// manifold's nullspace.
    ///
    /// Pass `0.0` for unregularised behaviour (matches the paper's
    /// "F_κ alone" objective — the deep-interior CPs are then
    /// unbounded in magnitude, only their smoothness is penalised).
    /// Pass a small positive value (e.g. `1e-3`) to bound the
    /// "extras" — free columns the rank-deficient RREF leaves over
    /// after seeding tensor border + paper's P. The regulariser
    /// is paper-orthogonal: it doesn't affect constraint
    /// satisfaction (`‖M · x‖` stays at machine zero) and is a
    /// no-op on the paper-P-only RREF (where extras default to
    /// zero anyway).
    ///
    /// Returns the per-sector grids and the achieved per-coord
    /// F_κ values `(fx, fy, fz)`.
    pub fn populate_cap_grids_from_input_with_f_kappa(
        &self,
        sectors: &[super::kp15_bi6_circulant::CirculantSectorRhs],
        star: Point3,
        q_local: &[[f64; 49]; 49],
        ridge_lambda: f64,
    ) -> CapGridsWithFKappa {
        let n = self.n;
        assert_eq!(sectors.len(), n);

        // Build full Q at 49n × 49n (block-diagonal of q_local).
        // We store it as a dense Vec<f64> row-major.
        let total = 49 * n;
        // Step 1: T (49n × |free|).
        let t = self.build_linear_map();
        let n_free = self.free_cols.len();

        // Step 2: H = T^T Q T. Exploit block-diagonal structure:
        // H = Σ_k T_k^T q_local T_k where T_k is rows
        // [49k .. 49(k+1)] of T.
        let mut h = vec![0.0_f64; n_free * n_free];
        for k in 0..n {
            let row_base = 49 * k;
            // q_local · T_k (49 × |free|).
            let mut qtk = vec![0.0_f64; 49 * n_free];
            for r in 0..49 {
                for c in 0..n_free {
                    let mut s = 0.0;
                    for inner in 0..49 {
                        s += q_local[r][inner] * t[(row_base + inner) * n_free + c];
                    }
                    qtk[r * n_free + c] = s;
                }
            }
            // H += T_k^T · qtk (|free| × |free|).
            for r in 0..n_free {
                for c in 0..n_free {
                    let mut s = 0.0;
                    for inner in 0..49 {
                        s += t[(row_base + inner) * n_free + r] * qtk[inner * n_free + c];
                    }
                    h[r * n_free + c] += s;
                }
            }
        }

        // Step 3: classify free cols → FIXED or VARIABLE indices
        // into the free-column ordering.
        let mut fixed_idx: Vec<usize> = Vec::new();
        let mut variable_idx: Vec<usize> = Vec::new();
        let mut seed_x: Vec<f64> = Vec::new();
        let mut seed_y: Vec<f64> = Vec::new();
        let mut seed_z: Vec<f64> = Vec::new();
        for (j, &fc) in self.free_cols.iter().enumerate() {
            // Decode (k, i, j_pos).
            let sector = fc / 49;
            let pos = fc % 49;
            let i = pos / 7;
            let j_pos = pos % 7;
            if let Some(p) = seeded_value(sectors, star, sector, i, j_pos) {
                fixed_idx.push(j);
                seed_x.push(p.x);
                seed_y.push(p.y);
                seed_z.push(p.z);
            } else {
                variable_idx.push(j);
            }
        }
        let n_fixed = fixed_idx.len();
        let n_var = variable_idx.len();

        // Step 4: per coord, solve H_vv · v = − H_vf · f.
        let pull_h = |rows: &[usize], cols: &[usize]| -> Vec<Vec<f64>> {
            rows.iter()
                .map(|&r| cols.iter().map(|&c| h[r * n_free + c]).collect())
                .collect()
        };
        let mut h_vv = pull_h(&variable_idx, &variable_idx);
        let h_vf = pull_h(&variable_idx, &fixed_idx);

        // Ridge regularisation: H_vv ← H_vv + λI. Pulls the SPD
        // solution toward the minimum-norm point in the nullspace.
        // Affects only variable-column magnitudes; constraint
        // satisfaction (`‖M · x‖`) is preserved because the table's
        // `apply` step enforces M·x = 0 by RREF construction
        // regardless of which free-value vector we plug in.
        if ridge_lambda != 0.0 {
            for (r, row) in h_vv.iter_mut().enumerate().take(n_var) {
                row[r] += ridge_lambda;
            }
        }

        let mut var_vals: [Vec<f64>; 3] = [vec![0.0; n_var], vec![0.0; n_var], vec![0.0; n_var]];
        for (d, seed) in [&seed_x, &seed_y, &seed_z].iter().enumerate() {
            // rhs = -H_vf · seed.
            let rhs: Vec<f64> = (0..n_var)
                .map(|r| {
                    let mut s = 0.0;
                    for c in 0..n_fixed {
                        s -= h_vf[r][c] * seed[c];
                    }
                    s
                })
                .collect();
            let solution = solve_dense(h_vv.clone(), rhs)
                .ok_or_else(|| format!("F_κ SPD solve singular at coord {d}"))?;
            var_vals[d] = solution;
        }

        // Step 5: assemble free-value vectors and apply.
        let assemble_free = |coord: usize, seed_coord: &[f64]| -> Vec<f64> {
            let mut free = vec![0.0_f64; n_free];
            for (k, &j) in fixed_idx.iter().enumerate() {
                free[j] = seed_coord[k];
            }
            for (k, &j) in variable_idx.iter().enumerate() {
                free[j] = var_vals[coord][k];
            }
            free
        };
        let xs = self.apply(&assemble_free(0, &seed_x));
        let ys = self.apply(&assemble_free(1, &seed_y));
        let zs = self.apply(&assemble_free(2, &seed_z));

        // Step 6: reshape into per-sector grids.
        let mut grids = Vec::with_capacity(n);
        for k in 0..n {
            let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
            for (i, row) in grid.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    let col = global_grid_col(n, k, i, j);
                    *cell = Point3::new(xs[col], ys[col], zs[col]);
                }
            }
            grids.push(grid);
        }

        // Step 7: compute achieved F_κ per coord.
        let f_coord = |x: &[f64]| -> f64 {
            let mut acc = 0.0;
            for k in 0..n {
                let row_base = 49 * k;
                for r in 0..49 {
                    for c in 0..49 {
                        acc += x[row_base + r] * q_local[r][c] * x[row_base + c];
                    }
                }
            }
            acc
        };
        let f_vals = [f_coord(&xs), f_coord(&ys), f_coord(&zs)];
        let _ = total;
        Ok((grids, f_vals))
    }

    /// F_κ table population for the paper-source split: source
    /// coefficients with `i + j <= 2` are central quadratic expansion
    /// variables owned by KP15 §3.2 S1/P, while the remaining
    /// `CirculantSectorRhs` slots are fixed tensor-border data.
    pub fn populate_cap_grids_from_input_with_f_kappa_central_quadratic_free(
        &self,
        sectors: &[super::kp15_bi6_circulant::CirculantSectorRhs],
        star: Point3,
        q_local: &[[f64; 49]; 49],
        ridge_lambda: f64,
    ) -> CapGridsWithFKappa {
        let n = self.n;
        assert_eq!(sectors.len(), n);
        let n_free = self.free_cols.len();
        let t = self.build_linear_map();
        let mut h = vec![0.0_f64; n_free * n_free];
        for k in 0..n {
            let row_base = 49 * k;
            let mut qtk = vec![0.0_f64; 49 * n_free];
            for r in 0..49 {
                for c in 0..n_free {
                    let mut s = 0.0;
                    for inner in 0..49 {
                        s += q_local[r][inner] * t[(row_base + inner) * n_free + c];
                    }
                    qtk[r * n_free + c] = s;
                }
            }
            for r in 0..n_free {
                for c in 0..n_free {
                    let mut s = 0.0;
                    for inner in 0..49 {
                        s += t[(row_base + inner) * n_free + r] * qtk[inner * n_free + c];
                    }
                    h[r * n_free + c] += s;
                }
            }
        }

        let mut fixed_idx = Vec::new();
        let mut variable_idx = Vec::new();
        let mut seed_x = Vec::new();
        let mut seed_y = Vec::new();
        let mut seed_z = Vec::new();
        for (j, &fc) in self.free_cols.iter().enumerate() {
            let sector = fc / 49;
            let pos = fc % 49;
            let i = pos / 7;
            let j_pos = pos % 7;
            if let Some(p) = seeded_value(sectors, star, sector, i, j_pos) {
                if i + j_pos > 2 {
                    fixed_idx.push(j);
                    seed_x.push(p.x);
                    seed_y.push(p.y);
                    seed_z.push(p.z);
                } else {
                    variable_idx.push(j);
                }
            } else {
                variable_idx.push(j);
            }
        }

        let pull_h = |rows: &[usize], cols: &[usize]| -> Vec<Vec<f64>> {
            rows.iter()
                .map(|&r| cols.iter().map(|&c| h[r * n_free + c]).collect())
                .collect()
        };
        let mut h_vv = pull_h(&variable_idx, &variable_idx);
        let h_vf = pull_h(&variable_idx, &fixed_idx);
        if ridge_lambda != 0.0 {
            for (r, row) in h_vv.iter_mut().enumerate() {
                row[r] += ridge_lambda;
            }
        }
        let mut var_vals: [Vec<f64>; 3] = [
            vec![0.0; variable_idx.len()],
            vec![0.0; variable_idx.len()],
            vec![0.0; variable_idx.len()],
        ];
        for (d, seed) in [&seed_x, &seed_y, &seed_z].iter().enumerate() {
            let rhs = (0..variable_idx.len())
                .map(|r| {
                    let mut s = 0.0;
                    for c in 0..fixed_idx.len() {
                        s -= h_vf[r][c] * seed[c];
                    }
                    s
                })
                .collect::<Vec<_>>();
            var_vals[d] = solve_dense(h_vv.clone(), rhs)
                .ok_or_else(|| format!("F_κ central-quadratic-free solve singular at coord {d}"))?;
        }

        let assemble_free = |coord: usize, seed_coord: &[f64]| -> Vec<f64> {
            let mut free = vec![0.0_f64; n_free];
            for (k, &j) in fixed_idx.iter().enumerate() {
                free[j] = seed_coord[k];
            }
            for (k, &j) in variable_idx.iter().enumerate() {
                free[j] = var_vals[coord][k];
            }
            free
        };
        let xs = self.apply(&assemble_free(0, &seed_x));
        let ys = self.apply(&assemble_free(1, &seed_y));
        let zs = self.apply(&assemble_free(2, &seed_z));

        let mut grids = Vec::with_capacity(n);
        for k in 0..n {
            let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
            for (i, row) in grid.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    let col = global_grid_col(n, k, i, j);
                    *cell = Point3::new(xs[col], ys[col], zs[col]);
                }
            }
            grids.push(grid);
        }
        let f_coord = |x: &[f64]| -> f64 {
            let mut acc = 0.0;
            for k in 0..n {
                let row_base = 49 * k;
                for r in 0..49 {
                    for c in 0..49 {
                        acc += x[row_base + r] * q_local[r][c] * x[row_base + c];
                    }
                }
            }
            acc
        };
        Ok((grids, [f_coord(&xs), f_coord(&ys), f_coord(&zs)]))
    }

    /// F_κ population for the source-closed KP15 Fig. 2b table.
    ///
    /// `tensor_source` must contain the full two-sided Fig. 12
    /// tensor-border source values in global sector grids. Only columns
    /// returned by [`kp15_two_sided_tensor_border_columns`] are consumed.
    /// Central quadratic and paper `P` columns remain variable, as
    /// required by KP15 §3.6.
    pub fn populate_cap_grids_from_two_sided_tensor_source_with_f_kappa(
        &self,
        tensor_source: &[[[Point3; 7]; 7]],
        q_local: &[[f64; 49]; 49],
        ridge_lambda: f64,
    ) -> CapGridsWithFKappa {
        let n = self.n;
        assert_eq!(tensor_source.len(), n);
        let source_constraints: Vec<(usize, Point3)> = kp15_two_sided_tensor_border_columns(n)
            .into_iter()
            .map(|col| {
                let sector = col / 49;
                let pos = col % 49;
                let i = pos / 7;
                let j = pos % 7;
                (col, tensor_source[sector][i][j])
            })
            .collect();
        self.populate_cap_grids_from_column_constraints_with_f_kappa(
            &source_constraints,
            q_local,
            ridge_lambda,
            "two-sided tensor-source",
        )
    }

    /// F_κ population subject to arbitrary paper-source column
    /// identities plus this table's KP15 row equations.
    ///
    /// This is the common substrate for source projection: the caller
    /// supplies only source values that are truly prescribed by the
    /// paper construction, and the remaining BB coefficients are chosen
    /// by the same constrained F_κ minimization.
    pub fn populate_cap_grids_from_column_constraints_with_f_kappa(
        &self,
        source_constraints: &[(usize, Point3)],
        q_local: &[[f64; 49]; 49],
        ridge_lambda: f64,
        source_label: &str,
    ) -> CapGridsWithFKappa {
        let n = self.n;
        let total = 49 * n;
        for &(col, _) in source_constraints {
            if col >= total {
                return Err(format!(
                    "F_κ {source_label} source column {col} outside 0..{total}"
                ));
            }
        }
        let mut q_global = vec![vec![0.0_f64; total]; total];
        for k in 0..n {
            let row_base = 49 * k;
            for r in 0..49 {
                for c in 0..49 {
                    q_global[row_base + r][row_base + c] = q_local[r][c];
                }
            }
        }
        if ridge_lambda != 0.0 {
            for (i, row) in q_global.iter_mut().enumerate() {
                row[i] += ridge_lambda;
            }
        }

        let mut constraint_rows = Vec::<Vec<f64>>::new();
        let mut rhs = [Vec::<f64>::new(), Vec::<f64>::new(), Vec::<f64>::new()];
        for row in 0..self.rank() {
            let mut constraint = vec![0.0_f64; total];
            let base = row * self.cols;
            constraint[..total].copy_from_slice(&self.rref[base..(base + total)]);
            constraint_rows.push(constraint);
            for rhs_coord in &mut rhs {
                rhs_coord.push(0.0);
            }
        }
        for &(col, source_point) in source_constraints {
            let mut constraint = vec![0.0_f64; total];
            constraint[col] = 1.0;
            constraint_rows.push(constraint);
            rhs[0].push(source_point.x);
            rhs[1].push(source_point.y);
            rhs[2].push(source_point.z);
        }

        let mut coords: [Vec<f64>; 3] = [
            vec![0.0_f64; total],
            vec![0.0_f64; total],
            vec![0.0_f64; total],
        ];
        for coord in 0..3 {
            let (active_constraints, active_rhs) =
                independent_affine_constraint_rows(total, &constraint_rows, &rhs[coord], 1.0e-10)?;
            let kkt_size = total + active_constraints.len();
            let mut kkt = vec![vec![0.0_f64; kkt_size]; kkt_size];
            for r in 0..total {
                for c in 0..total {
                    kkt[r][c] = q_global[r][c];
                }
            }
            for (constraint_index, constraint) in active_constraints.iter().enumerate() {
                let lagrange_col = total + constraint_index;
                for (col, &value) in constraint.iter().enumerate() {
                    kkt[col][lagrange_col] = value;
                    kkt[lagrange_col][col] = value;
                }
            }
            let mut kkt_rhs = vec![0.0_f64; kkt_size];
            kkt_rhs[total..].copy_from_slice(&active_rhs);
            let solution = solve_dense(kkt, kkt_rhs).ok_or_else(|| {
                format!("F_κ {source_label} constrained solve singular at coord {coord}")
            })?;
            coords[coord][..total].copy_from_slice(&solution[..total]);
        }

        let xs = &coords[0];
        let ys = &coords[1];
        let zs = &coords[2];

        let mut max_source_residual: f64 = 0.0;
        for &(col, source_point) in source_constraints {
            max_source_residual = max_source_residual.max((xs[col] - source_point.x).abs());
            max_source_residual = max_source_residual.max((ys[col] - source_point.y).abs());
            max_source_residual = max_source_residual.max((zs[col] - source_point.z).abs());
        }
        if max_source_residual > 1.0e-8 {
            return Err(format!(
                "F_κ {source_label} solve moved source constraints by {max_source_residual:e}"
            ));
        }

        let mut max_row_residual: f64 = 0.0;
        for row in 0..self.rank() {
            let base = row * self.cols;
            for coord_values in [&coords[0], &coords[1], &coords[2]] {
                let mut residual = 0.0;
                for (col, &value) in coord_values.iter().enumerate() {
                    residual += self.rref[base + col] * value;
                }
                max_row_residual = max_row_residual.max(residual.abs());
            }
        }
        if max_row_residual > 1.0e-8 {
            return Err(format!(
                "F_κ {source_label} solve violated KP15 rows by {max_row_residual:e}"
            ));
        }

        let mut grids = Vec::with_capacity(n);
        for k in 0..n {
            let mut grid = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
            for (i, row) in grid.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    let col = global_grid_col(n, k, i, j);
                    *cell = Point3::new(xs[col], ys[col], zs[col]);
                }
            }
            grids.push(grid);
        }
        let f_coord = |x: &[f64]| -> f64 {
            let mut acc = 0.0;
            for k in 0..n {
                let row_base = 49 * k;
                for r in 0..49 {
                    for c in 0..49 {
                        acc += x[row_base + r] * q_local[r][c] * x[row_base + c];
                    }
                }
            }
            acc
        };
        Ok((grids, [f_coord(xs), f_coord(ys), f_coord(zs)]))
    }
}

fn independent_affine_constraint_rows(
    cols: usize,
    rows: &[Vec<f64>],
    rhs: &[f64],
    tol: f64,
) -> Result<(Vec<Vec<f64>>, Vec<f64>), String> {
    assert_eq!(rows.len(), rhs.len());
    let mut basis = Vec::<Vec<f64>>::new();
    let mut targets = Vec::<f64>::new();
    for (row, &target) in rows.iter().zip(rhs) {
        assert_eq!(row.len(), cols);
        let mut candidate = row.clone();
        let mut candidate_target = target;
        for (basis_row, &basis_target) in basis.iter().zip(&targets) {
            let dot = candidate
                .iter()
                .zip(basis_row)
                .map(|(a, b)| a * b)
                .sum::<f64>();
            for (value, basis_value) in candidate.iter_mut().zip(basis_row) {
                *value -= dot * basis_value;
            }
            candidate_target -= dot * basis_target;
        }
        let norm = candidate.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm <= tol {
            if candidate_target.abs() > 1.0e-8 {
                return Err(format!(
                    "inconsistent KP15 two-sided tensor-source constraint residual {candidate_target:e}"
                ));
            }
            continue;
        }
        for value in &mut candidate {
            *value /= norm;
        }
        basis.push(candidate);
        targets.push(candidate_target / norm);
    }
    Ok((basis, targets))
}

/// Global column indices of [KP15 §3.6 step 4]'s "P" set — the
/// cap's free-parameter set in the paper's notation:
///
/// ```text
///   P := { p^0_{ij} : i + j ≤ 2 } ∪ { p^0_{22} }
///         ∪ { p^k_{23} : k ∈ 0..n } ∪ { p^k_{33} : k ∈ 0..n }
///         ∪ { p^0_{21} }                          [if n = 6]
/// ```
///
/// Size: `2n + 7` (or `2n + 8` for `n = 6`).
///
/// Consumed by
/// [`Kp15Bi6GeneratingTable::build_with_kp15_p_set_as_free`] to
/// arrange the RREF so P falls into the free column set of the
/// resulting generating table — making P the natural F_5
/// minimisation domain on the downstream call.
pub fn kp15_p_set_columns(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    // Sector 0's quadratic-near-star (i + j ≤ 2): 6 positions.
    for i in 0..=2 {
        for j in 0..=2 {
            if i + j <= 2 {
                out.push(global_grid_col(n, 0, i, j));
            }
        }
    }
    // p^0_{22}.
    out.push(global_grid_col(n, 0, 2, 2));
    // p^k_{23} per sector.
    for k in 0..n {
        out.push(global_grid_col(n, k, 2, 3));
    }
    // p^k_{33} per sector.
    for k in 0..n {
        out.push(global_grid_col(n, k, 3, 3));
    }
    // p^0_{21} at n = 6.
    if n == 6 {
        out.push(global_grid_col(n, 0, 2, 1));
    }
    out
}

/// Fixed tensor-border seed columns supplied by `CirculantSectorRhs`.
pub fn kp15_seed_columns(n: usize) -> Vec<usize> {
    let mut out = Vec::with_capacity(12 * n);
    for k in 0..n {
        for (i, j) in [
            (0, 0),
            (1, 0),
            (2, 0),
            (4, 0),
            (5, 0),
            (6, 0),
            (0, 1),
            (1, 1),
            (0, 2),
            (4, 2),
            (5, 2),
            (6, 2),
        ] {
            out.push(global_grid_col(n, k, i, j));
        }
    }
    out
}

/// Every column in the two outer tensor-product strips of the sector
/// chart.
///
/// This helper is retained for older diagnostics that compare the
/// closure effect of prescribing only a subset of the tensor-border
/// against prescribing the full depth-2 strip required by KP15 §3.6
/// step 3.
pub(crate) fn kp15_two_sided_tensor_l_band_columns(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for k in 0..n {
        for i in 0..=6 {
            for j in 0..=6 {
                if i >= 4 || j >= 4 {
                    out.push(global_grid_col(n, k, i, j));
                }
            }
        }
    }
    out
}

/// KP15 §3.6 step 3 full two-sided tensor-border source columns in
/// the Fig. 2b sector chart.
///
/// Appendix Fig. 12 defines the depth-2 strip `q := b o beta`; after
/// degree elevation to the bi-6 cap patch, step 3 fixes the full
/// tensor-border L-band of each sector. With the paper free set `P`
/// pinned separately, this full strip closes the Step-4 row system.
pub(crate) fn kp15_two_sided_tensor_border_columns(n: usize) -> Vec<usize> {
    kp15_two_sided_tensor_l_band_columns(n)
}

/// KP15 §3.2 S1 central quadratic expansion columns across all
/// sectors. These are derived from the paper's sector-0 source data
/// and are not independent fixed tensor-border pins.
pub fn kp15_central_quadratic_columns(n: usize) -> Vec<usize> {
    let mut out = Vec::with_capacity(6 * n);
    for k in 0..n {
        for i in 0..=2 {
            for j in 0..=2 {
                if i + j <= 2 {
                    out.push(global_grid_col(n, k, i, j));
                }
            }
        }
    }
    out
}

/// Bi-3 → bi-6 degree elevation for a 4-CP curve. Three iterated
/// degree raises, each preserving the underlying curve.
fn elevate_bi3_to_bi6(b: &[Point3; 4]) -> [Point3; 7] {
    let mut cur: Vec<Point3> = b.to_vec();
    for _ in 0..3 {
        let d = cur.len() - 1;
        let dp1 = (d + 1) as f64;
        let mut next = Vec::with_capacity(d + 2);
        next.push(cur[0]);
        for i in 1..=d {
            let frac = (i as f64) / dp1;
            let prev = cur[i - 1];
            let this = cur[i];
            next.push(Point3::new(
                frac * prev.x + (1.0 - frac) * this.x,
                frac * prev.y + (1.0 - frac) * this.y,
                frac * prev.z + (1.0 - frac) * this.z,
            ));
        }
        next.push(cur[d]);
        cur = next;
    }
    core::array::from_fn(|i| cur[i])
}

/// Weighted-sum helper: combine `(coef, point)` pairs into a single
/// `Point3` (coordinate-wise linear combination).
fn weighted_point_sum(pairs: &[(f64, Point3)]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for &(w, p) in pairs {
        out.x += w * p.x;
        out.y += w * p.y;
        out.z += w * p.z;
    }
    out
}

/// [KP15 §3.3 eq (17)] β reparameterization of the bi-3 tensor
/// border `b`. Returns the bi-6 reparameterized tensor border
/// `q := b ∘ β`: `q[j][i]` for `j ∈ 0..=2, i ∈ 0..=6`. By the
/// paper's construction:
/// - `q[0][i]` equals the bi-3 → bi-6 elevation of `b[0]`.
/// - `q[1][i]` preserves the first cross-derivative of the bi-3
///   tensor-border by cross-degree elevation:
///   `q[1] = (E(b[0]) + E(b[1])) / 2`.
/// - `q[2][i]` for `i ∈ 0..=1` use the ordinary cross-degree
///   elevation blend `(b[0], b[1], b[2]) -> (0.2, 0.6, 0.2)`.
/// - `q[2][i]` for `i ∈ 2..=6` are given by eq (17)'s γ-formula
///   weighted sums of the bi-3 BB-coefs.
///
/// Input `b` is indexed `b[j][i]` for `j ∈ 0..=2` (depth), `i ∈
/// 0..=3` (along-edge). Output `q[j][i]` follows the bi-6 grid
/// indexing the rest of the module uses.
pub fn apply_beta_reparameterization(
    b: &[[Point3; 4]; 3],
    params: &Bi6Reparams,
) -> [[Point3; 7]; 3] {
    let c = params.c;
    let tau = params.tau;
    let d2 = params.d2;

    // Step 1: cross-degree elevate the tensor-border BB coefficients.
    // Fig. 12 lists b_ij and q_ij as BB coefficients. Preserving the
    // first cross-derivative therefore makes q_1 the midpoint blend of
    // the elevated b_0 and b_1 rows, not the direct elevation of b_1.
    let mut q = [[Point3::new(0.0, 0.0, 0.0); 7]; 3];
    let raised0 = elevate_bi3_to_bi6(&b[0]);
    let raised1 = elevate_bi3_to_bi6(&b[1]);
    let raised2 = elevate_bi3_to_bi6(&b[2]);
    for i in 0..=6 {
        q[0][i] = raised0[i];
        q[1][i] = weighted_point_sum(&[(0.5, raised0[i]), (0.5, raised1[i])]);
        q[2][i] = weighted_point_sum(&[(0.2, raised0[i]), (0.6, raised1[i]), (0.2, raised2[i])]);
    }

    // Step 2: overwrite q[2][i] for i ∈ 2..=6 per eq (17).

    // γ²_{ij} (symmetric 3×3), independent of τ; involves d₂.
    let g2_00 = 1.0 - d2 / 25.0;
    let g2_10 = (6.0 + d2) / 50.0; // = γ²_{01}
    let g2_20 = 1.0 / 25.0; // = γ²_{02}
    let g2_11 = 9.0 / 25.0;
    let g2_21 = 3.0 / 25.0; // = γ²_{12}
    let g2_22 = 1.0 / 25.0;
    let g2: [[f64; 3]; 3] = [
        [g2_00, g2_10, g2_20],
        [g2_10, g2_11, g2_21],
        [g2_20, g2_21, g2_22],
    ];
    let mut pairs: Vec<(f64, Point3)> = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            pairs.push((g2[i][j], b[j][i]));
        }
    }
    q[2][2] = weighted_point_sum(&pairs);

    // γ³_{i, j} for i ∈ 0..=3, j ∈ 0..=2.
    let g3_0: [f64; 4] = [
        ((3.0 + 5.0 * c - 3.0 * d2) * tau - 6.0 * c) / (300.0 * tau),
        ((9.0 - 2.0 * c - 7.0 * d2) * tau + 2.0 * c) / (100.0 * tau),
        (9.0 + 3.0 * d2) / 100.0,
        1.0 / 100.0,
    ];
    let g3_1: [f64; 4] = [
        (18.0 + 2.0 * c + 3.0 * d2) / 600.0,
        (54.0 + 9.0 * d2) / 200.0,
        27.0 / 100.0,
        3.0 / 100.0,
    ];
    let g3_2: [f64; 4] = [1.0 / 100.0, 9.0 / 100.0, 9.0 / 100.0, 1.0 / 100.0];
    pairs.clear();
    for i in 0..4 {
        pairs.push((g3_0[i], b[0][i]));
        pairs.push((g3_1[i], b[1][i]));
        pairs.push((g3_2[i], b[2][i]));
    }
    q[2][3] = weighted_point_sum(&pairs);

    // γ⁴_{i, j} for i ∈ 1..=3, j ∈ 0..=2.
    let g4_0: [f64; 3] = [
        ((6.0 + 6.0 * c - 5.0 * d2) * tau - 8.0 * c) / (150.0 * tau),
        ((9.0 - 4.0 * c - 5.0 * d2) * tau + 4.0 * c) / (75.0 * tau),
        (2.0 + d2) / 50.0,
    ];
    let g4_1: [f64; 3] = [
        (18.0 + 2.0 * c + 3.0 * d2) / 150.0,
        (18.0 + 3.0 * d2) / 50.0,
        3.0 / 25.0,
    ];
    let g4_2: [f64; 3] = [1.0 / 25.0, 3.0 / 25.0, 1.0 / 25.0];
    pairs.clear();
    for i in 0..3 {
        let ii = i + 1;
        pairs.push((g4_0[i], b[0][ii]));
        pairs.push((g4_1[i], b[1][ii]));
        pairs.push((g4_2[i], b[2][ii]));
    }
    q[2][4] = weighted_point_sum(&pairs);

    // γ⁵_{i, j} for i ∈ 2..=3, j ∈ 0..=2.
    let g5_20 = ((3.0 + c - 2.0 * d2) * tau - 2.0 * c) / (30.0 * tau);
    let g5_30 = ((3.0 - 2.0 * c - d2) * tau + 2.0 * c) / (30.0 * tau);
    let g5_21 = (18.0 + 2.0 * c + 3.0 * d2) / 60.0;
    let g5_31 = (6.0 + d2) / 20.0;
    let g5_22 = 1.0 / 10.0;
    let g5_32 = 1.0 / 10.0;
    q[2][5] = weighted_point_sum(&[
        (g5_20, b[0][2]),
        (g5_30, b[0][3]),
        (g5_21, b[1][2]),
        (g5_31, b[1][3]),
        (g5_22, b[2][2]),
        (g5_32, b[2][3]),
    ]);

    // γ⁶_{3, j} for j ∈ 0..=2.
    let g6_30 = (6.0 - 2.0 * c - 3.0 * d2) / 30.0;
    let g6_31 = (18.0 + 2.0 * c + 3.0 * d2) / 30.0;
    let g6_32 = 1.0 / 5.0;
    q[2][6] = weighted_point_sum(&[(g6_30, b[0][3]), (g6_31, b[1][3]), (g6_32, b[2][3])]);

    q
}

/// Returns the `Point3` value at global position `(k, i, j)` if
/// that position is "seeded" by the cap input, else `None`.
/// Seed positions: `(0, 0)` → star; the 12 `CirculantSectorRhs`
/// fixed-position slots → the corresponding field.
fn seeded_value(
    sectors: &[super::kp15_bi6_circulant::CirculantSectorRhs],
    star: Point3,
    k: usize,
    i: usize,
    j: usize,
) -> Option<Point3> {
    let s = &sectors[k];
    Some(match (i, j) {
        (0, 0) => star,
        (1, 0) => s.p_10,
        (2, 0) => s.p_20,
        (4, 0) => s.p_40,
        (5, 0) => s.p_50,
        (6, 0) => s.p_60,
        (0, 1) => s.p_01,
        (1, 1) => s.p_11,
        (0, 2) => s.p_02,
        (4, 2) => s.p_42,
        (5, 2) => s.p_52,
        (6, 2) => s.p_62,
        _ => return None,
    })
}

#[cfg(test)]
#[path = "kp15_bi6_local_constraint_tests.rs"]
mod tests;
