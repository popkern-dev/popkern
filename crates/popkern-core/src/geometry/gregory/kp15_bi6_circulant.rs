//! Phase I.8 — K-P 2015 bi-6 circulant system solvers.
//!
//! Paper: [KP15 §3.4] (transcribed in `docs/journal/polar-kp-math.md`
//! §2.2). After rules (S1) pins `ṕ_{00..02}` and (S2) gives the
//! *relation* between adjacent sectors for `ṕ_{21}`, `ṕ_{12}`,
//! `ṕ_{22}`, going around the ring closes those relations into
//! linear systems:
//!
//! - A **`2n × 2n` circulant** from eq (7) determining `p^k_{21}`
//!   and `p^k_{12}` for all k simultaneously.
//! - An **`n × n` circulant** from eq (8) determining `p^k_{22}`.
//!
//! # Free parameters
//!
//! Per [KP15 §3.6] step 4, the free-parameter set `P` includes
//! `p⁰_{00..02}`, `p⁰_{22}`, and `p^k_{23}`, `p^k_{33}`.
//! Additionally, `p⁰_{21}` becomes a free parameter at `n = 6`
//! (circulant rank-deficient per [KP15 §3.4]) and at `n = 3`
//! (paper §3.5; for `n = 3`, `p¹_{22}` and `p²_{22}` are also
//! free — our p22 solver skips the whole ring for `n = 3`).
//!
//! # Algorithm
//!
//! 1. Build the `2n × 2n` matrix from eq (7) with the µ̇ / µ̈
//!    coefficients in the LHS and per-sector row-0/row-1 known
//!    data contracted into a 3-coord RHS vector per axis.
//! 2. For `n ∈ {3, 6}`, replace row 0 with the identity row
//!    pinning `p⁰_{21}` to its free-parameter value.
//! 3. Solve via [`solve_dense`] three times (x, y, z).
//! 4. Extract `p^k_{21}`, `p^k_{12}` from the flat solution.
//!
//! The `n × n` p22 solve follows the same pattern with eq (8).
//!
//! # Complexity
//!
//! Naive `O(n³)` dense solve for each axis. Our target `n ≤ 12`
//! means 24×24 at worst, negligible. A DFT-based `O(n log n)`
//! solver is possible but not needed.

use super::direct_linalg::solve_dense;
use super::kp15_bi6_reparam::Bi6Reparams;
use crate::geometry::Point3;

/// Known-row data per sector, in the frame of that sector. The
/// circulant solver consumes this for every `k ∈ 0..n`.
#[derive(Debug, Clone, Copy)]
pub struct CirculantSectorRhs {
    /// `p^k_{00}` — star point (identical across sectors, but
    /// carried per-sector for API symmetry).
    pub p_00: Point3,
    /// `p^k_{10}` — shared edge.
    pub p_10: Point3,
    /// `p^k_{20}` — shared edge.
    pub p_20: Point3,
    /// `p^k_{40}` — shared edge.
    pub p_40: Point3,
    /// `p^k_{50}` — shared edge.
    pub p_50: Point3,
    /// `p^k_{60}` — shared edge.
    pub p_60: Point3,
    /// `p^k_{01}` (from (S1), or free at k = 0).
    pub p_01: Point3,
    /// `p^k_{11}` (from (S1), or free at k = 0).
    pub p_11: Point3,
    /// `p^k_{02}` (from (S1), or free at k = 0).
    pub p_02: Point3,
    /// `p^k_{42}` — row-2 tensor-border at i = 4.
    pub p_42: Point3,
    /// `p^k_{52}` — row-2 tensor-border at i = 5.
    pub p_52: Point3,
    /// `p^k_{62}` — row-2 tensor-border at i = 6.
    pub p_62: Point3,
}

/// Result of [`solve_bi6_circulant_s2_eq7`] — per-sector
/// assignments for `p^k_{21}` and `p^k_{12}`.
#[derive(Debug, Clone)]
pub struct S2Eq7Result {
    /// `p^k_{21}` for k = 0..n.
    pub p_21: Vec<Point3>,
    /// `p^k_{12}` for k = 0..n.
    pub p_12: Vec<Point3>,
}

/// Solve the `2n × 2n` circulant system from [KP15 eq (7)] for
/// `p^k_{21}` and `p^k_{12}`, k = 0..n − 1.
///
/// `sectors[k]` must carry the row-0 / row-1 BB-CPs for sector
/// `k` (the star point, shared-edge CPs, and the (S1)-derived
/// `p^k_{01, 11, 02}`).
///
/// For `n ∈ {3, 6}` the system is rank-deficient; `p0_21_pinned`
/// must be `Some(value)` providing the free-parameter choice
/// for `p⁰_{21}`. For other `n` the parameter is ignored.
///
/// Returns an error if the matrix is singular for `n` where it
/// shouldn't be.
pub fn solve_bi6_circulant_s2_eq7(
    params: &Bi6Reparams,
    mu_dot: [f64; 5],
    mu_ddot: [f64; 8],
    sectors: &[CirculantSectorRhs],
    p0_21_pinned: Option<Point3>,
) -> Result<S2Eq7Result, String> {
    let n = sectors.len();
    if n < 3 {
        return Err(format!("bi-6 circulant needs n ≥ 3, got {n}"));
    }
    if params.n != n {
        return Err(format!("params.n={} but sectors.len()={n}", params.n));
    }
    let needs_pin = n == 3 || n == 6;
    if needs_pin && p0_21_pinned.is_none() {
        return Err(format!(
            "n = {n}: p0_21_pinned must be Some(free parameter value)"
        ));
    }

    // Flat index layout: [x^0_{21}, ..., x^{n-1}_{21}, x^0_{12}, ..., x^{n-1}_{12}]
    //
    // Paper notation vs global indices [KP15 §3.2, p184–189 note]:
    // Adjacent sectors share a u↔v swap, so the (p̀, ṕ) local
    // frame maps to global CPs as:
    //   p̀_{21} = p^{k-1}_{12}
    //   p̀_{12} = p^{k-1}_{21}
    //   ṕ_{21} = p^k_{21}
    //   ṕ_{12} = p^k_{12}
    // (cf. paper line 186–189: "pk−1_{21} = p̀_{12}, pk−1_{12} =
    // p̀_{21}"). Substituting into eq (7) gives the coupled system:
    //   p^k_{21} + p^{k-1}_{12}                       = rhs_21_k
    //   p^k_{12} − p^{k-1}_{21} − µ̈_5·p^{k-1}_{12}   = rhs_12_k
    let nvars = 2 * n;
    let prev_of = |k: usize| -> usize { (k + n - 1) % n };

    let sum_mu_dot: f64 = mu_dot.iter().sum();
    let sum_mu_ddot: f64 = mu_ddot.iter().sum();

    // Build the matrix once (same for all three coords).
    let mut mat: Vec<Vec<f64>> = vec![vec![0.0; nvars]; nvars];
    for k in 0..n {
        let km1 = prev_of(k);

        // (21) eq k: p^k_{21} + p^{k-1}_{12} = rhs_21_k
        let row_21 = k;
        mat[row_21][k] = 1.0;
        mat[row_21][n + km1] += 1.0;

        // (12) eq k: p^k_{12} − p^{k-1}_{21} − µ̈_5·p^{k-1}_{12} = rhs_12_k
        let row_12 = n + k;
        mat[row_12][n + k] = 1.0;
        mat[row_12][km1] -= 1.0;
        mat[row_12][n + km1] -= mu_ddot[4]; // µ̈_5
    }

    // Handle pinned cases (n ∈ {3, 6}): replace row 0 (the
    // `x^0_{21}` equation) with the identity row pinning it to
    // `p0_21_pinned`.
    if needs_pin {
        mat[0].iter_mut().take(nvars).for_each(|cell| *cell = 0.0);
        mat[0][0] = 1.0;
    }

    // Build per-axis RHS and solve.
    let build_rhs = |coord: fn(Point3) -> f64| -> Vec<f64> {
        let mut b = vec![0.0; nvars];
        for k in 0..n {
            let km1 = prev_of(k);
            let s = &sectors[km1]; // p̀ = p^{k-1}
            // rhs_21_k
            b[k] = (2.0 - sum_mu_dot) * coord(s.p_00)
                + mu_dot[0] * coord(s.p_10)
                + mu_dot[1] * coord(s.p_20)
                + mu_dot[2] * coord(s.p_40)
                + mu_dot[3] * coord(s.p_50)
                + mu_dot[4] * coord(s.p_60);
            // rhs_12_k (µ̈_5·p̀_{21} moves to LHS; not here).
            b[n + k] = -sum_mu_ddot * coord(s.p_00)
                + mu_ddot[0] * coord(s.p_01)
                + mu_ddot[1] * coord(s.p_10)
                + mu_ddot[2] * coord(s.p_11)
                + mu_ddot[3] * coord(s.p_20)
                + mu_ddot[5] * coord(s.p_40)
                + mu_ddot[6] * coord(s.p_50)
                + mu_ddot[7] * coord(s.p_60);
        }
        if needs_pin && let Some(pinned) = p0_21_pinned {
            b[0] = coord(pinned);
        }
        b
    };

    let sol_x = solve_dense(mat.clone(), build_rhs(|p| p.x))
        .ok_or_else(|| format!("bi-6 (S2) eq (7) system singular for n = {n}"))?;
    let sol_y = solve_dense(mat.clone(), build_rhs(|p| p.y))
        .ok_or_else(|| format!("bi-6 (S2) eq (7) system singular for n = {n}"))?;
    let sol_z = solve_dense(mat, build_rhs(|p| p.z))
        .ok_or_else(|| format!("bi-6 (S2) eq (7) system singular for n = {n}"))?;

    let p_21 = (0..n)
        .map(|k| Point3::new(sol_x[k], sol_y[k], sol_z[k]))
        .collect();
    let p_12 = (0..n)
        .map(|k| Point3::new(sol_x[n + k], sol_y[n + k], sol_z[n + k]))
        .collect();
    Ok(S2Eq7Result { p_21, p_12 })
}

/// Solve the `n × n` circulant system from [KP15 eq (8)] for
/// `p^k_{22}`, k = 0..n − 1.
///
/// `p^k_{22}`, k = 0..n − 1. `p0_22_pinned` pins `p⁰_{22}` —
/// per [KP15 §3.6 step 4] this is one of the free parameters in
/// the set `P`.
///
/// `p_12` from a prior [`solve_bi6_circulant_s2_eq7`] call is
/// used for the `⃛µ_5 · p̀_{21}` term. Per the [KP15 §3.2
/// p184–189] u-v swap convention, `p̀_{21}` in the pair-frame
/// maps to global `q^{k-1}_{12}` (= caller's `p_12[k-1]`), not
/// `q^{k-1}_{21}`. See module-level docs for the full swap
/// mapping.
///
/// For `n = 3` the entire `p^k_{22}` ring is free (per §3.5);
/// callers should not invoke this solver when `n = 3`.
pub fn solve_bi6_circulant_s2_eq8(
    params: &Bi6Reparams,
    mu_tdot: [f64; 11],
    sectors: &[CirculantSectorRhs],
    p_12: &[Point3],
    p0_22_pinned: Point3,
) -> Result<Vec<Point3>, String> {
    let n = sectors.len();
    if n < 3 {
        return Err(format!("bi-6 circulant p22 needs n ≥ 3, got {n}"));
    }
    if n == 3 {
        return Err(
            "bi-6 circulant p22 not solvable at n = 3; all p^k_{22} are free (§3.5)".to_string(),
        );
    }
    if params.n != n {
        return Err(format!("params.n={} but sectors.len()={n}", params.n));
    }
    if p_12.len() != n {
        return Err(format!("p_12 length {} ≠ n = {n}", p_12.len()));
    }

    let nvars = n;
    let prev_of = |k: usize| -> usize { (k + n - 1) % n };
    let sum_mu_1_to_8: f64 = mu_tdot.iter().take(8).sum();

    // Matrix: eq k is x^k_{22} - x^{prev}_{22} = rhs_22_k.
    // Pin row 0 as x^0_{22} = p0_22_pinned.
    let mut mat: Vec<Vec<f64>> = vec![vec![0.0; nvars]; nvars];
    mat[0][0] = 1.0;
    // `mat[k][k]` and `mat[k][km1]` index a different row via
    // km1 = (k-1)%n, so the `needless_range_loop` rewrite to
    // iter_mut().enumerate() doesn't apply cleanly.
    #[allow(clippy::needless_range_loop)]
    for k in 1..n {
        let km1 = prev_of(k);
        mat[k][k] = 1.0;
        mat[k][km1] -= 1.0;
    }

    let build_rhs = |coord: fn(Point3) -> f64| -> Vec<f64> {
        let mut b = vec![0.0; nvars];
        b[0] = coord(p0_22_pinned);
        for k in 1..n {
            let km1 = prev_of(k);
            let s = &sectors[km1];
            let s_cur = &sectors[k]; // ṕ = p^k for the diff terms.
            // rhs_22_k per eq (8):
            //   (-Σ⃛µ_{1..8})·p̀_{00}
            //   + ⃛µ_1·p̀_{01} + ⃛µ_2·p̀_{10} + ⃛µ_3·p̀_{11}
            //   + ⃛µ_4·p̀_{20} + ⃛µ_5·p̀_{21} + ⃛µ_6·p̀_{40}
            //   + ⃛µ_7·p̀_{50} + ⃛µ_8·p̀_{60}
            //   + ⃛µ_9·(p̀_{42} − ṕ_{42}) + ⃛µ_{10}·(p̀_{52} − ṕ_{52})
            //   + ⃛µ_{11}·(p̀_{62} − ṕ_{62})
            b[k] = -sum_mu_1_to_8 * coord(s.p_00)
                + mu_tdot[0] * coord(s.p_01)
                + mu_tdot[1] * coord(s.p_10)
                + mu_tdot[2] * coord(s.p_11)
                + mu_tdot[3] * coord(s.p_20)
                + mu_tdot[4] * coord(p_12[km1])
                + mu_tdot[5] * coord(s.p_40)
                + mu_tdot[6] * coord(s.p_50)
                + mu_tdot[7] * coord(s.p_60)
                + mu_tdot[8] * (coord(s.p_42) - coord(s_cur.p_42))
                + mu_tdot[9] * (coord(s.p_52) - coord(s_cur.p_52))
                + mu_tdot[10] * (coord(s.p_62) - coord(s_cur.p_62));
        }
        b
    };

    let sol_x = solve_dense(mat.clone(), build_rhs(|p| p.x))
        .ok_or_else(|| format!("bi-6 (S2) eq (8) system singular for n = {n}"))?;
    let sol_y = solve_dense(mat.clone(), build_rhs(|p| p.y))
        .ok_or_else(|| format!("bi-6 (S2) eq (8) system singular for n = {n}"))?;
    let sol_z = solve_dense(mat, build_rhs(|p| p.z))
        .ok_or_else(|| format!("bi-6 (S2) eq (8) system singular for n = {n}"))?;

    Ok((0..n)
        .map(|k| Point3::new(sol_x[k], sol_y[k], sol_z[k]))
        .collect())
}

#[cfg(test)]
#[path = "kp15_bi6_circulant_tests.rs"]
mod tests;
