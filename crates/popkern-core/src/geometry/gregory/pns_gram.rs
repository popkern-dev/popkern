//! PnS cap — F_κ Gram matrix + constrained minimiser for bi-4
//! patches.
//!
//! Builds on [`pns_f_kappa`](super::pns_f_kappa) to express the
//! quadratic form
//!
//! `F_κ(p) = ∫∫ Σ_{i+j=κ} (κ!/(i!·j!)) · (∂^i_s ∂^j_t p)² ds dt`
//!
//! explicitly as `F_κ(p) = p^T · G^{(κ)} · p` (per coordinate),
//! where `G^{(κ)}` is a `25 × 25` symmetric PSD matrix independent
//! of `p`. Flattening index: `(a, b) → 5·a + b`, where `a ∈ 0..=4`
//! is the `s`-BB index and `b ∈ 0..=4` is the `t`-BB index.
//!
//! # Separability
//!
//! `G^{(κ)}` decomposes as
//!
//! `G^{(κ)}_{(a,b),(c,d)} = Σ_{i=0}^{κ} (κ!/(i!·(κ−i)!)) · H_i[a,c] · H_{κ−i}[b,d]`
//!
//! where `H_i[a, c] = ∫_0^1 (d^i B_a^4/du^i)(u) · (d^i B_c^4/du^i)(u) du`
//! is the 1D Gram matrix of the `i`-th-order Bernstein derivatives.
//! Each `H_i` is `5 × 5`, symmetric PSD.
//!
//! # Constrained minimiser
//!
//! The PnS cap's F_κ-minimisation problem (eq 6 of `14bi234.pdf`)
//! partitions the `25` BB-coefs into:
//! - **Fixed** coefs (from G¹ constraints, b∘β, central expansion,
//!   eq 4 p_22)
//! - **Free** coefs (`n + 5` of them per the paper)
//!
//! Minimising `F_κ` over free coefs with fixed coefs pinned is the
//! standard quadratic KKT problem: setting `∂F_κ/∂p_free = 0`
//! yields the linear system
//!
//! `G_vv · p_free = −G_vf · p_fixed`
//!
//! where `G_vv` is the free-block submatrix of G, `G_vf` maps
//! fixed coefs to their contribution on the free gradient. Solved
//! via [`super::direct_linalg::solve_dense`].
//!
//! Each coordinate (x, y, z) gets its own RHS but shares the same
//! system matrix `G_vv` — we factorise once (conceptually; this
//! implementation re-solves for simplicity, still fast at ≤25×25).

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::direct_linalg::solve_dense;
use crate::geometry::Point3;

/// Flatten `(a, b) ∈ {0..=4}²` to a single index `0..=24` in
/// row-major order: `idx = 5·a + b`.
#[inline]
pub fn bi4_flat_index(a: usize, b: usize) -> usize {
    debug_assert!(a < 5 && b < 5);
    5 * a + b
}

/// Unflatten `idx ∈ 0..=24` back to `(a, b)`.
#[inline]
pub fn bi4_unflat(idx: usize) -> (usize, usize) {
    debug_assert!(idx < 25);
    (idx / 5, idx % 5)
}

/// Full `25 × 25` Gram matrix of the F_κ bilinear form on a bi-4
/// patch, for `κ ∈ 1..=4`.
///
/// The matrix is symmetric positive semi-definite. The entry
/// `G[idx(a, b), idx(c, d)]` equals the coefficient of
/// `p_{a,b}·p_{c,d}` in `F_κ(p)` (per coordinate).
pub fn gram_bi4_f_kappa(kappa: usize) -> [[f64; 25]; 25] {
    assert!(
        (1..=4).contains(&kappa),
        "gram_bi4_f_kappa: kappa must be in 1..=4, got {kappa}"
    );
    let mut g = [[0.0_f64; 25]; 25];
    // Precompute the 1D Gram matrices H_i for i = 0..=κ.
    let h: Vec<[[f64; 5]; 5]> = (0..=kappa).map(h_matrix_bi4_deriv).collect();
    for i in 0..=kappa {
        let j = kappa - i;
        let coef = factorial(kappa) as f64 / (factorial(i) as f64 * factorial(j) as f64);
        let hi = &h[i];
        let hj = &h[j];
        for a in 0..5 {
            for b in 0..5 {
                for c in 0..5 {
                    for d in 0..5 {
                        let idx_ab = bi4_flat_index(a, b);
                        let idx_cd = bi4_flat_index(c, d);
                        g[idx_ab][idx_cd] += coef * hi[a][c] * hj[b][d];
                    }
                }
            }
        }
    }
    g
}

/// Minimise `F_κ` over a subset of free BB-coefs, keeping the rest
/// fixed at their current values.
///
/// `fixed_mask[idx]` is `true` when that BB-coef is pinned (input
/// value preserved), `false` when it's a free DoF to be optimised.
///
/// On success, `grid` is updated with the optimal free values.
///
/// Returns `Err` if the free-block submatrix `G_vv` is singular to
/// machine precision (indicates ill-posed problem — usually means
/// the constraint set leaves fewer DoFs than expected).
pub fn minimize_f_kappa_bi4(
    grid: &mut [[Point3; 5]; 5],
    kappa: usize,
    fixed_mask: &[bool; 25],
) -> Result<(), String> {
    let g = gram_bi4_f_kappa(kappa);

    // Collect free indices + extract fixed values.
    let mut free_indices: Vec<usize> = Vec::new();
    let mut fixed_indices: Vec<usize> = Vec::new();
    for idx in 0..25 {
        if fixed_mask[idx] {
            fixed_indices.push(idx);
        } else {
            free_indices.push(idx);
        }
    }
    if free_indices.is_empty() {
        return Ok(());
    }
    let n_free = free_indices.len();

    // Build G_vv (free × free submatrix) once — same for all coords.
    let mut g_vv: Vec<Vec<f64>> = vec![vec![0.0; n_free]; n_free];
    for (i, &row_idx) in free_indices.iter().enumerate() {
        for (j, &col_idx) in free_indices.iter().enumerate() {
            g_vv[i][j] = g[row_idx][col_idx];
        }
    }

    // Solve per coordinate.
    for coord in 0..3 {
        // RHS: -G_vf · p_fixed  (n_free vector)
        let mut rhs: Vec<f64> = vec![0.0; n_free];
        for (i, &row_idx) in free_indices.iter().enumerate() {
            for &fix_idx in &fixed_indices {
                let (a, b) = bi4_unflat(fix_idx);
                let p_fixed = match coord {
                    0 => grid[a][b].x,
                    1 => grid[a][b].y,
                    _ => grid[a][b].z,
                };
                rhs[i] -= g[row_idx][fix_idx] * p_fixed;
            }
        }

        let solved = solve_dense(g_vv.clone(), rhs).ok_or_else(|| {
            format!(
                "minimize_f_kappa_bi4: G_vv singular (coord={coord}, kappa={kappa}, \
                 n_free={n_free})"
            )
        })?;
        // Write free values back.
        for (i, &free_idx) in free_indices.iter().enumerate() {
            let (a, b) = bi4_unflat(free_idx);
            match coord {
                0 => grid[a][b].x = solved[i],
                1 => grid[a][b].y = solved[i],
                _ => grid[a][b].z = solved[i],
            }
        }
    }
    Ok(())
}

// ── 1D Bernstein Gram matrices ───────────────────────────────────

/// `H_i[a, c] := ∫_0^1 (d^i B_a^4/du^i)(u) · (d^i B_c^4/du^i)(u) du`
/// for `a, c ∈ 0..=4`, `i ∈ 0..=4`.
fn h_matrix_bi4_deriv(i: usize) -> [[f64; 5]; 5] {
    // 5-point Gauss-Legendre is exact for polynomials up to degree 9.
    // (d^i B_a^4)² has degree 2·(4 − i), max 8 at i=0. Always ≤ 9. ✓
    let (nodes, weights) = gauss_legendre_5_unit_interval();
    let mut h = [[0.0_f64; 5]; 5];
    for (k, &u) in nodes.iter().enumerate() {
        let w = weights[k];
        let b = bernstein_bi4_deriv_order(i, u);
        for a in 0..5 {
            for c in 0..5 {
                h[a][c] += w * b[a] * b[c];
            }
        }
    }
    h
}

/// Bernstein bi-4 derivatives (all 5 basis values, up to 4th order).
fn bernstein_bi4_deriv_order(k: usize, u: f64) -> [f64; 5] {
    match k {
        0 => {
            let omu = 1.0 - u;
            let c: [f64; 5] = [1.0, 4.0, 6.0, 4.0, 1.0];
            let mut out = [0.0_f64; 5];
            let mut tp = 1.0_f64;
            for i in 0..5 {
                let om_pow = omu.powi((4 - i) as i32);
                out[i] = c[i] * tp * om_pow;
                tp *= u;
            }
            out
        }
        1 => {
            let b3 = bi3_all(u);
            [
                -4.0 * b3[0],
                4.0 * (b3[0] - b3[1]),
                4.0 * (b3[1] - b3[2]),
                4.0 * (b3[2] - b3[3]),
                4.0 * b3[3],
            ]
        }
        2 => {
            let b2 = bi2_all(u);
            let bm = |k: i32| -> f64 {
                if (0..=2).contains(&k) {
                    b2[k as usize]
                } else {
                    0.0
                }
            };
            let mut out = [0.0_f64; 5];
            for i in 0..5 {
                let ii = i as i32;
                out[i] = 12.0 * (bm(ii - 2) - 2.0 * bm(ii - 1) + bm(ii));
            }
            out
        }
        3 => {
            let b1 = [1.0 - u, u];
            let bm = |k: i32| -> f64 {
                if (0..=1).contains(&k) {
                    b1[k as usize]
                } else {
                    0.0
                }
            };
            let mut out = [0.0_f64; 5];
            for i in 0..5 {
                let ii = i as i32;
                out[i] = 24.0 * (bm(ii - 3) - 3.0 * bm(ii - 2) + 3.0 * bm(ii - 1) - bm(ii));
            }
            out
        }
        4 => [24.0, -96.0, 144.0, -96.0, 24.0],
        _ => [0.0; 5],
    }
}

fn bi3_all(u: f64) -> [f64; 4] {
    let omu = 1.0 - u;
    [
        omu * omu * omu,
        3.0 * u * omu * omu,
        3.0 * u * u * omu,
        u * u * u,
    ]
}

fn bi2_all(u: f64) -> [f64; 3] {
    let omu = 1.0 - u;
    [omu * omu, 2.0 * u * omu, u * u]
}

fn gauss_legendre_5_unit_interval() -> ([f64; 5], [f64; 5]) {
    let a1 = (5.0 - 2.0 * (10.0_f64 / 7.0).sqrt()).sqrt() / 3.0;
    let a2 = (5.0 + 2.0 * (10.0_f64 / 7.0).sqrt()).sqrt() / 3.0;
    let w0 = 128.0 / 225.0;
    let w1 = (322.0 + 13.0 * 70.0_f64.sqrt()) / 900.0;
    let w2 = (322.0 - 13.0 * 70.0_f64.sqrt()) / 900.0;
    let nodes = [
        (-a2 + 1.0) * 0.5,
        (-a1 + 1.0) * 0.5,
        0.5,
        (a1 + 1.0) * 0.5,
        (a2 + 1.0) * 0.5,
    ];
    let weights = [w2 * 0.5, w1 * 0.5, w0 * 0.5, w1 * 0.5, w2 * 0.5];
    (nodes, weights)
}

fn factorial(n: usize) -> usize {
    (1..=n).product::<usize>().max(1)
}

#[cfg(test)]
#[path = "pns_gram_tests.rs"]
mod tests;
