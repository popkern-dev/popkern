//! Phase H.5 — least-squares interior CP fit.
//!
//! Given a bidegree-9 patch with pinned boundary + row-1/row-2
//! frame CPs (64 CPs), fit the 36 interior CPs (rows 3..=6, cols
//! 3..=6) to minimise the L² distance to the bicubic guide:
//!
//! ```text
//!   min_{P[3..=6][3..=6]} ∫₀¹∫₀¹ ‖S(u, v) − G(u, v)‖² du dv
//! ```
//!
//! where `S(u, v) = Σ P[i][j] · B_i^9(u) · B_j^9(v)` (10×10
//! bidegree-9 tensor-product Bezier) and `G(u, v)` is the
//! bicubic guide (4×4 tensor-product Bezier).
//!
//! # Derivation
//!
//! Setting the gradient w.r.t. each unknown `P[i][j]` to zero
//! gives the normal equations: for `(i*, j*) ∈ U × U` where
//! `U = {3, 4, 5, 6}`:
//!
//! ```text
//!   Σ_{i, j ∈ 0..=9} P[i][j] · M[i][i*] · M[j][j*]
//!     = Σ_{k, l ∈ 0..=3} G[k][l] · H[k][i*] · H[l][j*]
//! ```
//!
//! where
//!
//! ```text
//!   M[a][b] = ∫ B_a^9(t) B_b^9(t) dt  = C(9,a) C(9,b) / (19 C(18, a+b))
//!   H[k][i] = ∫ B_k^3(t) B_i^9(t) dt  = C(3,k) C(9,i) / (13 C(12, k+i))
//! ```
//!
//! Splitting the LHS into `U × U` unknowns and pinned contributions:
//!
//! ```text
//!   M_UU · X · M_UU^T = RHS - K
//! ```
//!
//! where `M_UU` is the 4×4 sub-block `M[3..=6][3..=6]`, `X` is
//! the 4×4 unknown CPs (per coordinate), `RHS` is the guide's
//! projection, and `K` is the pinned-CP contribution. Solved via
//! two 4×4 Gaussian-elimination passes (symmetric positive-
//! definite Gram → stable).

#![allow(clippy::needless_range_loop)]

use super::direct_linalg::solve_dense;
use super::kp_guide_patch::KpGuidePatch;
use super::kp_guide_patch_binonic::KpGuidePatchBiNonic;
use super::kp_guide_patch_biquintic::KpGuidePatchBiQuintic;
use super::rectangular_binonic::RectangularGregoryBiNonic;

/// Indices of the unknown interior block along each axis.
const U_INDICES: [usize; 4] = [3, 4, 5, 6];

/// Fit the 36 interior CPs (rows 3..=6, cols 3..=6) of `patch`
/// to the `guide` via least-squares. Mutates `patch.grid` in
/// place — only the 36 interior entries are written; the 64
/// frame entries are untouched.
pub fn fit_interior_cps_to_guide(patch: &mut RectangularGregoryBiNonic, guide: &KpGuidePatch) {
    let m = bernstein_gram_9();
    let h = bernstein_cross_gram_3_9();
    // M restricted to U × U.
    let mut m_uu = [[0.0_f64; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            m_uu[a][b] = m[U_INDICES[a]][U_INDICES[b]];
        }
    }

    // Per-coordinate solve. Each coord is independent because
    // the Bernstein basis scales linearly with position.
    for coord in 0..3 {
        // Guide CPs for this coordinate (4×4 f64).
        let mut g_coord = [[0.0_f64; 4]; 4];
        for k in 0..4 {
            for l in 0..4 {
                let p = guide.grid[k][l];
                g_coord[k][l] = match coord {
                    0 => p.x,
                    1 => p.y,
                    _ => p.z,
                };
            }
        }
        // Pinned CPs for this coordinate (10×10 with holes filled
        // later). Only positions in U × U are unknown; everywhere
        // else the value is already set on `patch.grid`.
        let mut p_coord = [[0.0_f64; 10]; 10];
        for i in 0..10 {
            for j in 0..10 {
                let p = patch.grid[i][j];
                p_coord[i][j] = match coord {
                    0 => p.x,
                    1 => p.y,
                    _ => p.z,
                };
            }
        }

        // RHS[i*][j*] = Σ_{k, l} G[k][l] · H[k][i*] · H[l][j*],
        // restricted to (i*, j*) ∈ U × U.
        let mut rhs = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                let i_star = U_INDICES[a];
                let j_star = U_INDICES[b];
                let mut sum = 0.0_f64;
                for k in 0..4 {
                    for l in 0..4 {
                        sum += g_coord[k][l] * h[k][i_star] * h[l][j_star];
                    }
                }
                rhs[a][b] = sum;
            }
        }

        // K[i*][j*] = Σ_{(i, j) ∈ all \ U×U} P[i][j] · M[i][i*] · M[j][j*].
        // Equivalently: (Σ over all i, j) minus (Σ over i, j ∈ U × U).
        // Since the unknown U × U entries of p_coord are zero at
        // the time of this call (the caller hasn't written them),
        // subtracting nothing from the full sum gives exactly K.
        let mut k_mat = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                let i_star = U_INDICES[a];
                let j_star = U_INDICES[b];
                let mut sum = 0.0_f64;
                for i in 0..10 {
                    for j in 0..10 {
                        if U_INDICES.contains(&i) && U_INDICES.contains(&j) {
                            continue;
                        }
                        sum += p_coord[i][j] * m[i][i_star] * m[j][j_star];
                    }
                }
                k_mat[a][b] = sum;
            }
        }

        // Effective RHS for the U × U system.
        let mut eff_rhs = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                eff_rhs[a][b] = rhs[a][b] - k_mat[a][b];
            }
        }

        // Solve M_UU · X · M_UU = eff_rhs.
        // Step 1: Y = M_UU^{-1} · eff_rhs (solve 4 RHS columns).
        // The Bernstein Gram is a constant symmetric positive-
        // definite matrix — `solve_dense` is guaranteed to
        // succeed, so `unreachable!` guards the algebraic
        // contract without violating the no-unwrap structural
        // gate.
        let mut y = [[0.0_f64; 4]; 4];
        for col in 0..4 {
            let matrix: Vec<Vec<f64>> = m_uu.iter().map(|row| row.to_vec()).collect();
            let b_vec: Vec<f64> = (0..4).map(|r| eff_rhs[r][col]).collect();
            let Some(x) = solve_dense(matrix, b_vec) else {
                unreachable!("Bernstein Gram M_UU is constant and SPD")
            };
            for r in 0..4 {
                y[r][col] = x[r];
            }
        }
        // Step 2: X · M_UU = Y → X = Y · M_UU^{-1}. Since M_UU is
        // symmetric, this is equivalent to M_UU · X^T = Y^T.
        let mut x_final = [[0.0_f64; 4]; 4];
        for row in 0..4 {
            let matrix: Vec<Vec<f64>> = m_uu.iter().map(|r| r.to_vec()).collect();
            let b_vec: Vec<f64> = (0..4).map(|c| y[row][c]).collect();
            let Some(x) = solve_dense(matrix, b_vec) else {
                unreachable!("Bernstein Gram M_UU is constant and SPD")
            };
            x_final[row].copy_from_slice(&x[..4]);
        }

        // Write the solved interior CPs back to the patch.
        for a in 0..4 {
            for b in 0..4 {
                let i = U_INDICES[a];
                let j = U_INDICES[b];
                let val = x_final[a][b];
                match coord {
                    0 => patch.grid[i][j].x = val,
                    1 => patch.grid[i][j].y = val,
                    _ => patch.grid[i][j].z = val,
                }
            }
        }
    }
}

/// Fit the 36 interior CPs of `patch` to a **bi-quintic** guide
/// via least-squares. Same algorithm as [`fit_interior_cps_to_guide`]
/// but with a 6×6 source grid and a degree-5 vs degree-9
/// cross-Gram matrix. Mutates only the 36 interior entries of
/// `patch.grid`.
pub fn fit_interior_cps_to_biquintic_guide(
    patch: &mut RectangularGregoryBiNonic,
    guide: &KpGuidePatchBiQuintic,
) {
    let m = bernstein_gram_9();
    let h = bernstein_cross_gram_5_9();
    let mut m_uu = [[0.0_f64; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            m_uu[a][b] = m[U_INDICES[a]][U_INDICES[b]];
        }
    }

    for coord in 0..3 {
        // Guide CPs for this coordinate (6×6 f64).
        let mut g_coord = [[0.0_f64; 6]; 6];
        for k in 0..6 {
            for l in 0..6 {
                let p = guide.grid[k][l];
                g_coord[k][l] = match coord {
                    0 => p.x,
                    1 => p.y,
                    _ => p.z,
                };
            }
        }
        let mut p_coord = [[0.0_f64; 10]; 10];
        for i in 0..10 {
            for j in 0..10 {
                let p = patch.grid[i][j];
                p_coord[i][j] = match coord {
                    0 => p.x,
                    1 => p.y,
                    _ => p.z,
                };
            }
        }

        // RHS[i*][j*] = Σ_{k, l in 0..=5} G[k][l] · H[k][i*] · H[l][j*].
        let mut rhs = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                let i_star = U_INDICES[a];
                let j_star = U_INDICES[b];
                let mut sum = 0.0_f64;
                for k in 0..6 {
                    for l in 0..6 {
                        sum += g_coord[k][l] * h[k][i_star] * h[l][j_star];
                    }
                }
                rhs[a][b] = sum;
            }
        }

        let mut k_mat = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                let i_star = U_INDICES[a];
                let j_star = U_INDICES[b];
                let mut sum = 0.0_f64;
                for i in 0..10 {
                    for j in 0..10 {
                        if U_INDICES.contains(&i) && U_INDICES.contains(&j) {
                            continue;
                        }
                        sum += p_coord[i][j] * m[i][i_star] * m[j][j_star];
                    }
                }
                k_mat[a][b] = sum;
            }
        }

        let mut eff_rhs = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                eff_rhs[a][b] = rhs[a][b] - k_mat[a][b];
            }
        }

        let mut y = [[0.0_f64; 4]; 4];
        for col in 0..4 {
            let matrix: Vec<Vec<f64>> = m_uu.iter().map(|row| row.to_vec()).collect();
            let b_vec: Vec<f64> = (0..4).map(|r| eff_rhs[r][col]).collect();
            let Some(x) = solve_dense(matrix, b_vec) else {
                unreachable!("Bernstein Gram M_UU is constant and SPD")
            };
            for r in 0..4 {
                y[r][col] = x[r];
            }
        }
        let mut x_final = [[0.0_f64; 4]; 4];
        for row in 0..4 {
            let matrix: Vec<Vec<f64>> = m_uu.iter().map(|r| r.to_vec()).collect();
            let b_vec: Vec<f64> = (0..4).map(|c| y[row][c]).collect();
            let Some(x) = solve_dense(matrix, b_vec) else {
                unreachable!("Bernstein Gram M_UU is constant and SPD")
            };
            x_final[row].copy_from_slice(&x[..4]);
        }

        for a in 0..4 {
            for b in 0..4 {
                let i = U_INDICES[a];
                let j = U_INDICES[b];
                let val = x_final[a][b];
                match coord {
                    0 => patch.grid[i][j].x = val,
                    1 => patch.grid[i][j].y = val,
                    _ => patch.grid[i][j].z = val,
                }
            }
        }
    }
}

/// Fit the 36 interior CPs of `patch` to a **bi-nonic** guide.
///
/// Note: even though guide and output are at the same degree, the
/// LSQ is **non-trivial** — the boundary rows/cols of the output
/// come from input data (shared-spoke curves) rather than from the
/// guide, so the interior CPs must compensate for any boundary
/// mismatch. Algorithm is identical to the bicubic/bi-quintic
/// paths; the cross-Gram is just the regular degree-9 Bernstein
/// Gram (degree-9 vs degree-9 is the trivial cross-Gram = gram_9).
pub fn fit_interior_cps_to_binonic_guide(
    patch: &mut RectangularGregoryBiNonic,
    guide: &KpGuidePatchBiNonic,
) {
    let m = bernstein_gram_9();
    let mut m_uu = [[0.0_f64; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            m_uu[a][b] = m[U_INDICES[a]][U_INDICES[b]];
        }
    }
    for coord in 0..3 {
        // Guide CPs for this coordinate (10×10 f64).
        let mut g_coord = [[0.0_f64; 10]; 10];
        for k in 0..10 {
            for l in 0..10 {
                let p = guide.grid[k][l];
                g_coord[k][l] = match coord {
                    0 => p.x,
                    1 => p.y,
                    _ => p.z,
                };
            }
        }
        let mut p_coord = [[0.0_f64; 10]; 10];
        for i in 0..10 {
            for j in 0..10 {
                let p = patch.grid[i][j];
                p_coord[i][j] = match coord {
                    0 => p.x,
                    1 => p.y,
                    _ => p.z,
                };
            }
        }

        // RHS[i*][j*] = Σ_{k, l ∈ 0..=9} G[k][l] · M[k][i*] · M[l][j*].
        let mut rhs = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                let i_star = U_INDICES[a];
                let j_star = U_INDICES[b];
                let mut sum = 0.0_f64;
                for k in 0..10 {
                    for l in 0..10 {
                        sum += g_coord[k][l] * m[k][i_star] * m[l][j_star];
                    }
                }
                rhs[a][b] = sum;
            }
        }

        let mut k_mat = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                let i_star = U_INDICES[a];
                let j_star = U_INDICES[b];
                let mut sum = 0.0_f64;
                for i in 0..10 {
                    for j in 0..10 {
                        if U_INDICES.contains(&i) && U_INDICES.contains(&j) {
                            continue;
                        }
                        sum += p_coord[i][j] * m[i][i_star] * m[j][j_star];
                    }
                }
                k_mat[a][b] = sum;
            }
        }

        let mut eff_rhs = [[0.0_f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                eff_rhs[a][b] = rhs[a][b] - k_mat[a][b];
            }
        }

        let mut y = [[0.0_f64; 4]; 4];
        for col in 0..4 {
            let matrix: Vec<Vec<f64>> = m_uu.iter().map(|row| row.to_vec()).collect();
            let b_vec: Vec<f64> = (0..4).map(|r| eff_rhs[r][col]).collect();
            let Some(x) = solve_dense(matrix, b_vec) else {
                unreachable!("Bernstein Gram M_UU is constant and SPD")
            };
            for r in 0..4 {
                y[r][col] = x[r];
            }
        }
        let mut x_final = [[0.0_f64; 4]; 4];
        for row in 0..4 {
            let matrix: Vec<Vec<f64>> = m_uu.iter().map(|r| r.to_vec()).collect();
            let b_vec: Vec<f64> = (0..4).map(|c| y[row][c]).collect();
            let Some(x) = solve_dense(matrix, b_vec) else {
                unreachable!("Bernstein Gram M_UU is constant and SPD")
            };
            x_final[row].copy_from_slice(&x[..4]);
        }

        for a in 0..4 {
            for b in 0..4 {
                let i = U_INDICES[a];
                let j = U_INDICES[b];
                let val = x_final[a][b];
                match coord {
                    0 => patch.grid[i][j].x = val,
                    1 => patch.grid[i][j].y = val,
                    _ => patch.grid[i][j].z = val,
                }
            }
        }
    }
}

/// 10×10 degree-9 Bernstein Gram matrix.
/// `M[a][b] = C(9, a) · C(9, b) / (19 · C(18, a + b))`.
fn bernstein_gram_9() -> [[f64; 10]; 10] {
    let c9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let c18: [f64; 19] = [
        1.0, 18.0, 153.0, 816.0, 3060.0, 8568.0, 18564.0, 31824.0, 43758.0, 48620.0, 43758.0,
        31824.0, 18564.0, 8568.0, 3060.0, 816.0, 153.0, 18.0, 1.0,
    ];
    let mut m = [[0.0_f64; 10]; 10];
    for a in 0..10 {
        for b in 0..10 {
            m[a][b] = c9[a] * c9[b] / (19.0 * c18[a + b]);
        }
    }
    m
}

/// 4×10 cross-Gram matrix for degree 3 against degree 9.
/// `H[k][i] = C(3, k) · C(9, i) / (13 · C(12, k + i))`.
fn bernstein_cross_gram_3_9() -> [[f64; 10]; 4] {
    let c3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    let c9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let c12: [f64; 13] = [
        1.0, 12.0, 66.0, 220.0, 495.0, 792.0, 924.0, 792.0, 495.0, 220.0, 66.0, 12.0, 1.0,
    ];
    let mut h = [[0.0_f64; 10]; 4];
    for k in 0..4 {
        for i in 0..10 {
            h[k][i] = c3[k] * c9[i] / (13.0 * c12[k + i]);
        }
    }
    h
}

/// 6×10 cross-Gram matrix for degree 5 against degree 9.
/// `H[k][i] = ∫₀¹ B_k^5(t) · B_i^9(t) dt = C(5, k) · C(9, i) / (15 · C(14, k + i))`.
fn bernstein_cross_gram_5_9() -> [[f64; 10]; 6] {
    let c5: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];
    let c9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let c14: [f64; 15] = [
        1.0, 14.0, 91.0, 364.0, 1001.0, 2002.0, 3003.0, 3432.0, 3003.0, 2002.0, 1001.0, 364.0,
        91.0, 14.0, 1.0,
    ];
    let mut h = [[0.0_f64; 10]; 6];
    for k in 0..6 {
        for i in 0..10 {
            h[k][i] = c5[k] * c9[i] / (15.0 * c14[k + i]);
        }
    }
    h
}

#[cfg(test)]
#[path = "kp_interior_fit_tests.rs"]
mod tests;
