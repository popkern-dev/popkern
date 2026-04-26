//! 4×4 symmetric eigendecomposition via cyclic Jacobi rotations.
//!
//! Used by `decompose.rs` to determine the rank of a degenerate
//! quadric fiber from the pencil. Jacobi is overkill for a 4×4
//! matrix in raw flop terms but it's a few-dozen-line robust
//! algorithm with no exotic edge cases.
//!
//! Output: eigenvalues in `eigenvalues[k]` paired with the unit
//! eigenvector forming the `k`-th column of `q`. `q` is orthogonal
//! to ~machine precision.
//!
//! Convergence: cyclic sweeps until the sum of squared off-diagonal
//! entries falls below `TOL`. Bounded at `MAX_SWEEPS = 50` iterations
//! (4×4 typically converges in 5–8 sweeps).

use super::quadric::Mat4Sym;

const TOL_OFF: f64 = 1e-28;
const MAX_SWEEPS: usize = 50;

/// Eigendecomposition of a 4×4 real symmetric matrix.
///
/// Returns `(eigenvalues, q)` where `eigenvalues[i]` corresponds to
/// the unit eigenvector `q[*][i]` (column `i`). Order is the
/// natural Jacobi-sweep order — *not* sorted by magnitude. Callers
/// that need sorted output should sort after.
pub(crate) fn jacobi_eigh_4x4(input: &Mat4Sym) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut a = input.0;
    let mut q = [[0.0_f64; 4]; 4];
    for (i, row) in q.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    for _sweep in 0..MAX_SWEEPS {
        let mut off = 0.0_f64;
        for (i, row) in a.iter().enumerate() {
            for &val in &row[(i + 1)..] {
                off += val * val;
            }
        }
        if off < TOL_OFF {
            break;
        }
        // Cyclic sweep over (p, r) pairs with p < r.
        for p in 0..3 {
            for r in (p + 1)..4 {
                let app = a[p][p];
                let arr = a[r][r];
                let apr = a[p][r];
                if apr.abs() < 1e-300 {
                    continue;
                }
                // Compute the rotation tangent stably.
                let theta = (arr - app) / (2.0 * apr);
                let t = if theta.abs() > 1e10 {
                    0.5 / theta
                } else {
                    let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                    sign / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                // Update diagonal.
                a[p][p] = app - t * apr;
                a[r][r] = arr + t * apr;
                a[p][r] = 0.0;
                a[r][p] = 0.0;
                // Update other rows/cols affected by the (p, r) rotation.
                for i in (0..4).filter(|i| *i != p && *i != r) {
                    let aip = a[i][p];
                    let air = a[i][r];
                    a[i][p] = c * aip - s * air;
                    a[i][r] = s * aip + c * air;
                    a[p][i] = a[i][p];
                    a[r][i] = a[i][r];
                }
                // Update accumulated rotation Q.
                for row in q.iter_mut() {
                    let qp = row[p];
                    let qr = row[r];
                    row[p] = c * qp - s * qr;
                    row[r] = s * qp + c * qr;
                }
            }
        }
    }

    let eigenvalues = [a[0][0], a[1][1], a[2][2], a[3][3]];
    (eigenvalues, q)
}

/// Test-only: maximum off-diagonal magnitude in the 4×4 matrix.
#[cfg(test)]
pub(crate) fn max_off_diagonal(m: &[[f64; 4]; 4]) -> f64 {
    let mut max = 0.0_f64;
    for (i, row) in m.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if i != j {
                max = max.max(v.abs());
            }
        }
    }
    max
}
