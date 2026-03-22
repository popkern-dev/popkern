//! Linear algebra helpers for surface-surface intersection marching.

/// Solve a 4×4 linear system `A * x = b` via Gaussian elimination with partial pivoting.
///
/// Returns `None` if the matrix is singular (pivot smaller than `1e-30`).
#[allow(clippy::needless_range_loop)]
pub(super) fn solve_4x4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    // Forward elimination with partial pivoting.
    for col in 0..4 {
        // Find the row with the largest absolute value in this column.
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for row in (col + 1)..4 {
            if a[row][col].abs() > max_val {
                max_val = a[row][col].abs();
                max_row = row;
            }
        }

        if max_val < 1e-30 {
            return None;
        }

        // Swap the pivot row into position.
        if max_row != col {
            a.swap(col, max_row);
            b.swap(col, max_row);
        }

        // Eliminate entries below the pivot.
        for row in (col + 1)..4 {
            let factor = a[row][col] / a[col][col];
            for k in col..4 {
                a[row][k] -= factor * a[col][k];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution.
    let mut x = [0.0_f64; 4];
    for i in (0..4).rev() {
        x[i] = b[i];
        for j in (i + 1)..4 {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Basic sanity check: solving against the identity matrix returns the
    /// right-hand side unchanged.
    #[test]
    fn solve_4x4_identity() {
        let a = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let b = [1.0, 2.0, 3.0, 4.0];
        let x = solve_4x4(a, b).expect("identity system must be solvable");
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 2.0).abs() < 1e-10);
        assert!((x[2] - 3.0).abs() < 1e-10);
        assert!((x[3] - 4.0).abs() < 1e-10);
    }

    /// A truly singular matrix must return None.
    #[test]
    fn solve_4x4_singular_returns_none() {
        let a = [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 8.0], // row 1 = 2 * row 0
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let b = [1.0, 2.0, 3.0, 4.0];
        assert!(solve_4x4(a, b).is_none());
    }

    /// Verify that the solver produces a correct answer for a known 4×4 system.
    #[test]
    fn solve_4x4_known_solution() {
        // System: x0=1, x1=2, x2=3, x3=4 encoded as a non-trivial matrix.
        let a = [
            [2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 3.0],
        ];
        // b = A * [1, 2, 3, 4]^T
        let b = [
            2.0 * 1.0 + 1.0 * 2.0,
            1.0 * 1.0 + 3.0 * 2.0 + 1.0 * 3.0,
            1.0 * 2.0 + 4.0 * 3.0 + 1.0 * 4.0,
            1.0 * 3.0 + 3.0 * 4.0,
        ];
        let x = solve_4x4(a, b).expect("system must be solvable");
        assert!((x[0] - 1.0).abs() < 1e-10, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1] = {}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-10, "x[2] = {}", x[2]);
        assert!((x[3] - 4.0).abs() < 1e-10, "x[3] = {}", x[3]);
    }
}
