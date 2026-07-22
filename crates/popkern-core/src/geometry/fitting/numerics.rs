//! Numerical kernels for B-spline fitting: Cox–de Boor basis
//! evaluation and the dense linear solve.
//!
//! Split from `fitting.rs` to honour the 500-line structural cap.

use crate::error::CoreError;

/// Compute the i-th B-spline basis function of given degree at parameter t.
///
/// Uses the de Boor-Cox recursion formula. Zero denominators yield 0/0 = 0.
pub(crate) fn basis_function(i: usize, degree: usize, t: f64, knots: &[f64]) -> f64 {
    if degree == 0 {
        // For the last knot span, use <= for the right endpoint so that
        // t == knots[last] is captured.
        if knots[i] <= t && t < knots[i + 1] {
            return 1.0;
        }
        // Special case: t == knots[i+1] and this is the last valid span
        if (t - knots[i + 1]).abs() < 1e-14 && knots[i] < knots[i + 1] {
            // Check if there are no more non-degenerate spans after this
            let mut last_span = true;
            for j in (i + 1)..knots.len().saturating_sub(1) {
                if knots[j] < knots[j + 1] {
                    last_span = false;
                    break;
                }
            }
            if last_span {
                return 1.0;
            }
        }
        return 0.0;
    }

    let mut result = 0.0;

    let denom1 = knots[i + degree] - knots[i];
    if denom1.abs() > 1e-14 {
        result += ((t - knots[i]) / denom1) * basis_function(i, degree - 1, t, knots);
    }

    let denom2 = knots[i + degree + 1] - knots[i + 1];
    if denom2.abs() > 1e-14 {
        result +=
            ((knots[i + degree + 1] - t) / denom2) * basis_function(i + 1, degree - 1, t, knots);
    }

    result
}

/// Solve A*x = b for a symmetric positive-definite matrix A using Gaussian
/// elimination with partial pivoting. A is n x n, b is n x 1.
/// Returns the solution vector x, or an error if the system is singular.
#[allow(clippy::needless_range_loop)]
pub(crate) fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, CoreError> {
    let n = a.len();
    if n == 0 {
        return Ok(vec![]);
    }

    // Build augmented matrix
    let mut aug: Vec<Vec<f64>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = a[i].clone();
        row.push(b[i]);
        aug.push(row);
    }

    // Forward elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let mut max_val = aug[col][col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let val = aug[row][col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }

        if max_val < 1e-14 {
            return Err(CoreError::PreconditionFailed(
                "singular matrix in curve fitting".to_string(),
            ));
        }

        if max_row != col {
            aug.swap(col, max_row);
        }

        let pivot = aug[col][col];
        for row in (col + 1)..n {
            let factor = aug[row][col] / pivot;
            for j in col..=n {
                let val = aug[col][j];
                aug[row][j] -= factor * val;
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum -= aug[i][j] * x[j];
        }
        x[i] = sum / aug[i][i];
    }

    Ok(x)
}
