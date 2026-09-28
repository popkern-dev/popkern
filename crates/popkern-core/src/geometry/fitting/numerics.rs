//! Numerical kernels for B-spline fitting: Cox–de Boor basis
//! evaluation and the banded linear solve (with its dense reference).
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
///
/// The dense reference for [`solve_banded_system`], kept to pin the
/// banded path's bit-identity (production fits are banded, PK-572).
#[cfg(test)]
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
            // Exact-zero skip: banded fitting systems (PK-329) leave most
            // sub-pivot entries structurally 0.0; eliminating them applies
            // `x -= ±0.0` across the row, which is bit-neutral here (no
            // entry in these chains is ever -0.0), so skipping is
            // bit-identical and turns the dense O(n³) elimination into
            // O(n²·bandwidth) on banded inputs.
            if aug[row][col] == 0.0 {
                continue;
            }
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

/// One row of a banded square matrix: the stored window starts at
/// column `start` and holds `vals.len()` consecutive entries; every
/// column outside the window is exactly `0.0`.
#[derive(Debug, Clone)]
pub(crate) struct BandRow {
    /// First stored column.
    pub(crate) start: usize,
    /// Stored entries for columns `start..start + vals.len()`.
    pub(crate) vals: Vec<f64>,
}

impl BandRow {
    fn get(&self, col: usize) -> f64 {
        if col < self.start {
            return 0.0;
        }
        self.vals.get(col - self.start).copied().unwrap_or(0.0)
    }

    fn end(&self) -> usize {
        self.start + self.vals.len()
    }
}

/// Solve `A·X = B` for a banded `A` (row windows) and several
/// right-hand sides at once, by Gaussian elimination with partial
/// pivoting.
///
/// Memory is O(n·bandwidth) and work O(n·bandwidth²), where the dense
/// [`solve_linear_system`] costs O(n²) memory: a 198k-sample marched
/// SSI trace fitted near full interpolation needed a 12.8 GB normal
/// matrix per copy (PK-572, OOM-killed in the test sweep).
///
/// Bit-identical to the dense solver run once per right-hand side:
/// the pivot search, row swaps, eliminations, and back substitution
/// perform the same floating-point operations in the same order and
/// skip only terms the dense path applies as `x -= ±0.0` (bit-neutral
/// on these chains, which never hold `-0.0`). The candidate rows at
/// column `c` are those within the lower bandwidth `w` below it: a
/// row whose window starts at or before `c` either never moved (its
/// original index is ≤ `c + w`) or was swapped down to a pivot
/// candidate position (≤ `c' + w` for an earlier `c' < c`).
/// Elimination fill extends a row's window rightward only — the
/// standard LU-with-pivoting upper-band growth to `2w`.
#[allow(clippy::needless_range_loop)]
pub(crate) fn solve_banded_system(
    mut rows: Vec<BandRow>,
    rhs: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, CoreError> {
    let n = rows.len();
    let k = rhs.len();
    if rhs.iter().any(|b| b.len() != n) {
        return Err(CoreError::PreconditionFailed(
            "banded solve: right-hand side length mismatch".to_string(),
        ));
    }
    if n == 0 {
        return Ok(vec![Vec::new(); k]);
    }
    // Row-major right-hand sides so swaps move a row's values together.
    let mut b: Vec<Vec<f64>> = (0..n).map(|i| rhs.iter().map(|c| c[i]).collect()).collect();
    let w = rows
        .iter()
        .enumerate()
        .map(|(i, r)| i.saturating_sub(r.start))
        .max()
        .unwrap_or(0);

    for col in 0..n {
        let last = (col + w).min(n - 1);
        let mut max_val = rows[col].get(col).abs();
        let mut max_row = col;
        for row in (col + 1)..=last {
            let val = rows[row].get(col).abs();
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
            rows.swap(col, max_row);
            b.swap(col, max_row);
        }
        let (head, tail) = rows.split_at_mut(col + 1);
        let pivot_row = &head[col];
        let pivot = pivot_row.get(col);
        let pivot_end = pivot_row.end();
        let (b_head, b_tail) = b.split_at_mut(col + 1);
        let b_pivot = &b_head[col];
        for (offset, row) in tail.iter_mut().enumerate().take(last - col) {
            let a = row.get(col);
            if a == 0.0 {
                continue;
            }
            let factor = a / pivot;
            if row.end() < pivot_end {
                row.vals.resize(pivot_end - row.start, 0.0);
            }
            for j in col..pivot_end {
                let val = pivot_row.vals[j - pivot_row.start];
                row.vals[j - row.start] -= factor * val;
            }
            let b_row = &mut b_tail[offset];
            for c in 0..k {
                b_row[c] -= factor * b_pivot[c];
            }
        }
    }

    let mut x = vec![vec![0.0; n]; k];
    for i in (0..n).rev() {
        let row = &rows[i];
        let diag = row.get(i);
        for c in 0..k {
            let mut sum = b[i][c];
            for j in (i + 1)..row.end() {
                sum -= row.vals[j - row.start] * x[c][j];
            }
            x[c][i] = sum / diag;
        }
    }
    Ok(x)
}
