//! Dense linear-system solver used by
//! [`super::spoke_poly_g1_solver`].
//!
//! Extracted out of the main solver file to keep it under the
//! 500-line limit. Gaussian elimination with partial pivoting,
//! dedicated to the small (≤ 60 × 60) systems the polynomial-
//! identity G1 solver assembles per spoke. For larger or ill-
//! conditioned systems use a proper library (LAPACK / nalgebra).

/// Solve `m · x = b` via Gaussian elimination with partial
/// pivoting. Returns `None` if the matrix is singular (pivot
/// below `1e-14`).
pub(super) fn solve_linear_system(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = m.len();
    for i in 0..n {
        let mut pivot = i;
        let mut best = m[i][i].abs();
        for (r, row) in m.iter().enumerate().skip(i + 1) {
            let v = row[i].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if pivot != i {
            m.swap(i, pivot);
            b.swap(i, pivot);
        }
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            let pivot_row = m[i].clone();
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                m[r][c] -= factor * pc;
            }
            b[r] -= factor * b[i];
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for (c, xc) in x.iter().enumerate().skip(i + 1) {
            s -= m[i][c] * xc;
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}
