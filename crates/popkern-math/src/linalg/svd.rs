//! One-sided Jacobi singular value decomposition (Hestenes), moved here
//! from the variational cap fit so every caller shares one
//! deterministic SVD.
//!
//! The fit's systems are exactly rank deficient: spoke blocks and
//! seam-relation blocks carry many exact zero singular values. faer
//! 0.20's SVD returned factors that do not reproduce such matrices
//! (measured: a 32 × 20 rank-8 seam-relation block reconstructed with
//! error 2.65). One-sided Jacobi orthogonalises the columns of `A V` by
//! plane rotations until every pair is orthogonal to working precision;
//! zero columns stay zero, so rank deficiency costs nothing, and the
//! singular values come out with high relative accuracy. It is a fixed
//! sequence of rotations, so the result depends only on the input.

use super::DenseMat;

/// Sweeps before the iteration gives up (each sweep visits every column
/// pair; convergence is quadratic, typically 6-10 sweeps).
const MAX_SWEEPS: usize = 60;

/// `A V = W`: `sigma` the column norms of `W` (nonincreasing, one per
/// column of `A`), `u` the columns of `W` normalised (a zero column left
/// zero), `v` square and orthogonal.
#[derive(Debug, Clone)]
pub struct Svd {
    /// Singular values, nonincreasing.
    pub sigma: Vec<f64>,
    /// `rows × cols` left vectors (zero for zero singular values).
    pub u: DenseMat,
    /// `cols × cols` right vectors.
    pub v: DenseMat,
}

/// The SVD of `a`; `None` when the sweeps do not converge. A wide `a`
/// leaves `cols − rows` singular values at zero, whose `v` columns span
/// its null space.
pub fn jacobi_svd(a: &DenseMat) -> Option<Svd> {
    let (m, n) = (a.nrows(), a.ncols());
    let mut w: Vec<Vec<f64>> = (0..n).map(|j| a.col(j).to_vec()).collect();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|j| (0..n).map(|i| f64::from(u8::from(i == j))).collect())
        .collect();
    let eps = f64::EPSILON;
    // A column at rounding level of the whole matrix is a zero singular
    // value: it has no direction to orthogonalise, so it takes no more
    // rotations.
    let frobenius2: f64 = w.iter().flatten().map(|x| x * x).sum();
    let negligible = (eps * eps) * frobenius2 * (m.max(n) as f64);
    let mut converged = n < 2;
    for _ in 0..MAX_SWEEPS {
        let mut rotated = false;
        for p in 0..n {
            for q in p + 1..n {
                let (alpha, beta, gamma) = {
                    let (wp, wq) = (&w[p], &w[q]);
                    let (mut alpha, mut beta, mut gamma) = (0.0, 0.0, 0.0);
                    for i in 0..m {
                        alpha += wp[i] * wp[i];
                        beta += wq[i] * wq[i];
                        gamma += wp[i] * wq[i];
                    }
                    (alpha, beta, gamma)
                };
                if gamma == 0.0
                    || alpha <= negligible
                    || beta <= negligible
                    || gamma.abs() <= eps * (alpha * beta).sqrt()
                {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                rotate(&mut w, p, q, c, s);
                rotate(&mut v, p, q, c, s);
            }
        }
        if !rotated {
            converged = true;
            break;
        }
    }
    if !converged {
        return None;
    }
    let norms: Vec<f64> = w
        .iter()
        .map(|col| col.iter().map(|x| x * x).sum::<f64>().sqrt())
        .collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| norms[y].total_cmp(&norms[x]));
    let sigma: Vec<f64> = order.iter().map(|&j| norms[j]).collect();
    let u = DenseMat::from_fn(m, n, |i, k| {
        let j = order[k];
        if norms[j] > 0.0 {
            w[j][i] / norms[j]
        } else {
            0.0
        }
    });
    let v = DenseMat::from_fn(n, n, |i, k| v[order[k]][i]);
    Some(Svd { sigma, u, v })
}

/// Rotate columns `p` and `q`: `(c x_p − s x_q, s x_p + c x_q)`.
fn rotate(cols: &mut [Vec<f64>], p: usize, q: usize, c: f64, s: f64) {
    let (head, tail) = cols.split_at_mut(q);
    let (xp, xq) = (&mut head[p], &mut tail[0]);
    for (a, b) in xp.iter_mut().zip(xq.iter_mut()) {
        let (x, y) = (*a, *b);
        *a = c * x - s * y;
        *b = s * x + c * y;
    }
}
