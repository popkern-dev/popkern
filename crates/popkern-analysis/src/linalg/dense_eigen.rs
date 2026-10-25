//! Small dense symmetric eigenproblems by cyclic Jacobi rotations, in a
//! fixed sweep order (deterministic).

use super::dense_chol::Cholesky;

/// Eigenvalues (ascending) and eigenvectors (columns) of a symmetric
/// N×N matrix by cyclic Jacobi rotations, sweeping pairs in fixed order.
pub fn symmetric_eigen<const N: usize>(mut a: [[f64; N]; N]) -> ([f64; N], [[f64; N]; N]) {
    let mut v = [[0.0; N]; N];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..N)
            .flat_map(|p| (0..N).filter(move |&q| q != p).map(move |q| (p, q)))
            .map(|(p, q)| a[p][q] * a[p][q])
            .sum();
        let total: f64 = a.iter().flatten().map(|x| x * x).sum();
        if off <= 1e-30 * total {
            break;
        }
        for p in 0..N {
            for q in p + 1..N {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for row in a.iter_mut() {
                    let (akp, akq) = (row[p], row[q]);
                    row[p] = c * akp - s * akq;
                    row[q] = s * akp + c * akq;
                }
                let (rp, rq) = (a[p], a[q]);
                a[p] = std::array::from_fn(|k| c * rp[k] - s * rq[k]);
                a[q] = std::array::from_fn(|k| s * rp[k] + c * rq[k]);
                for row in &mut v {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..N).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = std::array::from_fn(|k| a[order[k]][order[k]]);
    let vectors = std::array::from_fn(|r| std::array::from_fn(|k| v[r][order[k]]));
    (values, vectors)
}

/// The largest λ with a·x = λ b·x for symmetric n×n `a` and symmetric
/// positive definite `b` (row-major; None when `b` does not factor), by
/// `iterations` steps of power iteration on L⁻¹ a L⁻ᵀ from a fixed start
/// (deterministic): the Rayleigh quotient, a lower estimate that
/// approaches λ from below.
pub fn max_generalized_eigenvalue_dense(
    a: &[f64],
    b: &[f64],
    n: usize,
    iterations: usize,
) -> Option<f64> {
    if n == 0 || a.len() != n * n {
        return None;
    }
    let l = Cholesky::factor(b, n, 0.0).ok()?;
    let l = l.factor_l();
    let forward = |x: &mut [f64]| {
        for i in 0..n {
            let s: f64 = (0..i).map(|k| l[i * n + k] * x[k]).sum();
            x[i] = (x[i] - s) / l[i * n + i];
        }
    };
    let backward = |x: &mut [f64]| {
        for i in (0..n).rev() {
            let s: f64 = (i + 1..n).map(|k| l[k * n + i] * x[k]).sum();
            x[i] = (x[i] - s) / l[i * n + i];
        }
    };
    // y = L⁻¹ a L⁻ᵀ x.
    let apply = |x: &[f64]| {
        let mut z = x.to_vec();
        backward(&mut z);
        let mut y: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| a[i * n + j] * z[j]).sum())
            .collect();
        forward(&mut y);
        y
    };
    let norm = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut x: Vec<f64> = (0..n)
        .map(|i| 1.0 + 0.1 * ((i * 7919) % 13) as f64)
        .collect();
    let mut lambda = 0.0;
    for _ in 0..iterations {
        let nx = norm(&x);
        if nx == 0.0 {
            return Some(0.0);
        }
        x.iter_mut().for_each(|v| *v /= nx);
        let y = apply(&x);
        lambda = x.iter().zip(&y).map(|(p, q)| p * q).sum::<f64>();
        x = y;
    }
    Some(lambda)
}

/// Eigenvalues (ascending) and eigenvectors (columns, row-major n×n) of
/// the symmetric n×n row-major `a`, by cyclic Jacobi rotations in a
/// fixed sweep order: [`symmetric_eigen`] for a size known at run time.
pub fn symmetric_eigen_dense(mut a: Vec<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in 0..n {
                if p != q {
                    off += a[p * n + q] * a[p * n + q];
                }
            }
        }
        let total: f64 = a.iter().map(|x| x * x).sum();
        if off <= 1e-30 * total {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i * n + i].total_cmp(&a[j * n + j]));
    let values = order.iter().map(|&i| a[i * n + i]).collect();
    let mut vectors = vec![0.0; n * n];
    for (col, &i) in order.iter().enumerate() {
        for r in 0..n {
            vectors[r * n + col] = v[r * n + i];
        }
    }
    (values, vectors)
}
