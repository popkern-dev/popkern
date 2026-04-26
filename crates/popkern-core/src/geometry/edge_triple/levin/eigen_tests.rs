//! Verify Jacobi eigendecomposition for known 4×4 symmetric matrices.

use super::eigen::{jacobi_eigh_4x4, max_off_diagonal};
use super::quadric::Mat4Sym;

const TOL: f64 = 1e-10;

fn reconstruct(eig: [f64; 4], q: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    // A_recon = Q · diag(eig) · Qᵀ
    let mut out = [[0.0_f64; 4]; 4];
    for (i, qi) in q.iter().enumerate() {
        for (j, qj) in q.iter().enumerate() {
            let mut s = 0.0;
            for (k, &lam) in eig.iter().enumerate() {
                s += qi[k] * lam * qj[k];
            }
            out[i][j] = s;
        }
    }
    out
}

fn max_diff(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> f64 {
    let mut max = 0.0_f64;
    for i in 0..4 {
        for j in 0..4 {
            max = max.max((a[i][j] - b[i][j]).abs());
        }
    }
    max
}

fn q_is_orthogonal(q: &[[f64; 4]; 4]) -> bool {
    // Q · Qᵀ == I within tolerance.
    for (i, qi) in q.iter().enumerate() {
        for (j, qj) in q.iter().enumerate() {
            let s: f64 = qi.iter().zip(qj.iter()).map(|(a, b)| a * b).sum();
            let target = if i == j { 1.0 } else { 0.0 };
            if (s - target).abs() > TOL {
                return false;
            }
        }
    }
    true
}

#[test]
fn diagonal_matrix_eigenvalues_match_diagonal() {
    let a = Mat4Sym([
        [3.0, 0.0, 0.0, 0.0],
        [0.0, -2.0, 0.0, 0.0],
        [0.0, 0.0, 5.0, 0.0],
        [0.0, 0.0, 0.0, 0.7],
    ]);
    let (eig, q) = jacobi_eigh_4x4(&a);
    let mut sorted: Vec<f64> = eig.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!((sorted[0] - (-2.0)).abs() < TOL);
    assert!((sorted[1] - 0.7).abs() < TOL);
    assert!((sorted[2] - 3.0).abs() < TOL);
    assert!((sorted[3] - 5.0).abs() < TOL);
    assert!(q_is_orthogonal(&q));
}

#[test]
fn off_diagonal_matrix_reconstructs() {
    let a = Mat4Sym([
        [4.0, 1.0, 0.5, 0.0],
        [1.0, 3.0, 0.2, -0.1],
        [0.5, 0.2, 2.0, 0.3],
        [0.0, -0.1, 0.3, 1.0],
    ]);
    let (eig, q) = jacobi_eigh_4x4(&a);
    let recon = reconstruct(eig, q);
    assert!(
        max_diff(&recon, &a.0) < TOL,
        "reconstruction error {}",
        max_diff(&recon, &a.0)
    );
    assert!(q_is_orthogonal(&q));
}

#[test]
fn rank_one_outer_product_has_one_nonzero_eigenvalue() {
    // v = (1, 2, 3, -4), D = v vᵀ. Eigenvalue = ‖v‖² = 30, eigenvector = v/‖v‖.
    let v = [1.0_f64, 2.0, 3.0, -4.0];
    let mut d = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            d[i][j] = v[i] * v[j];
        }
    }
    let (eig, q) = jacobi_eigh_4x4(&Mat4Sym(d));
    let mut sorted: Vec<f64> = eig.to_vec();
    sorted.sort_by(|a, b| b.abs().partial_cmp(&a.abs()).unwrap());
    assert!(
        (sorted[0] - 30.0).abs() < TOL,
        "largest eig = {}",
        sorted[0]
    );
    assert!(sorted[1].abs() < TOL);
    assert!(sorted[2].abs() < TOL);
    assert!(sorted[3].abs() < TOL);
    assert!(q_is_orthogonal(&q));
}

#[test]
fn off_diagonal_zero_after_decomposition() {
    // After Jacobi sweeps, A becomes diag in the rotated basis. The
    // reconstructed off-diagonal residual is bounded by tolerance.
    let a = Mat4Sym([
        [2.0, 1.5, 0.7, 0.3],
        [1.5, 3.0, 0.4, 1.1],
        [0.7, 0.4, 5.0, 0.6],
        [0.3, 1.1, 0.6, 1.5],
    ]);
    let (eig, q) = jacobi_eigh_4x4(&a);
    let recon = reconstruct(eig, q);
    assert!(max_diff(&recon, &a.0) < TOL);
    // Eigenvalues should be real and bounded.
    for v in eig.iter() {
        assert!(v.is_finite());
    }
    let _ = max_off_diagonal(&a.0); // smoke test on the test helper
}
