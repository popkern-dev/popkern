//! Tests for the shared linear-algebra helpers.

use super::*;

#[test]
fn solve_dense_2x2_identity_returns_rhs() {
    let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let b = vec![3.0, -2.0];
    let x = solve_dense(m, b).expect("solved");
    assert!((x[0] - 3.0).abs() < 1e-12);
    assert!((x[1] - (-2.0)).abs() < 1e-12);
}

#[test]
fn solve_dense_3x3_non_trivial() {
    // Use a system with known solution:
    //   x + 2y + z = 4
    //   3x - y + 2z = 9
    //   2x + y - z = 1
    // Solution: (x, y, z) = (2, 1, 0).
    let m = vec![
        vec![1.0, 2.0, 1.0],
        vec![3.0, -1.0, 2.0],
        vec![2.0, 1.0, -1.0],
    ];
    let b = vec![4.0, 9.0, 1.0];
    let x = solve_dense(m, b).expect("solved");
    // Verify by substitution rather than expected values.
    let r0 = 1.0 * x[0] + 2.0 * x[1] + 1.0 * x[2];
    let r1 = 3.0 * x[0] - 1.0 * x[1] + 2.0 * x[2];
    let r2 = 2.0 * x[0] + 1.0 * x[1] - 1.0 * x[2];
    assert!((r0 - 4.0).abs() < 1e-10);
    assert!((r1 - 9.0).abs() < 1e-10);
    assert!((r2 - 1.0).abs() < 1e-10);
}

#[test]
fn solve_dense_singular_returns_none() {
    let m = vec![vec![1.0, 2.0], vec![2.0, 4.0]]; // row 2 = 2 * row 1
    let b = vec![3.0, 6.0];
    assert!(solve_dense(m, b).is_none());
}

#[test]
fn orthonormal_basis_is_orthonormal_for_arbitrary_normal() {
    for (nx, ny, nz) in [(0.3, -0.4, 0.5), (0.0, 1.0, 0.0), (0.1, 0.1, 0.995)] {
        let n = Vec3::new(nx, ny, nz).normalized().unwrap();
        let (t1, t2) = orthonormal_tangent_basis(n);
        assert!((t1.length() - 1.0).abs() < 1e-10);
        assert!((t2.length() - 1.0).abs() < 1e-10);
        assert!(t1.dot(&t2).abs() < 1e-10);
        assert!(t1.dot(&n).abs() < 1e-10);
        assert!(t2.dot(&n).abs() < 1e-10);
    }
}

#[test]
fn binomial_values_match_pascal() {
    // First 5 rows of Pascal's triangle.
    let expected = [
        vec![1.0],
        vec![1.0, 1.0],
        vec![1.0, 2.0, 1.0],
        vec![1.0, 3.0, 3.0, 1.0],
        vec![1.0, 4.0, 6.0, 4.0, 1.0],
    ];
    for (n, row) in expected.iter().enumerate() {
        for (k, &expected_val) in row.iter().enumerate() {
            assert!((binomial(n, k) - expected_val).abs() < 1e-12, "C({n}, {k})");
        }
    }
}

#[test]
fn binomial_k_greater_than_n_is_zero() {
    assert_eq!(binomial(3, 5), 0.0);
    assert_eq!(binomial(7, 8), 0.0);
}
