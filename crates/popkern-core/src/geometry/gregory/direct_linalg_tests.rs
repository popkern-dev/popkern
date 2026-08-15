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

/// Frozen pre-PK-289 clone-based elimination loop: the
/// arithmetic-order oracle for `solve_dense`'s bit-identity
/// contract (see the module docs). Do not "fix" or modernize this
/// body — its value is that it performs exactly the operations the
/// original solver performed, in the same order.
fn solve_dense_reference(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
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

/// Deterministic linear congruential generator mapped to [-1, 1) so
/// the oracle systems are reproducible without external randomness.
fn lcg_next(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
}

#[test]
fn solve_dense_is_bit_identical_to_the_clone_based_reference() {
    let mut state = 0x5eed_u64;
    for n in [1usize, 2, 3, 7, 24, 61, 128] {
        // Dense systems with magnitudes spread across rows so the
        // partial-pivot search actually swaps rows, exercising the
        // same elimination paths the KKT solves take.
        let mut m = vec![vec![0.0_f64; n]; n];
        for row in m.iter_mut() {
            for entry in row.iter_mut() {
                *entry = lcg_next(&mut state);
            }
        }
        for (r, row) in m.iter_mut().enumerate() {
            let scale = if r % 3 == 0 { 100.0 } else { 1.0 };
            for entry in row.iter_mut() {
                *entry *= scale;
            }
        }
        let b: Vec<f64> = (0..n).map(|_| lcg_next(&mut state)).collect();
        let reference = solve_dense_reference(m.clone(), b.clone());
        let solved = solve_dense(m, b);
        match (reference, solved) {
            (Some(reference), Some(solved)) => {
                assert_eq!(reference.len(), solved.len(), "n={n}");
                for (i, (r, s)) in reference.iter().zip(&solved).enumerate() {
                    assert_eq!(
                        r.to_bits(),
                        s.to_bits(),
                        "n={n} x[{i}]: reference {r:e} vs solved {s:e}"
                    );
                }
            }
            (None, None) => {}
            (reference, solved) => {
                panic!("solvability diverged at n={n}: reference={reference:?} solved={solved:?}")
            }
        }
    }
}

#[test]
fn finite_guarded_solve_matches_reference_on_well_posed_system() {
    let m = vec![
        vec![1.0, 2.0, 1.0],
        vec![3.0, -1.0, 2.0],
        vec![2.0, 1.0, -1.0],
    ];
    let b = vec![4.0, 9.0, 1.0];
    let x = solve_dense_finite_guarded(m.clone(), b.clone()).expect("solved");
    // Verify by substitution.
    for (row, &target) in m.iter().zip(&b) {
        let value: f64 = row.iter().zip(&x).map(|(&c, &xi)| c * xi).sum();
        assert!((value - target).abs() < 1e-10);
    }
}

#[test]
fn finite_guarded_solve_rejects_non_finite_input() {
    let m = vec![vec![1.0, f64::NAN], vec![0.0, 1.0]];
    let b = vec![1.0, 1.0];
    assert!(solve_dense_finite_guarded(m, b).is_none());
    let m = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let b = vec![1.0, f64::INFINITY];
    assert!(solve_dense_finite_guarded(m, b).is_none());
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
