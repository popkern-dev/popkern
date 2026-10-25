use super::dense_chol::{Cholesky, CholeskyError};
use super::spd_test_helpers::{Lcg, bits_hash, frobenius, matmul, spectrum, with_spectrum};

const SIZES: [usize; 6] = [1, 2, 3, 8, 24, 33];
const CONDITIONS: [f64; 4] = [1.0, 1e2, 1e5, 1e10];

fn factor(a: &[f64], n: usize) -> Cholesky {
    Cholesky::factor(a, n, 0.0).unwrap_or_else(|e| panic!("n {n}: {e}"))
}

#[test]
fn the_factor_reconstructs_the_matrix() {
    let mut worst = 0.0f64;
    for (s, &n) in SIZES.iter().enumerate() {
        for (c, &kappa) in CONDITIONS.iter().enumerate() {
            let a = with_spectrum(n, &spectrum(n, kappa), (s * 10 + c) as u64 + 1);
            let l = factor(&a, n);
            let mut lt = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    lt[j * n + i] = l.factor_l()[i * n + j];
                }
            }
            let llt = matmul(l.factor_l(), &lt, n);
            let diff: Vec<f64> = llt.iter().zip(&a).map(|(x, y)| x - y).collect();
            let ratio = frobenius(&diff) / frobenius(&a) / (n as f64 * f64::EPSILON);
            assert!(
                ratio <= 4.0,
                "n {n} kappa {kappa}: ||LLt - A|| / ||A|| = {ratio} n eps"
            );
            worst = worst.max(ratio);
        }
    }
    // Measured worst case 0.354 n eps (2026-10-10); a regression past
    // 1 n eps shows here before it reaches the 4 n eps bound.
    assert!(worst < 1.0, "worst reconstruction {worst} n eps");
}

#[test]
fn the_solve_recovers_the_known_solution() {
    let mut worst = 0.0f64;
    for (s, &n) in SIZES.iter().enumerate() {
        for &kappa in &CONDITIONS {
            let a = with_spectrum(n, &spectrum(n, kappa), 100 + s as u64);
            let mut rng = Lcg(7 + n as u64);
            let x_true: Vec<f64> = (0..n).map(|_| rng.next_signed()).collect();
            let b: Vec<f64> = (0..n)
                .map(|i| (0..n).map(|j| a[i * n + j] * x_true[j]).sum())
                .collect();
            let x = factor(&a, n).solve(&b).unwrap_or_else(|e| panic!("{e}"));
            let err: Vec<f64> = x.iter().zip(&x_true).map(|(p, q)| p - q).collect();
            let rel = frobenius(&err) / frobenius(&x_true);
            let ratio = rel / (kappa * n as f64 * f64::EPSILON);
            assert!(
                ratio <= 1.0,
                "n {n} kappa {kappa}: error {rel} = {ratio} kappa n eps"
            );
            worst = worst.max(ratio);
        }
    }
    // Measured worst case 0.424 kappa n eps (2026-10-10). Non-zero error
    // shows the check is not vacuous.
    assert!(worst > 0.0, "the bound check must see some error");
}

#[test]
fn indefinite_singular_and_nan_matrices_fail_loudly() {
    assert_eq!(
        Cholesky::factor(&[1.0, 2.0, 2.0, 1.0], 2, 0.0),
        Err(CholeskyError::NonPositivePivot {
            index: 1,
            value: -3.0,
            diagonal: 1.0
        })
    );
    assert_eq!(
        Cholesky::factor(&[1.0, 1.0, 1.0, 1.0], 2, 0.0),
        Err(CholeskyError::NonPositivePivot {
            index: 1,
            value: 0.0,
            diagonal: 1.0
        })
    );
    assert_eq!(
        Cholesky::factor(&[-2.0], 1, 0.0),
        Err(CholeskyError::NonPositivePivot {
            index: 0,
            value: -2.0,
            diagonal: -2.0
        })
    );
    assert_eq!(
        Cholesky::factor(&[4.0, 0.0, f64::NAN, 1.0], 2, 0.0),
        Err(CholeskyError::NonFinite { row: 1, col: 0 })
    );
    // An indefinite matrix with a generated eigenbasis fails at a
    // pivot, whichever one it is.
    let n = 8;
    let mut lambda = spectrum(n, 1e3);
    lambda[5] = -0.5;
    let a = with_spectrum(n, &lambda, 42);
    assert!(matches!(
        Cholesky::factor(&a, n, 0.0),
        Err(CholeskyError::NonPositivePivot { .. })
    ));
}

#[test]
fn a_pivot_ratio_rejects_numerically_singular_matrices() {
    let n = 24;
    let mut lambda = spectrum(n, 1e3);
    lambda[n - 1] = 0.0;
    let a = with_spectrum(n, &lambda, 9);
    let ratio = n as f64 * f64::EPSILON;
    let Err(CholeskyError::NonPositivePivot {
        value, diagonal, ..
    }) = Cholesky::factor(&a, n, ratio)
    else {
        panic!("a rank-deficient matrix factored with ratio {ratio}");
    };
    assert!(value <= ratio * diagonal);
    // The same matrix at a healthy spectrum passes with that ratio.
    let healthy = with_spectrum(n, &spectrum(n, 1e3), 9);
    assert!(Cholesky::factor(&healthy, n, ratio).is_ok());
}

#[test]
fn shape_and_length_mismatches_are_typed() {
    assert_eq!(
        Cholesky::factor(&[1.0, 0.0, 1.0], 2, 0.0),
        Err(CholeskyError::Shape { len: 3, n: 2 })
    );
    let l = factor(&[4.0, 0.0, 0.0, 9.0], 2);
    assert_eq!(
        l.solve(&[1.0]),
        Err(CholeskyError::RhsLength { len: 1, n: 2 })
    );
    assert_eq!(l.solve(&[4.0, 9.0]), Ok(vec![1.0, 1.0]));
    assert_eq!(l.factor_l(), &[2.0, 0.0, 0.0, 3.0]);
    assert_eq!(l.order(), 2);
}

/// The documented loop order, written out independently.
fn reference_factor(a: &[f64], n: usize) -> Vec<f64> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = if i == j { s.sqrt() } else { s / l[j * n + j] };
        }
    }
    l
}

#[test]
fn the_factor_bits_are_fixed() {
    let n = 33;
    let a = with_spectrum(n, &spectrum(n, 1e5), 2024);
    let first = factor(&a, n);
    let second = factor(&a, n);
    let first_bits: Vec<u64> = first.factor_l().iter().map(|x| x.to_bits()).collect();
    let second_bits: Vec<u64> = second.factor_l().iter().map(|x| x.to_bits()).collect();
    assert_eq!(first_bits, second_bits);
    let reference: Vec<u64> = reference_factor(&a, n)
        .iter()
        .map(|x| x.to_bits())
        .collect();
    assert_eq!(first_bits, reference);
    assert_eq!(bits_hash(first.factor_l()), PINNED_FACTOR_HASH);
}

/// FNV-1a of L's bits for the n = 33, κ = 1e5, seed 2024 matrix. A
/// change means the factor's arithmetic order changed.
const PINNED_FACTOR_HASH: u64 = 9605497310090428802;
