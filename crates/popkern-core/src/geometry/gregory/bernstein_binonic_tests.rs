//! Tests for [`super`] degree-9 Bernstein basis + derivatives.
//!
//! Philosophy: cross-check the specialised implementations
//! against the generic [`super::super::super::bernstein_poly::bernstein_basis`]
//! and against finite-difference references. This catches
//! transcription errors (binomial table off-by-one, factor-72
//! mistake, etc.) with machine-precision comparisons.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::bernstein_poly::{BernsteinPoly, bernstein_basis};

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

#[test]
fn partition_of_unity_at_many_ts() {
    // Σ_i B_i^9(t) ≡ 1 for every t.
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=9).map(|i| bernstein_binonic(i, t)).sum();
        assert!((sum - 1.0).abs() < 1e-10, "t = {t}, sum = {sum}");
    }
}

#[test]
fn derivative_partition_of_zero() {
    // The derivative of the constant 1 is 0, so Σ_i B_i^9'(t) ≡ 0.
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=9).map(|i| bernstein_binonic_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-10, "t = {t}, sum = {sum}");
    }
}

#[test]
fn second_derivative_partition_of_zero() {
    // The second derivative of 1 is 0, so Σ_i B_i^9''(t) ≡ 0.
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=9).map(|i| bernstein_binonic_second_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-9, "t = {t}, sum = {sum}");
    }
}

#[test]
fn endpoint_evaluation_matches_delta() {
    // B_i^9(0) = 1 iff i = 0, else 0.
    // B_i^9(1) = 1 iff i = 9, else 0.
    for i in 0..=9 {
        let v0 = bernstein_binonic(i, 0.0);
        let v1 = bernstein_binonic(i, 1.0);
        let expected_v0 = if i == 0 { 1.0 } else { 0.0 };
        let expected_v1 = if i == 9 { 1.0 } else { 0.0 };
        assert!(approx_eq(v0, expected_v0), "B_{i}^9(0) = {v0}");
        assert!(approx_eq(v1, expected_v1), "B_{i}^9(1) = {v1}");
    }
}

#[test]
fn specialised_agrees_with_generic_bernstein_basis() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        for i in 0..=9 {
            let specialised = bernstein_binonic(i, t);
            let generic = bernstein_basis(9, i, t);
            assert!(
                approx_eq(specialised, generic),
                "i = {i}, t = {t}: specialised = {specialised}, generic = {generic}"
            );
        }
    }
}

#[test]
fn derivative_matches_finite_difference() {
    // B_i^9'(t) ≈ (B_i^9(t + h) − B_i^9(t − h)) / (2h).
    let h = 1e-6;
    for k in 1..20 {
        let t = (k as f64) / 20.0;
        for i in 0..=9 {
            let analytic = bernstein_binonic_deriv(i, t);
            let fd = (bernstein_binonic(i, t + h) - bernstein_binonic(i, t - h)) / (2.0 * h);
            assert!(
                (analytic - fd).abs() < 1e-4,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn second_derivative_matches_finite_difference() {
    // B_i^9''(t) ≈ (B_i^9(t+h) − 2·B_i^9(t) + B_i^9(t-h)) / h².
    let h = 1e-4;
    for k in 2..18 {
        let t = (k as f64) / 20.0;
        for i in 0..=9 {
            let analytic = bernstein_binonic_second_deriv(i, t);
            let fd = (bernstein_binonic(i, t + h) - 2.0 * bernstein_binonic(i, t)
                + bernstein_binonic(i, t - h))
                / (h * h);
            assert!(
                (analytic - fd).abs() < 1e-2,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn bernstein_binonic_all_matches_individual() {
    // `bernstein_binonic_all(t)` must produce the same array as
    // 10 individual calls to `bernstein_binonic(i, t)`.
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let batched = bernstein_binonic_all(t);
        for i in 0..10 {
            let individual = bernstein_binonic(i, t);
            assert!(
                approx_eq(batched[i], individual),
                "i = {i}, t = {t}: batch = {}, individual = {}",
                batched[i],
                individual
            );
        }
    }
}

#[test]
fn bernstein_binonic_deriv_all_matches_individual() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let batched = bernstein_binonic_deriv_all(t);
        for i in 0..10 {
            let individual = bernstein_binonic_deriv(i, t);
            assert!(
                approx_eq(batched[i], individual),
                "i = {i}, t = {t}: batch = {}, individual = {}",
                batched[i],
                individual
            );
        }
    }
}

#[test]
fn bernstein_binonic_second_deriv_all_matches_individual() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let batched = bernstein_binonic_second_deriv_all(t);
        for i in 0..10 {
            let individual = bernstein_binonic_second_deriv(i, t);
            assert!(
                approx_eq(batched[i], individual),
                "i = {i}, t = {t}: batch = {}, individual = {}",
                batched[i],
                individual
            );
        }
    }
}

#[test]
fn out_of_range_returns_zero() {
    // i > 9 must return 0 for all three flavors.
    assert_eq!(bernstein_binonic(10, 0.5), 0.0);
    assert_eq!(bernstein_binonic(100, 0.5), 0.0);
    assert_eq!(bernstein_binonic_deriv(10, 0.5), 0.0);
    assert_eq!(bernstein_binonic_second_deriv(42, 0.5), 0.0);
}

#[test]
fn linear_polynomial_via_bernstein_binonic_sums_correctly() {
    // Express `f(t) = 2 + 3t` as a degree-9 Bernstein poly via
    // degree elevation from Bernstein (2, 5). Verify the
    // specialised basis sums correctly.
    let p = BernsteinPoly::<f64>::new(vec![2.0, 5.0]).degree_elevate(8);
    assert_eq!(p.degree(), 9);
    let coeffs = p.coeffs();
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let expected = 2.0 + 3.0 * t;
        let sum_via_specialised: f64 = (0..=9).map(|i| coeffs[i] * bernstein_binonic(i, t)).sum();
        assert!(
            (sum_via_specialised - expected).abs() < 1e-10,
            "t = {t}, specialised = {sum_via_specialised}, expected = {expected}"
        );
    }
}

#[test]
fn endpoint_derivatives_match_classical_formulas() {
    // At t = 0: B_i^9'(0) is nonzero only for i = 0 (= −9) and
    // i = 1 (= +9). All other i give 0.
    assert!(approx_eq(bernstein_binonic_deriv(0, 0.0), -9.0));
    assert!(approx_eq(bernstein_binonic_deriv(1, 0.0), 9.0));
    for i in 2..=9 {
        assert!(
            approx_eq(bernstein_binonic_deriv(i, 0.0), 0.0),
            "B_{i}^9'(0) = {}",
            bernstein_binonic_deriv(i, 0.0)
        );
    }
    // At t = 1: B_i^9'(1) is nonzero only for i = 8 (= −9) and
    // i = 9 (= +9).
    assert!(approx_eq(bernstein_binonic_deriv(8, 1.0), -9.0));
    assert!(approx_eq(bernstein_binonic_deriv(9, 1.0), 9.0));
    for i in 0..=7 {
        assert!(
            approx_eq(bernstein_binonic_deriv(i, 1.0), 0.0),
            "B_{i}^9'(1) = {}",
            bernstein_binonic_deriv(i, 1.0)
        );
    }
}

#[test]
fn endpoint_second_derivatives_match_classical_formulas() {
    // At t = 0: B_i^9''(0) is nonzero only for i ∈ {0, 1, 2}.
    //   B_0^9''(0) = 9·8 = 72
    //   B_1^9''(0) = −2 · 72 = −144
    //   B_2^9''(0) = 72
    assert!(approx_eq(bernstein_binonic_second_deriv(0, 0.0), 72.0));
    assert!(approx_eq(bernstein_binonic_second_deriv(1, 0.0), -144.0));
    assert!(approx_eq(bernstein_binonic_second_deriv(2, 0.0), 72.0));
    for i in 3..=9 {
        assert!(
            approx_eq(bernstein_binonic_second_deriv(i, 0.0), 0.0),
            "B_{i}^9''(0) = {}",
            bernstein_binonic_second_deriv(i, 0.0)
        );
    }
    // At t = 1: nonzero only for i ∈ {7, 8, 9}.
    //   B_7^9''(1) = 72, B_8^9''(1) = −144, B_9^9''(1) = 72.
    assert!(approx_eq(bernstein_binonic_second_deriv(7, 1.0), 72.0));
    assert!(approx_eq(bernstein_binonic_second_deriv(8, 1.0), -144.0));
    assert!(approx_eq(bernstein_binonic_second_deriv(9, 1.0), 72.0));
}

#[test]
fn symmetric_at_half_peak_index_4_and_5() {
    // B_i^9(0.5) is symmetric about i = 4.5:
    //   B_0^9(0.5) = B_9^9(0.5), B_1^9 = B_8^9, etc.
    let t = 0.5;
    let values: Vec<f64> = (0..=9).map(|i| bernstein_binonic(i, t)).collect();
    for i in 0..=4 {
        assert!(
            approx_eq(values[i], values[9 - i]),
            "asymmetry at i = {i}: {} vs {}",
            values[i],
            values[9 - i]
        );
    }
}
