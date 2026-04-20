//! Tests for the degree-5 Bernstein basis + derivatives.
//!
//! Cross-checks the specialised implementations against the generic
//! [`super::super::super::bernstein_poly::bernstein_basis`] and
//! against finite-difference references.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::bernstein_poly::bernstein_basis;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

#[test]
fn partition_of_unity() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=5).map(|i| bernstein_bi5(i, t)).sum();
        assert!((sum - 1.0).abs() < 1e-12, "t = {t}, sum = {sum}");
    }
}

#[test]
fn derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=5).map(|i| bernstein_bi5_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-11, "t = {t}, sum = {sum}");
    }
}

#[test]
fn second_derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=5).map(|i| bernstein_bi5_second_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-11, "t = {t}, sum = {sum}");
    }
}

#[test]
fn third_derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=5).map(|i| bernstein_bi5_third_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-10, "t = {t}, sum = {sum}");
    }
}

#[test]
fn endpoint_evaluation_matches_delta() {
    for i in 0..=5 {
        let v0 = bernstein_bi5(i, 0.0);
        let v1 = bernstein_bi5(i, 1.0);
        let expected_v0 = if i == 0 { 1.0 } else { 0.0 };
        let expected_v1 = if i == 5 { 1.0 } else { 0.0 };
        assert!(approx_eq(v0, expected_v0), "B_{i}^5(0) = {v0}");
        assert!(approx_eq(v1, expected_v1), "B_{i}^5(1) = {v1}");
    }
}

#[test]
fn specialised_agrees_with_generic_bernstein_basis() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        for i in 0..=5 {
            let specialised = bernstein_bi5(i, t);
            let generic = bernstein_basis(5, i, t);
            assert!(
                approx_eq(specialised, generic),
                "i = {i}, t = {t}: specialised = {specialised}, generic = {generic}"
            );
        }
    }
}

#[test]
fn derivative_matches_finite_difference() {
    let h = 1e-6;
    for k in 1..20 {
        let t = (k as f64) / 20.0;
        for i in 0..=5 {
            let analytic = bernstein_bi5_deriv(i, t);
            let fd = (bernstein_bi5(i, t + h) - bernstein_bi5(i, t - h)) / (2.0 * h);
            assert!(
                (analytic - fd).abs() < 1e-5,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn second_derivative_matches_finite_difference() {
    let h = 1e-4;
    for k in 2..18 {
        let t = (k as f64) / 20.0;
        for i in 0..=5 {
            let analytic = bernstein_bi5_second_deriv(i, t);
            let fd = (bernstein_bi5(i, t + h) - 2.0 * bernstein_bi5(i, t)
                + bernstein_bi5(i, t - h))
                / (h * h);
            assert!(
                (analytic - fd).abs() < 1e-3,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn third_derivative_matches_finite_difference() {
    // Central FD for 3rd derivative:
    //   f'''(t) ≈ (f(t+2h) − 2·f(t+h) + 2·f(t-h) − f(t-2h)) / (2·h³).
    let h = 1e-3;
    for k in 4..16 {
        let t = (k as f64) / 20.0;
        for i in 0..=5 {
            let analytic = bernstein_bi5_third_deriv(i, t);
            let fd = (bernstein_bi5(i, t + 2.0 * h) - 2.0 * bernstein_bi5(i, t + h)
                + 2.0 * bernstein_bi5(i, t - h)
                - bernstein_bi5(i, t - 2.0 * h))
                / (2.0 * h * h * h);
            assert!(
                (analytic - fd).abs() < 1e-2,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn all_matches_individual() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let batched = bernstein_bi5_all(t);
        let batched_d = bernstein_bi5_deriv_all(t);
        let batched_d2 = bernstein_bi5_second_deriv_all(t);
        let batched_d3 = bernstein_bi5_third_deriv_all(t);
        for i in 0..6 {
            assert!(approx_eq(batched[i], bernstein_bi5(i, t)));
            assert!(approx_eq(batched_d[i], bernstein_bi5_deriv(i, t)));
            assert!(approx_eq(batched_d2[i], bernstein_bi5_second_deriv(i, t)));
            assert!(approx_eq(batched_d3[i], bernstein_bi5_third_deriv(i, t)));
        }
    }
}

#[test]
fn out_of_range_returns_zero() {
    assert_eq!(bernstein_bi5(6, 0.5), 0.0);
    assert_eq!(bernstein_bi5(100, 0.5), 0.0);
}

#[test]
fn endpoint_derivative_classical() {
    // B_i^5'(0) is nonzero only for i = 0 (= −5) and i = 1 (= +5).
    assert!(approx_eq(bernstein_bi5_deriv(0, 0.0), -5.0));
    assert!(approx_eq(bernstein_bi5_deriv(1, 0.0), 5.0));
    for i in 2..=5 {
        assert!(approx_eq(bernstein_bi5_deriv(i, 0.0), 0.0));
    }
    // B_i^5'(1): nonzero only for i = 4 (= −5) and i = 5 (= +5).
    assert!(approx_eq(bernstein_bi5_deriv(4, 1.0), -5.0));
    assert!(approx_eq(bernstein_bi5_deriv(5, 1.0), 5.0));
    for i in 0..=3 {
        assert!(approx_eq(bernstein_bi5_deriv(i, 1.0), 0.0));
    }
}

#[test]
fn endpoint_second_derivative_classical() {
    // At t = 0, B_i^5''(0) nonzero only for i ∈ {0, 1, 2}.
    //   B_0^5''(0) = n·(n−1) = 20
    //   B_1^5''(0) = −2 · 20 = −40
    //   B_2^5''(0) = 20
    assert!(approx_eq(bernstein_bi5_second_deriv(0, 0.0), 20.0));
    assert!(approx_eq(bernstein_bi5_second_deriv(1, 0.0), -40.0));
    assert!(approx_eq(bernstein_bi5_second_deriv(2, 0.0), 20.0));
    // At t = 1, B_i^5''(1) nonzero only for i ∈ {3, 4, 5}.
    assert!(approx_eq(bernstein_bi5_second_deriv(3, 1.0), 20.0));
    assert!(approx_eq(bernstein_bi5_second_deriv(4, 1.0), -40.0));
    assert!(approx_eq(bernstein_bi5_second_deriv(5, 1.0), 20.0));
}

#[test]
fn endpoint_third_derivative_classical() {
    // B_i^5'''(0) nonzero for i ∈ {0, 1, 2, 3}.
    //   B_0^5'''(0) = −60
    //   B_1^5'''(0) =  180
    //   B_2^5'''(0) = −180
    //   B_3^5'''(0) =   60
    assert!(approx_eq(bernstein_bi5_third_deriv(0, 0.0), -60.0));
    assert!(approx_eq(bernstein_bi5_third_deriv(1, 0.0), 180.0));
    assert!(approx_eq(bernstein_bi5_third_deriv(2, 0.0), -180.0));
    assert!(approx_eq(bernstein_bi5_third_deriv(3, 0.0), 60.0));
    for i in 4..=5 {
        assert!(approx_eq(bernstein_bi5_third_deriv(i, 0.0), 0.0));
    }
}

#[test]
fn symmetric_at_half() {
    let t = 0.5;
    let values: Vec<f64> = (0..=5).map(|i| bernstein_bi5(i, t)).collect();
    for i in 0..=2 {
        assert!(approx_eq(values[i], values[5 - i]));
    }
}
