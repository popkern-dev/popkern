//! Tests for the degree-3 Bernstein basis + derivatives.
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
        let sum: f64 = (0..=3).map(|i| bernstein_bi3(i, t)).sum();
        assert!((sum - 1.0).abs() < 1e-12, "t = {t}, sum = {sum}");
    }
}

#[test]
fn derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=3).map(|i| bernstein_bi3_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-11, "t = {t}, sum = {sum}");
    }
}

#[test]
fn second_derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=3).map(|i| bernstein_bi3_second_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-11, "t = {t}, sum = {sum}");
    }
}

#[test]
fn third_derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=3).map(|i| bernstein_bi3_third_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-10, "t = {t}, sum = {sum}");
    }
}

#[test]
fn endpoint_evaluation_matches_delta() {
    for i in 0..=3 {
        let v0 = bernstein_bi3(i, 0.0);
        let v1 = bernstein_bi3(i, 1.0);
        let expected_v0 = if i == 0 { 1.0 } else { 0.0 };
        let expected_v1 = if i == 3 { 1.0 } else { 0.0 };
        assert!(approx_eq(v0, expected_v0), "B_{i}^3(0) = {v0}");
        assert!(approx_eq(v1, expected_v1), "B_{i}^3(1) = {v1}");
    }
}

#[test]
fn specialised_agrees_with_generic_bernstein_basis() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        for i in 0..=3 {
            let specialised = bernstein_bi3(i, t);
            let generic = bernstein_basis(3, i, t);
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
        for i in 0..=3 {
            let analytic = bernstein_bi3_deriv(i, t);
            let fd = (bernstein_bi3(i, t + h) - bernstein_bi3(i, t - h)) / (2.0 * h);
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
        for i in 0..=3 {
            let analytic = bernstein_bi3_second_deriv(i, t);
            let fd = (bernstein_bi3(i, t + h) - 2.0 * bernstein_bi3(i, t)
                + bernstein_bi3(i, t - h))
                / (h * h);
            assert!(
                (analytic - fd).abs() < 1e-3,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn third_derivative_is_piecewise_constant() {
    // B^3'''(t) is a constant on [0, 1] for each basis function.
    // Specifically B^3'''(t) = 6 · (δ_{i,0} − 3δ_{i,1} + 3δ_{i,2} − δ_{i,3})
    // when t ∈ (0, 1), which gives:
    //   B_0^3''' = 6,  B_1^3''' = −18,  B_2^3''' = 18,  B_3^3''' = −6
    //
    // BUT our formula uses B_j^0(t) which is 1 iff j==0, so at a
    // point t ∈ (0, 1):
    //   B_i^3'''(t) = 6·(B_{i-3}^0 − 3·B_{i-2}^0 + 3·B_{i-1}^0 − B_i^0)
    //              = 6·((i==3) − 3·(i==2) + 3·(i==1) − (i==0))
    //              = 6 · [−1, 3, −3, 1][i]
    //              = [−6, 18, −18, 6][i]
    //
    // Note the sign convention — B_i^0 is zero for i ≠ 0 so only
    // the non-trivial index contributes.
    let expected = [-6.0, 18.0, -18.0, 6.0];
    for k in 1..20 {
        let t = (k as f64) / 20.0;
        for i in 0..=3 {
            let analytic = bernstein_bi3_third_deriv(i, t);
            assert!(
                (analytic - expected[i]).abs() < 1e-10,
                "i = {i}, t = {t}: analytic = {analytic}, expected = {}",
                expected[i]
            );
        }
    }
}

#[test]
fn all_at_once_matches_individual_calls() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let all = bernstein_bi3_all(t);
        for i in 0..=3 {
            assert!(approx_eq(all[i], bernstein_bi3(i, t)), "i = {i}, t = {t}");
        }

        let dall = bernstein_bi3_deriv_all(t);
        for i in 0..=3 {
            assert!(
                (dall[i] - bernstein_bi3_deriv(i, t)).abs() < EPS,
                "deriv i = {i}, t = {t}"
            );
        }

        let ddall = bernstein_bi3_second_deriv_all(t);
        for i in 0..=3 {
            assert!(
                (ddall[i] - bernstein_bi3_second_deriv(i, t)).abs() < EPS,
                "second deriv i = {i}, t = {t}"
            );
        }
    }
}

#[test]
fn out_of_range_returns_zero() {
    for i in 4..10 {
        assert_eq!(bernstein_bi3(i, 0.3), 0.0);
    }
}
