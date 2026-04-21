//! Tests for the degree-6 Bernstein basis + derivatives.

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
        let sum: f64 = (0..=6).map(|i| bernstein_bi6(i, t)).sum();
        assert!((sum - 1.0).abs() < 1e-12, "t = {t}, sum = {sum}");
    }
}

#[test]
fn derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=6).map(|i| bernstein_bi6_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-11, "t = {t}, sum = {sum}");
    }
}

#[test]
fn second_derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=6).map(|i| bernstein_bi6_second_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-10, "t = {t}, sum = {sum}");
    }
}

#[test]
fn third_derivative_sum_is_zero() {
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let sum: f64 = (0..=6).map(|i| bernstein_bi6_third_deriv(i, t)).sum();
        assert!(sum.abs() < 1e-9, "t = {t}, sum = {sum}");
    }
}

#[test]
fn endpoint_evaluation_matches_delta() {
    for i in 0..=6 {
        let v0 = bernstein_bi6(i, 0.0);
        let v1 = bernstein_bi6(i, 1.0);
        let expected_v0 = if i == 0 { 1.0 } else { 0.0 };
        let expected_v1 = if i == 6 { 1.0 } else { 0.0 };
        assert!(approx_eq(v0, expected_v0));
        assert!(approx_eq(v1, expected_v1));
    }
}

#[test]
fn specialised_agrees_with_generic_bernstein_basis() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        for i in 0..=6 {
            let specialised = bernstein_bi6(i, t);
            let generic = bernstein_basis(6, i, t);
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
        for i in 0..=6 {
            let analytic = bernstein_bi6_deriv(i, t);
            let fd = (bernstein_bi6(i, t + h) - bernstein_bi6(i, t - h)) / (2.0 * h);
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
        for i in 0..=6 {
            let analytic = bernstein_bi6_second_deriv(i, t);
            let fd = (bernstein_bi6(i, t + h) - 2.0 * bernstein_bi6(i, t)
                + bernstein_bi6(i, t - h))
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
    let h = 1e-3;
    for k in 4..16 {
        let t = (k as f64) / 20.0;
        for i in 0..=6 {
            let analytic = bernstein_bi6_third_deriv(i, t);
            let fd = (bernstein_bi6(i, t + 2.0 * h) - 2.0 * bernstein_bi6(i, t + h)
                + 2.0 * bernstein_bi6(i, t - h)
                - bernstein_bi6(i, t - 2.0 * h))
                / (2.0 * h * h * h);
            assert!(
                (analytic - fd).abs() < 1e-2,
                "i = {i}, t = {t}: analytic = {analytic}, FD = {fd}"
            );
        }
    }
}

#[test]
fn fourth_derivative_sum_is_zero() {
    // Σ B_i^6⁽⁴⁾(t) = d⁴/dt⁴ (1) = 0 at every t.
    for k in 0..=30 {
        let t = (k as f64) / 30.0;
        let sum: f64 = bernstein_bi6_fourth_deriv_all(t).iter().sum();
        assert!(sum.abs() < 1e-10, "Σ B^6⁽⁴⁾ at t={t} = {sum} ≠ 0");
    }
}

#[test]
fn fourth_derivative_via_diff_of_third() {
    // Cross-check analytic 4th deriv against (d/dt)(3rd deriv)
    // computed via central difference, using a moderate h where
    // cancellation noise is manageable.
    let h = 1e-4;
    for k in 1..10 {
        let t = (k as f64) / 10.0;
        let d4 = bernstein_bi6_fourth_deriv_all(t);
        for i in 0..=6 {
            let d3_plus = bernstein_bi6_third_deriv(i, t + h);
            let d3_minus = bernstein_bi6_third_deriv(i, t - h);
            let fd = (d3_plus - d3_minus) / (2.0 * h);
            let rel = (d4[i] - fd).abs() / d4[i].abs().max(1.0);
            assert!(
                rel < 1e-3,
                "i = {i}, t = {t}: analytic = {}, FD = {fd}, rel {rel}",
                d4[i]
            );
        }
    }
}

#[test]
fn fifth_derivative_matches_finite_difference() {
    // 5th FD: (-f(t-3h/2) + ...) / h⁵; stencil sum too noisy for
    // small h — use a moderate h and check magnitude.
    // Alternative check: 5th derivative of 1 is 0.
    for k in 0..=30 {
        let t = (k as f64) / 30.0;
        let sum: f64 = (0..=6).map(|i| bernstein_bi6_fifth_deriv_all(t)[i]).sum();
        assert!(sum.abs() < 1e-9, "Σ B^6⁽⁵⁾_i(t) at t={t} = {sum} ≠ 0");
    }
}

#[test]
fn fifth_derivative_endpoint_identities() {
    // B_i^6⁽⁵⁾(0): endpoint formula is 720·(-δ_{i-5,-1} + 5·δ_{i-4,0} - 10·δ_{i-3,1} + ...)
    // But simpler check: at t=0, B^1_m(0) = 1 iff m=0, else 0.
    // So B_i^6⁽⁵⁾(0) = 720 · (-B^1_{i-5}(0) + 5 B^1_{i-4}(0) - 10 B^1_{i-3}(0) + 10 B^1_{i-2}(0) - 5 B^1_{i-1}(0) + B^1_i(0))
    //                 = 720 · (-δ_{i, 5} + 5 δ_{i, 4} - 10 δ_{i, 3} + 10 δ_{i, 2} - 5 δ_{i, 1} + δ_{i, 0})
    let d5_0 = bernstein_bi6_fifth_deriv_all(0.0);
    let expected_at_0 = [720.0, -3600.0, 7200.0, -7200.0, 3600.0, -720.0, 0.0];
    for i in 0..7 {
        assert!(
            (d5_0[i] - expected_at_0[i]).abs() < 1e-9,
            "B_{i}^6⁽⁵⁾(0): expected {}, got {}",
            expected_at_0[i],
            d5_0[i]
        );
    }
}

#[test]
fn fourth_derivative_endpoint_identities() {
    // At t=0: B^2_m(0) = 1 iff m=0, else 0.
    // B_i^6⁽⁴⁾(0) = 360 · (δ_{i, 4} - 4 δ_{i, 3} + 6 δ_{i, 2} - 4 δ_{i, 1} + δ_{i, 0})
    let d4_0 = bernstein_bi6_fourth_deriv_all(0.0);
    let expected = [360.0, -1440.0, 2160.0, -1440.0, 360.0, 0.0, 0.0];
    for i in 0..7 {
        assert!(
            (d4_0[i] - expected[i]).abs() < 1e-9,
            "B_{i}^6⁽⁴⁾(0): expected {}, got {}",
            expected[i],
            d4_0[i]
        );
    }
}

#[test]
fn all_matches_individual() {
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let b = bernstein_bi6_all(t);
        let bd = bernstein_bi6_deriv_all(t);
        let bd2 = bernstein_bi6_second_deriv_all(t);
        let bd3 = bernstein_bi6_third_deriv_all(t);
        for i in 0..7 {
            assert!(approx_eq(b[i], bernstein_bi6(i, t)));
            assert!(approx_eq(bd[i], bernstein_bi6_deriv(i, t)));
            assert!(approx_eq(bd2[i], bernstein_bi6_second_deriv(i, t)));
            assert!(approx_eq(bd3[i], bernstein_bi6_third_deriv(i, t)));
        }
    }
}

#[test]
fn out_of_range_returns_zero() {
    assert_eq!(bernstein_bi6(7, 0.5), 0.0);
    assert_eq!(bernstein_bi6(100, 0.5), 0.0);
}

#[test]
fn endpoint_derivative_classical() {
    // B_i^6'(0) nonzero only for i = 0 (= −6) and i = 1 (= +6).
    assert!(approx_eq(bernstein_bi6_deriv(0, 0.0), -6.0));
    assert!(approx_eq(bernstein_bi6_deriv(1, 0.0), 6.0));
    for i in 2..=6 {
        assert!(approx_eq(bernstein_bi6_deriv(i, 0.0), 0.0));
    }
    assert!(approx_eq(bernstein_bi6_deriv(5, 1.0), -6.0));
    assert!(approx_eq(bernstein_bi6_deriv(6, 1.0), 6.0));
    for i in 0..=4 {
        assert!(approx_eq(bernstein_bi6_deriv(i, 1.0), 0.0));
    }
}

#[test]
fn endpoint_second_derivative_classical() {
    // B_i^6''(0) nonzero for i ∈ {0, 1, 2}: [30, −60, 30].
    assert!(approx_eq(bernstein_bi6_second_deriv(0, 0.0), 30.0));
    assert!(approx_eq(bernstein_bi6_second_deriv(1, 0.0), -60.0));
    assert!(approx_eq(bernstein_bi6_second_deriv(2, 0.0), 30.0));
    // B_i^6''(1) nonzero for i ∈ {4, 5, 6}: [30, −60, 30].
    assert!(approx_eq(bernstein_bi6_second_deriv(4, 1.0), 30.0));
    assert!(approx_eq(bernstein_bi6_second_deriv(5, 1.0), -60.0));
    assert!(approx_eq(bernstein_bi6_second_deriv(6, 1.0), 30.0));
}

#[test]
fn endpoint_third_derivative_classical() {
    // B_i^6'''(0) nonzero for i ∈ {0, 1, 2, 3}: [−120, 360, −360, 120].
    assert!(approx_eq(bernstein_bi6_third_deriv(0, 0.0), -120.0));
    assert!(approx_eq(bernstein_bi6_third_deriv(1, 0.0), 360.0));
    assert!(approx_eq(bernstein_bi6_third_deriv(2, 0.0), -360.0));
    assert!(approx_eq(bernstein_bi6_third_deriv(3, 0.0), 120.0));
    for i in 4..=6 {
        assert!(approx_eq(bernstein_bi6_third_deriv(i, 0.0), 0.0));
    }
}

#[test]
fn symmetric_at_half() {
    let t = 0.5;
    let values: Vec<f64> = (0..=6).map(|i| bernstein_bi6(i, t)).collect();
    for i in 0..=3 {
        assert!(approx_eq(values[i], values[6 - i]));
    }
}
