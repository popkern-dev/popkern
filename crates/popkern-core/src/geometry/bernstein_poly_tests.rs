//! Tests for [`super::BernsteinPoly`] and its operations.
//!
//! Philosophy: prefer point-value equivalence tests over direct
//! coefficient comparison. Bernstein arithmetic should match
//! power-basis evaluation at arbitrary `t`; if it does, the
//! coefficient formulas are correct.

use super::*;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

fn approx_eq_vec3(a: Vec3, b: Vec3) -> bool {
    approx_eq(a.x, b.x) && approx_eq(a.y, b.y) && approx_eq(a.z, b.z)
}

#[test]
fn basis_partition_of_unity() {
    // Σ_i B_i^n(t) ≡ 1 on [0, 1] for every n. Cornerstone identity.
    for n in 0..=12 {
        for k in 0..=32 {
            let t = (k as f64) / 32.0;
            let sum: f64 = (0..=n).map(|i| bernstein_basis(n, i, t)).sum();
            assert!((sum - 1.0).abs() < 1e-10, "n={n}, t={t}, sum={sum}");
        }
    }
}

#[test]
fn binomial_matches_pascal() {
    // First 10 rows of Pascal's triangle.
    let pascal: Vec<Vec<f64>> = (0..=9)
        .map(|n| (0..=n).map(|k| binomial(n, k)).collect())
        .collect();
    assert_eq!(pascal[0], vec![1.0]);
    assert_eq!(pascal[1], vec![1.0, 1.0]);
    assert_eq!(pascal[2], vec![1.0, 2.0, 1.0]);
    assert_eq!(pascal[3], vec![1.0, 3.0, 3.0, 1.0]);
    assert_eq!(pascal[7], vec![1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0]);
}

#[test]
fn eval_linear_endpoints_and_midpoint() {
    // p(t) = 2·(1 − t) + 5·t = 2·B_0^1 + 5·B_1^1.
    let p = BernsteinPoly::new(vec![2.0, 5.0]);
    assert!(approx_eq(p.evaluate(0.0), 2.0));
    assert!(approx_eq(p.evaluate(1.0), 5.0));
    assert!(approx_eq(p.evaluate(0.5), 3.5));
}

#[test]
fn eval_matches_de_casteljau() {
    // For a handful of random-ish polynomials and t values, direct
    // Bernstein-sum evaluation and de Casteljau agree.
    let polys = [
        vec![1.0, 2.0, 3.0],
        vec![-1.5, 0.0, 2.7, -0.4],
        vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    ];
    for coeffs in polys {
        let p = BernsteinPoly::new(coeffs);
        for k in 0..=20 {
            let t = (k as f64) / 20.0;
            assert!(approx_eq(p.evaluate(t), p.de_casteljau(t)));
        }
    }
}

#[test]
fn differentiate_linear_matches_slope() {
    // p(t) = 2 + 3t as Bernstein-linear (c_0, c_1) = (2, 5).
    // p'(t) ≡ 3.
    let p = BernsteinPoly::new(vec![2.0, 5.0]);
    let dp = p.differentiate();
    assert_eq!(dp.degree(), 0);
    for k in 0..=5 {
        let t = (k as f64) / 5.0;
        assert!(approx_eq(dp.evaluate(t), 3.0));
    }
}

#[test]
fn differentiate_quadratic_matches_analytic() {
    // p(t) = 1 + 4t − 3t² = B_0^2·1 + B_1^2·3 + B_2^2·2.
    // Verification: p(0) = 1, p(1) = 2, p(½) = 1 + 2 − 0.75 = 2.25 ≠ ?
    // Let me expand: (1 − t)² + 6t(1 − t) + 2t² = 1 − 2t + t² + 6t − 6t² + 2t²
    //                                          = 1 + 4t − 3t². ✓
    // p'(t) = 4 − 6t.
    let p = BernsteinPoly::new(vec![1.0, 3.0, 2.0]);
    let dp = p.differentiate();
    assert_eq!(dp.degree(), 1);
    assert!(approx_eq(dp.evaluate(0.0), 4.0));
    assert!(approx_eq(dp.evaluate(1.0), -2.0));
    assert!(approx_eq(dp.evaluate(0.5), 1.0));
}

#[test]
fn degree_elevate_preserves_point_values() {
    // Random cubic. After elevating by k, point values must match.
    let p = BernsteinPoly::new(vec![1.0, 2.5, -1.0, 3.0]);
    for by in 1..=5 {
        let q = p.clone().degree_elevate(by);
        assert_eq!(q.degree(), p.degree() + by);
        for k in 0..=20 {
            let t = (k as f64) / 20.0;
            assert!(
                approx_eq(p.evaluate(t), q.evaluate(t)),
                "by={by}, t={t}, p={}, q={}",
                p.evaluate(t),
                q.evaluate(t)
            );
        }
    }
}

#[test]
fn product_of_basis_polys_hits_expected_degree_and_values() {
    // (1 − t) · t has classical Bernstein form B_1^2 / 2.
    let a = BernsteinPoly::new(vec![1.0, 0.0]); // B_0^1 = 1 − t
    let b = BernsteinPoly::new(vec![0.0, 1.0]); // B_1^1 = t
    let c = bernstein_product_scalar(&a, &b);
    assert_eq!(c.degree(), 2);
    assert!(approx_eq(c.evaluate(0.0), 0.0));
    assert!(approx_eq(c.evaluate(1.0), 0.0));
    assert!(approx_eq(c.evaluate(0.5), 0.25));
}

#[test]
fn product_scalar_scalar_matches_pointwise() {
    // For random-ish polynomials, (a · b)(t) = a(t) · b(t) at
    // every sampled t.
    let a = BernsteinPoly::new(vec![1.0, -2.0, 3.0]);
    let b = BernsteinPoly::new(vec![0.5, 1.5, -1.0, 2.0]);
    let c = bernstein_product_scalar(&a, &b);
    assert_eq!(c.degree(), 5);
    for k in 0..=32 {
        let t = (k as f64) / 32.0;
        let expected = a.evaluate(t) * b.evaluate(t);
        assert!(approx_eq(c.evaluate(t), expected));
    }
}

#[test]
fn product_scalar_vec3_matches_pointwise() {
    // β(t) · T(t) in the cross-tangent blending identity.
    let a = BernsteinPoly::new(vec![1.0, 2.0]); // linear scalar
    let b = BernsteinPoly::new(vec![
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]);
    let c = bernstein_product_scalar_vec3(&a, &b);
    assert_eq!(c.degree(), 3);
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let expected = b.evaluate(t).scale(a.evaluate(t));
        assert!(approx_eq_vec3(c.evaluate(t), expected));
    }
}

#[test]
fn add_elevated_matches_pointwise_across_degrees() {
    // (a + b)(t) = a(t) + b(t) when one gets degree-elevated.
    let a = BernsteinPoly::new(vec![1.0, 2.0]); // degree 1
    let b = BernsteinPoly::new(vec![3.0, 1.0, -1.0, 2.0]); // degree 3
    let c = a.add_elevated(&b);
    assert_eq!(c.degree(), 3);
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        assert!(approx_eq(c.evaluate(t), a.evaluate(t) + b.evaluate(t)));
    }
}

#[test]
fn differentiate_then_evaluate_matches_finite_difference_in_limit() {
    // Sanity check on derivative implementation against a central
    // finite difference for a non-trivial polynomial.
    let p = BernsteinPoly::new(vec![0.5, -1.3, 2.1, 0.7, -0.9, 1.5, -0.2]); // degree 6
    let dp = p.differentiate();
    let h = 1e-6;
    for k in 1..=10 {
        let t = (k as f64) / 12.0;
        let fd = (p.evaluate(t + h) - p.evaluate(t - h)) / (2.0 * h);
        assert!(
            (dp.evaluate(t) - fd).abs() < 1e-6,
            "t={t}, dp={}, fd={fd}",
            dp.evaluate(t)
        );
    }
}

#[test]
fn product_vector_scalar_evaluated_equals_pointwise() {
    // Cross-check the bernstein_product_scalar_vec3 path at
    // many degrees.
    for m in 1..=5 {
        for n in 1..=5 {
            let a_coeffs: Vec<f64> = (0..=m).map(|i| 1.0 + (i as f64) * 0.5).collect();
            let a = BernsteinPoly::new(a_coeffs);
            let b_coeffs: Vec<Vec3> = (0..=n)
                .map(|j| Vec3::new(j as f64, -(j as f64), (j as f64).powi(2)))
                .collect();
            let b = BernsteinPoly::new(b_coeffs);
            let c = bernstein_product_scalar_vec3(&a, &b);
            assert_eq!(c.degree(), m + n);
            for k in 0..=8 {
                let t = (k as f64) / 8.0;
                let expected = b.evaluate(t).scale(a.evaluate(t));
                assert!(
                    approx_eq_vec3(c.evaluate(t), expected),
                    "m={m}, n={n}, t={t}"
                );
            }
        }
    }
}
