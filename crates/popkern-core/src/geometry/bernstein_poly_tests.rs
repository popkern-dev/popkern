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
fn to_monomial_roundtrip_matches_evaluate_scalar() {
    // For random-ish scalar polynomials, Σ a_k · t^k must match
    // the BB-basis evaluation at every t.
    let polys = [
        vec![1.0, 2.0, 3.0],
        vec![-1.5, 0.0, 2.7, -0.4],
        vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        vec![1.0, -2.0, 3.0, 1.5, -0.7, 0.2, 4.0, -1.1],
    ];
    for coeffs in polys {
        let p = BernsteinPoly::new(coeffs);
        let mono = p.to_monomial_coefficients();
        for k in 0..=20 {
            let t = (k as f64) / 20.0;
            let mut expected = 0.0;
            let mut tk = 1.0;
            for &a in &mono {
                expected += a * tk;
                tk *= t;
            }
            assert!(
                approx_eq(expected, p.evaluate(t)),
                "deg={}, t={t}, mono-eval={expected}, bb-eval={}",
                p.degree(),
                p.evaluate(t)
            );
        }
    }
}

#[test]
fn to_monomial_endpoints_match_bb_endpoints() {
    // Endpoint identities: a_0 = c_0 (= p(0)) and Σ a_k = c_n (= p(1)).
    let p = BernsteinPoly::new(vec![1.0, 2.0, -3.0, 0.5, 4.1]);
    let mono = p.to_monomial_coefficients();
    assert!(approx_eq(mono[0], p.coeff(0)));
    let sum: f64 = mono.iter().sum();
    assert!(approx_eq(sum, p.coeff(p.degree())));
}

#[test]
fn to_monomial_vec3_roundtrip_matches_evaluate() {
    // Vector coefficients route through the same conversion via
    // the `PolyCoeff` trait. Useful for tangent-vector polynomials
    // in the K-P G¹ identity.
    let p = BernsteinPoly::new(vec![
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, -1.0, 0.5),
    ]);
    let mono = p.to_monomial_coefficients();
    for k in 0..=20 {
        let t = (k as f64) / 20.0;
        let mut expected = Vec3::zero();
        let mut tk = 1.0;
        for &a in &mono {
            expected = expected.add(a.scale(tk));
            tk *= t;
        }
        assert!(
            approx_eq_vec3(expected, p.evaluate(t)),
            "t={t}, mono={expected:?}, bb={:?}",
            p.evaluate(t)
        );
    }
}

#[test]
fn kp15_eq2_identity_zeros_at_u0_with_rule_s1_p01() {
    // Phase A cross-derivation [KP15 §3.4, eq (2)] with the
    // paper-prescribed choice `a := -1, b := 2c(1-u)²` (lines 116
    // and 152 of `docs/papers/15bi6g2.txt`). The identity then
    // reads, at the shared boundary v = 0:
    //
    //     ∂_v f̃(u, 0) + ∂_v f(u, 0) − 2c(1-u)² · ∂_u f(u, 0) ≡ 0.
    //
    // We construct this polynomial in `u` from BB-coefficient
    // primitives — differentiate, scale, multiply, add — using
    // unit basis values for the previous sector's p̀_{ij}, plug
    // rule_s1's ṕ_{01} formula into the current sector, and
    // confirm the resulting polynomial is zero at u = 0. This is
    // a positive control on the Phase A bernstein_poly toolkit:
    // any sign error in `differentiate`, `degree_elevate`,
    // `bernstein_product_scalar`, or `to_monomial_coefficients`
    // would surface here as a nonzero constant term.
    //
    // Three independent unit-basis loadings exercise the three
    // weights in `S1.p01`:
    //   (p̀_{00}=1, others 0)  ⇒ ṕ_{01} = 2(1−c)·1 = 2(1−c)
    //   (p̀_{10}=1, others 0)  ⇒ ṕ_{01} = 2c·1 = 2c
    //   (p̀_{01}=1, others 0)  ⇒ ṕ_{01} = (-1)·1 = -1
    let n = 6usize;
    let c = (2.0 * std::f64::consts::PI / n as f64).cos();
    let cases = [
        ("p̀_{00}=1", [1.0, 0.0, 0.0], 2.0 * (1.0 - c)),
        ("p̀_{10}=1", [0.0, 1.0, 0.0], 2.0 * c),
        ("p̀_{01}=1", [0.0, 0.0, 1.0], -1.0),
    ];
    for (label, [p00, p10, p01], expected_p01) in cases {
        // f(u, 0) row at v=0 for the previous sector. Loading just
        // the first two entries is enough: ∂_u f(u, 0) sees these
        // first, and the higher-degree BB coefs of f(u, 0) we set
        // to zero. The polynomial degree of the identity is then
        // bounded by the loaded support.
        let row_v0 = vec![p00, p10, 0.0, 0.0, 0.0, 0.0, 0.0];
        // f(u, 1) row at v=1: only p̀_{01} is loaded.
        let row_v1 = vec![p01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        // Strip dot-derivative: ∂_v f(u, 0) = 6·(f(u,1) − f(u,0))
        // at u along the BB row. As Bernstein polynomial in u:
        let strip_prev: Vec<f64> = (0..7).map(|i| 6.0 * (row_v1[i] - row_v0[i])).collect();
        let dv_prev = BernsteinPoly::new(strip_prev);

        // Current sector f̃: shared boundary ṕ_{i,0} = p̀_{i,0},
        // and ṕ_{0,1} from rule_s1; remaining ṕ_{i,1} zero.
        let p01_curr = 2.0 * (1.0 - c) * p00 + 2.0 * c * p10 + (-1.0) * p01;
        let row_v0_curr = row_v0.clone();
        let row_v1_curr = {
            let mut r = vec![0.0; 7];
            r[0] = p01_curr;
            r
        };
        let strip_curr: Vec<f64> = (0..7)
            .map(|i| 6.0 * (row_v1_curr[i] - row_v0_curr[i]))
            .collect();
        let dv_curr = BernsteinPoly::new(strip_curr);

        // ∂_u f(u, 0): degree-5 polynomial in u, BB coefs
        // 6·(p̀_{i+1,0} − p̀_{i,0}).
        let du_prev_coeffs: Vec<f64> = (0..6).map(|i| 6.0 * (row_v0[i + 1] - row_v0[i])).collect();
        let du_prev = BernsteinPoly::new(du_prev_coeffs);

        // b(u) = 2c·(1-u)² as Bernstein degree-2: [1, 0, 0] · 2c.
        let b = BernsteinPoly::new(vec![2.0 * c, 0.0, 0.0]);

        // b(u) · ∂_u f(u, 0): degree 7.
        let b_times_du = bernstein_product_scalar(&b, &du_prev);

        // Identity polynomial: ∂_v f̃ + ∂_v f − b · ∂_u f.
        // All three need to share a common degree before summing
        // BB-coefficients; degree-elevate dv_curr and dv_prev from
        // 6 → 7 to match b_times_du.
        let identity = dv_curr
            .clone()
            .degree_elevate(1)
            .add_elevated(&dv_prev.clone().degree_elevate(1))
            .add_elevated(&BernsteinPoly::new(
                b_times_du.coeffs().iter().map(|&c| -c).collect(),
            ));
        assert_eq!(identity.degree(), 7);

        // u^0 monomial coefficient = identity.evaluate(0) = c_0.
        let mono = identity.to_monomial_coefficients();
        assert!(
            mono[0].abs() < 1e-12,
            "{label}: rule_s1 ṕ_{{0,1}} = {p01_curr} (expected {expected_p01}); \
             identity u^0 coef = {} (should be 0)",
            mono[0]
        );
        assert!(
            approx_eq(p01_curr, expected_p01),
            "{label}: rule_s1 p01 = {p01_curr}, expected {expected_p01}"
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
