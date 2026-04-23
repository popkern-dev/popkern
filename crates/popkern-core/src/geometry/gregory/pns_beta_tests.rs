//! Tests for [`BetaReparam`] against the paper's closed-form beta
//! for bi-3 (n ∈ {3, 5}) and the eq-(3) formulas for bi-4 (n ≥ 5).

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-12;

fn approx(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

#[test]
fn bi4_defaults_match_paper_formulas() {
    // Paper eq (3) + Appendix B:
    //   a_0 = 1,  a_1 = (2 − c)/2 · a_2,  a_2^def = (n − 2)/(0.55n − 0.2)
    //   b_0 = 0,  b_1 = −2/3 + (2 − c)/3 · a_2,  b_2^def = 0,  b_3 = 0
    for n in [5_usize, 6, 7, 8, 10] {
        let beta = BetaReparam::bi4(n);
        let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();
        let a_2_expected = (n as f64 - 2.0) / (0.55 * (n as f64) - 0.2);
        let a_1_expected = (2.0 - c) / 2.0 * a_2_expected;
        let b_1_expected = -2.0 / 3.0 + (2.0 - c) / 3.0 * a_2_expected;
        assert!(approx(beta.a[0], 1.0, EPS), "n={n} a_0");
        assert!(approx(beta.a[1], a_1_expected, EPS), "n={n} a_1");
        assert!(approx(beta.a[2], a_2_expected, EPS), "n={n} a_2");
        assert!(approx(beta.b[0], 0.0, EPS), "n={n} b_0");
        assert!(approx(beta.b[1], b_1_expected, EPS), "n={n} b_1");
        assert!(approx(beta.b[2], 0.0, EPS), "n={n} b_2");
        assert!(approx(beta.b[3], 0.0, EPS), "n={n} b_3");
    }
}

#[test]
fn bi3_n3_matches_paper_closed_form() {
    // For n=3: c = cos(2π/3) = −1/2, so (1 − c) = 3/2, inv = 2/3.
    // a(u) = (1 − c(1 − u))/(1 − c) = 1 − c/(1 − c) + (c/(1 − c))·u
    //      = 1 − (−1/2)/(3/2) + ((−1/2)/(3/2))·u
    //      = 1 + 1/3 − (1/3)·u
    //      = 4/3 − u/3
    // At u=0: a(0) = 4/3.  WAIT. The paper says a(0) should be 1 at the boundary
    // — let me re-check.
    //
    // Actually reading the paper again, for bi-3 (eq 3³):
    //   a(u) := (1 − c(1 − u))/(1 − c)
    //   b(u) := c(1 − u)u/(1 − c)
    // At u=1: a(1) = 1/(1 − c). At u=0: a(0) = (1 − c)/(1 − c) = 1.
    // So a(0) = 1 indeed. My expansion above was off — let me redo:
    //   a(u) = (1/(1−c)) − (c/(1−c))·(1 − u)
    //        = (1/(1−c)) − c/(1−c) + (c/(1−c))·u
    //        = (1 − c)/(1 − c) + (c/(1 − c))·u
    //        = 1 + (c/(1 − c))·u
    // For n=3, c = −1/2, so c/(1−c) = (−1/2)/(3/2) = −1/3.
    // Thus a(u) = 1 − u/3.
    // At u=0: a=1. At u=1: a = 2/3. At u=0.5: a = 5/6.
    let beta = BetaReparam::bi3_for_n3_n5(3);
    assert!(approx(beta.eval_a(0.0), 1.0, EPS), "a(0) for n=3");
    assert!(approx(beta.eval_a(1.0), 2.0 / 3.0, EPS), "a(1) for n=3");
    assert!(approx(beta.eval_a(0.5), 5.0 / 6.0, EPS), "a(0.5) for n=3");

    // b(u) = c(1 − u)u/(1 − c) = (−1/2)(1 − u)u/(3/2) = −(1 − u)u/3
    // At u=0: 0. At u=0.5: −0.25/3 = −1/12. At u=1: 0.
    assert!(approx(beta.eval_b(0.0), 0.0, EPS), "b(0) for n=3");
    assert!(approx(beta.eval_b(1.0), 0.0, EPS), "b(1) for n=3");
    assert!(approx(beta.eval_b(0.5), -1.0 / 12.0, EPS), "b(0.5) for n=3");
}

#[test]
fn bi3_n5_matches_paper_closed_form() {
    // For n=5: c = cos(2π/5) = (√5 − 1)/4 ≈ 0.309017, (1 − c) ≈ 0.690983.
    let beta = BetaReparam::bi3_for_n3_n5(5);
    let c = (2.0 * std::f64::consts::PI / 5.0).cos();
    let inv = 1.0 / (1.0 - c);
    for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let a_expected = (1.0 - c * (1.0 - u)) * inv;
        let b_expected = c * (1.0 - u) * u * inv;
        assert!(
            approx(beta.eval_a(u), a_expected, 1e-12),
            "n=5 a({u}): got {} expected {a_expected}",
            beta.eval_a(u)
        );
        assert!(
            approx(beta.eval_b(u), b_expected, 1e-12),
            "n=5 b({u}): got {} expected {b_expected}",
            beta.eval_b(u)
        );
    }
}

#[test]
fn eval_consistency() {
    // eval(u, v) = (u + b(u) v, a(u) v).
    let beta = BetaReparam::bi4(6);
    for &u in &[0.0, 0.3, 0.5, 0.7, 1.0] {
        for &v in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            let (bx, by) = beta.eval(u, v);
            assert!(approx(bx, u + beta.eval_b(u) * v, EPS));
            assert!(approx(by, beta.eval_a(u) * v, EPS));
        }
    }
}

#[test]
fn bi3_preserves_boundary_endpoints() {
    // At v=0, beta(u, 0) = (u, 0) — beta preserves the boundary of the
    // spline complex (Section 4.1).
    for n in [3_usize, 5] {
        let beta = BetaReparam::bi3_for_n3_n5(n);
        for &u in &[0.0, 0.2, 0.4, 0.6, 0.8, 1.0] {
            let (bx, by) = beta.eval(u, 0.0);
            assert!(approx(bx, u, EPS), "beta_x(u,0) ≠ u at n={n}, u={u}");
            assert!(approx(by, 0.0, EPS), "beta_y(u,0) ≠ 0 at n={n}, u={u}");
        }
    }
}

#[test]
fn bi4_preserves_boundary_endpoints() {
    for n in [5_usize, 6, 7, 8] {
        let beta = BetaReparam::bi4(n);
        for &u in &[0.0, 0.2, 0.4, 0.6, 0.8, 1.0] {
            let (bx, by) = beta.eval(u, 0.0);
            assert!(approx(bx, u, EPS), "beta_x(u,0) ≠ u at n={n}, u={u}");
            assert!(approx(by, 0.0, EPS), "beta_y(u,0) ≠ 0 at n={n}, u={u}");
        }
    }
}

#[test]
fn bi4_n4_identity_limit() {
    // Appendix B: "consistency with the case n = 4 when beta is the
    // identity". At n = 4, c = cos(π/2) = 0, so:
    //   a_2^def = (4 − 2)/(2.2 − 0.2) = 2/2 = 1
    //   a_1 = (2 − 0)/2 · 1 = 1
    //   a_0 = 1
    //   b_1 = −2/3 + (2 − 0)/3 · 1 = 0
    //   All b's = 0
    // So a(u) ≡ 1, b(u) ≡ 0 → beta(u, v) = (u, v), the identity. ✓
    let beta = BetaReparam::bi4(4);
    for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        for &v in &[0.0, 0.5, 1.0] {
            let (bx, by) = beta.eval(u, v);
            assert!(approx(bx, u, EPS), "beta_x should be u at n=4");
            assert!(approx(by, v, EPS), "beta_y should be v at n=4");
        }
    }
}
