//! Tests for the K-P 2015 bi-6 reparameterisations ρ and β.
//!
//! Strategy: verify endpoint identities and spot-check γ table
//! entries against the published formulas (math summary §2.4,
//! transcribed from [KP15] Appendix).

use super::*;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

// ---- τ table ----

#[test]
fn tau_table_matches_published() {
    assert_eq!(tau_for_valence(3), 4.0 / 3.0);
    assert_eq!(tau_for_valence(5), 0.87);
    assert_eq!(tau_for_valence(6), 0.85);
    assert_eq!(tau_for_valence(7), 0.84);
    assert_eq!(tau_for_valence(8), 0.83);
    assert_eq!(tau_for_valence(9), 0.82);
    assert_eq!(tau_for_valence(10), 0.81);
    assert_eq!(tau_for_valence(11), 0.80);
    assert_eq!(tau_for_valence(20), 0.80);
    assert_eq!(tau_for_valence(100), 0.80);
}

#[test]
#[should_panic]
fn tau_for_valence_below_three_panics() {
    let _ = tau_for_valence(2);
}

// ---- Bi6Reparams constructor ----

#[test]
fn params_for_n5() {
    let p = Bi6Reparams::for_valence(5);
    assert_eq!(p.n, 5);
    assert!(approx_eq(p.c, (2.0 * std::f64::consts::PI / 5.0).cos()));
    assert_eq!(p.tau, 0.87);
    // d₂ = (1 − (3 + 2c)τ + (3 + c)τ²) / (2τ²)
    let c = p.c;
    let tau = p.tau;
    let expected_d2 = (1.0 - (3.0 + 2.0 * c) * tau + (3.0 + c) * tau * tau) / (2.0 * tau * tau);
    assert!(approx_eq(p.d2, expected_d2));
}

#[test]
fn params_for_n3_uses_0_point_1_default() {
    let p = Bi6Reparams::for_valence(3);
    assert_eq!(p.d2, 0.1);
    assert_eq!(p.tau, 4.0 / 3.0);
}

// ---- ρ (between-sector) ----

#[test]
fn rho_a_is_minus_one() {
    let p = Bi6Reparams::for_valence(5);
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        assert_eq!(p.rho_a(u), -1.0);
    }
}

#[test]
fn rho_b_endpoints() {
    let p = Bi6Reparams::for_valence(6);
    // b(0) = 2c
    assert!(approx_eq(p.rho_b(0.0), 2.0 * p.c));
    // b(1) = 0
    assert!(approx_eq(p.rho_b(1.0), 0.0));
}

#[test]
fn rho_b_prime_matches_fd() {
    let p = Bi6Reparams::for_valence(6);
    let h = 1e-6;
    for k in 1..10 {
        let u = k as f64 / 10.0;
        let analytic = p.rho_b_prime(u);
        let fd = (p.rho_b(u + h) - p.rho_b(u - h)) / (2.0 * h);
        assert!((analytic - fd).abs() < 1e-5);
    }
}

#[test]
fn rho_b_prime_endpoints() {
    // b'(u) = −4c(1 − u), so b'(0) = −4c, b'(1) = 0.
    let p = Bi6Reparams::for_valence(6);
    assert!(approx_eq(p.rho_b_prime(0.0), -4.0 * p.c));
    assert!(approx_eq(p.rho_b_prime(1.0), 0.0));
}

#[test]
fn rho_d_endpoints() {
    let p = Bi6Reparams::for_valence(6);
    // d(0) = δ · 1 / 1 = δ = 24c(τ − 1).
    let expected_delta = 24.0 * p.c * (p.tau - 1.0);
    assert!(approx_eq(p.rho_d(0.0), expected_delta));
    // d(1) = 0 (numerator vanishes).
    assert!(approx_eq(p.rho_d(1.0), 0.0));
}

#[test]
fn rho_d_is_zero_everywhere_when_tau_is_one() {
    // If τ = 1, δ = 0, so d(u) = 0 ∀ u.
    let mut p = Bi6Reparams::for_valence(6);
    p.tau = 1.0;
    // d₂ needs recomputing for consistency, but we're only testing
    // rho_d which depends on tau, not d2.
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        assert!(approx_eq(p.rho_d(u), 0.0));
    }
}

#[test]
fn rho_e_is_b_times_b_prime_minus_d_over_2() {
    let p = Bi6Reparams::for_valence(8);
    for k in 1..10 {
        let u = k as f64 / 10.0;
        let expected = p.rho_b(u) * (p.rho_b_prime(u) - p.rho_d(u) / 2.0);
        assert!(approx_eq(p.rho_e(u), expected));
    }
}

#[test]
fn rho_e_endpoints() {
    let p = Bi6Reparams::for_valence(6);
    // At u = 1: b(1) = 0 so e(1) = 0.
    assert!(approx_eq(p.rho_e(1.0), 0.0));
    // At u = 0: b(0) = 2c, b'(0) = −4c, d(0) = δ = 24c(τ-1).
    //   e(0) = 2c · (−4c − 12c(τ−1))
    //        = 2c · (−4c − 12cτ + 12c)
    //        = 2c · (8c − 12cτ)
    //        = 16c² − 24c²τ
    //        = 8c² · (2 − 3τ)
    let expected = 8.0 * p.c * p.c * (2.0 - 3.0 * p.tau);
    assert!(approx_eq(p.rho_e(0.0), expected));
}

// ---- β (tensor-border) ----

#[test]
fn beta_d_bb_coeffs_shape() {
    let p = Bi6Reparams::for_valence(6);
    let [d0, d1, d2, d3] = p.beta_d_bb_coeffs();
    assert_eq!(d0, 0.0);
    assert_eq!(d1, 0.0);
    assert!(approx_eq(d2, p.d2));
    assert!(approx_eq(d3 - d2, 2.0 * p.c / 3.0));
}

#[test]
fn beta_e_bb_coeffs_shape() {
    let p = Bi6Reparams::for_valence(6);
    let [e0, e1, e2, e3, e4] = p.beta_e_bb_coeffs();
    assert_eq!(e0, 0.0);
    assert_eq!(e1, 0.0);
    assert_eq!(e4, 0.0);
    assert!(approx_eq(e2, p.d2 / 2.0));
    let expected_e3 = (4.0 * p.c * (1.0 - p.tau) + p.tau * p.d2) / (4.0 * p.tau);
    assert!(approx_eq(e3, expected_e3));
}

#[test]
fn beta_bb_coeffs_for_n3_default() {
    let p = Bi6Reparams::for_valence(3);
    let [_, _, d2, d3] = p.beta_d_bb_coeffs();
    assert_eq!(d2, 0.1);
    assert!(approx_eq(d3, 2.0 * p.c / 3.0 + 0.1));
}

// ---- γ table ----

#[test]
fn gamma_k2_spot_check() {
    let p = Bi6Reparams::for_valence(6);
    assert!(approx_eq(gamma(&p, 2, 0, 0), 1.0 - p.d2 / 25.0));
    assert!(approx_eq(gamma(&p, 2, 1, 0), (6.0 + p.d2) / 50.0));
    assert!(approx_eq(gamma(&p, 2, 2, 0), 1.0 / 25.0));
    assert!(approx_eq(gamma(&p, 2, 1, 1), 9.0 / 25.0));
    assert!(approx_eq(gamma(&p, 2, 2, 1), 3.0 / 25.0));
    assert!(approx_eq(gamma(&p, 2, 2, 2), 1.0 / 25.0));
}

#[test]
fn gamma_k3_spot_check() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let tau = p.tau;
    let d2 = p.d2;
    assert!(approx_eq(
        gamma(&p, 3, 0, 0),
        ((3.0 + 5.0 * c - 3.0 * d2) * tau - 6.0 * c) / (300.0 * tau)
    ));
    assert!(approx_eq(gamma(&p, 3, 3, 2), 1.0 / 100.0));
    assert!(approx_eq(gamma(&p, 3, 0, 2), 1.0 / 100.0));
    assert!(approx_eq(gamma(&p, 3, 2, 1), 27.0 / 100.0));
}

#[test]
fn gamma_k6_spot_check() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let d2 = p.d2;
    assert!(approx_eq(
        gamma(&p, 6, 3, 0),
        (6.0 - 2.0 * c - 3.0 * d2) / 30.0
    ));
    assert!(approx_eq(
        gamma(&p, 6, 3, 1),
        (18.0 + 2.0 * c + 3.0 * d2) / 30.0
    ));
    assert!(approx_eq(gamma(&p, 6, 3, 2), 1.0 / 5.0));
}

#[test]
fn gamma_outside_sparsity_returns_zero() {
    let p = Bi6Reparams::for_valence(6);
    // Positions never exercised in the γ table.
    assert_eq!(gamma(&p, 2, 3, 0), 0.0);
    assert_eq!(gamma(&p, 2, 0, 1), 0.0);
    assert_eq!(gamma(&p, 5, 0, 0), 0.0);
    assert_eq!(gamma(&p, 5, 1, 0), 0.0);
    assert_eq!(gamma(&p, 6, 0, 0), 0.0);
    assert_eq!(gamma(&p, 6, 2, 0), 0.0);
    // k outside 2..=6.
    assert_eq!(gamma(&p, 0, 0, 0), 0.0);
    assert_eq!(gamma(&p, 1, 0, 0), 0.0);
    assert_eq!(gamma(&p, 7, 0, 0), 0.0);
    // j outside 0..=2.
    assert_eq!(gamma(&p, 3, 0, 3), 0.0);
}

#[test]
fn gamma_is_linear_in_parameters() {
    // γ is a rational function of (c, τ, d2). With two param
    // sets that agree on (c, τ) but differ in d2, the γ values
    // that depend on d2 should differ linearly in d2.
    let mut p1 = Bi6Reparams::for_valence(6);
    let mut p2 = p1;
    p1.d2 = 0.0;
    p2.d2 = 1.0;

    // γ²_{00} = 1 − d2 / 25: Δ = −(1 − 0)/25 = −1/25.
    let diff = gamma(&p2, 2, 0, 0) - gamma(&p1, 2, 0, 0);
    assert!(approx_eq(diff, -1.0 / 25.0));

    // γ²_{10} = (6 + d2)/50: Δ = 1/50.
    let diff = gamma(&p2, 2, 1, 0) - gamma(&p1, 2, 1, 0);
    assert!(approx_eq(diff, 1.0 / 50.0));

    // γ²_{11} = 9/25 (constant): Δ = 0.
    assert!(approx_eq(gamma(&p2, 2, 1, 1) - gamma(&p1, 2, 1, 1), 0.0));
}

// ---- Cross-valence smoke tests ----

#[test]
fn all_valence_params_finite() {
    for n in [3usize, 4, 5, 6, 7, 8, 9, 10, 11, 12, 20] {
        let p = Bi6Reparams::for_valence(n);
        assert!(p.c.is_finite());
        assert!(p.tau.is_finite() && p.tau > 0.0);
        assert!(p.d2.is_finite());
        // ρ at interior samples is finite.
        for k in 1..10 {
            let u = k as f64 / 10.0;
            assert!(p.rho_b(u).is_finite());
            assert!(p.rho_d(u).is_finite());
            assert!(p.rho_e(u).is_finite());
        }
    }
}

#[test]
fn gamma_all_entries_finite_for_standard_valences() {
    for n in [3usize, 5, 6, 7, 8, 9, 10, 11, 12] {
        let p = Bi6Reparams::for_valence(n);
        for k in 0..=7 {
            for i in 0..=4 {
                for j in 0..=3 {
                    assert!(
                        gamma(&p, k, i, j).is_finite(),
                        "n={n} γ^{k}_{{{i},{j}}} not finite"
                    );
                }
            }
        }
    }
}
