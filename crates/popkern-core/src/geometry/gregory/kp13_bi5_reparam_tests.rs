//! Tests for the KP13 bi-5 reparameterisations ρ and β.
//!
//! Strategy: verify endpoint identities and check formula-derived
//! values against the published constants (KP13 §4).

use super::*;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

// ---- constructor ----

#[test]
fn params_for_n5() {
    let p = Bi5Reparams::for_valence(5);
    assert_eq!(p.n, 5);
    assert!(approx_eq(p.c, (2.0 * std::f64::consts::PI / 5.0).cos()));
    p.validate().unwrap();
}

#[test]
fn params_for_n6() {
    let p = Bi5Reparams::for_valence(6);
    assert_eq!(p.n, 6);
    // c = cos(π/3) = 0.5
    assert!(approx_eq(p.c, 0.5));
}

#[test]
#[should_panic]
fn below_three_panics() {
    let _ = Bi5Reparams::for_valence(2);
}

// ---- ρ (between-sector) ----

#[test]
fn rho_a_is_minus_one() {
    let p = Bi5Reparams::for_valence(5);
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        assert_eq!(p.rho_a(u), -1.0);
    }
}

#[test]
fn rho_o_b_endpoints() {
    let p = Bi5Reparams::for_valence(6);
    // b(0) = 2c
    assert!(approx_eq(p.rho_o_b(0.0), 2.0 * p.c));
    // b(1) = 2c/3
    assert!(approx_eq(p.rho_o_b(1.0), 2.0 * p.c / 3.0));
}

#[test]
fn rho_o_b_prime_is_constant() {
    let p = Bi5Reparams::for_valence(6);
    let expected = -4.0 * p.c / 3.0;
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        assert!(approx_eq(p.rho_o_b_prime(u), expected));
    }
}

#[test]
fn rho_o_b_prime_matches_fd() {
    let p = Bi5Reparams::for_valence(7);
    let h = 1e-6;
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        let fd = (p.rho_o_b(u + h) - p.rho_o_b(u - h)) / (2.0 * h);
        assert!(
            (fd - p.rho_o_b_prime(u)).abs() < 1e-6,
            "FD mismatch at u = {u}"
        );
    }
}

#[test]
fn rho_o_d_is_zero() {
    let p = Bi5Reparams::for_valence(5);
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        assert_eq!(p.rho_o_d(u), 0.0);
    }
}

#[test]
fn rho_o_e_equals_b_times_bprime() {
    let p = Bi5Reparams::for_valence(6);
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        let b = p.rho_o_b(u);
        let bp = p.rho_o_b_prime(u);
        assert!(approx_eq(p.rho_o_e(u), b * bp));
    }
}

#[test]
fn rho_lr_b_endpoints() {
    let p = Bi5Reparams::for_valence(6);
    // b(0) = 2c/3
    assert!(approx_eq(p.rho_lr_b(0.0), 2.0 * p.c / 3.0));
    // b(1) = 0
    assert!(approx_eq(p.rho_lr_b(1.0), 0.0));
}

#[test]
fn rho_lr_b_prime_matches_fd() {
    let p = Bi5Reparams::for_valence(7);
    let h = 1e-6;
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        let fd = (p.rho_lr_b(u + h) - p.rho_lr_b(u - h)) / (2.0 * h);
        assert!(
            (fd - p.rho_lr_b_prime(u)).abs() < 1e-6,
            "FD mismatch at u = {u}"
        );
    }
}

#[test]
fn rho_lr_d_is_zero() {
    let p = Bi5Reparams::for_valence(5);
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        assert_eq!(p.rho_lr_d(u), 0.0);
    }
}

#[test]
fn rho_lr_e_equals_b_times_bprime() {
    let p = Bi5Reparams::for_valence(8);
    for k in 0..=10 {
        let u = k as f64 / 10.0;
        let b = p.rho_lr_b(u);
        let bp = p.rho_lr_b_prime(u);
        assert!(approx_eq(p.rho_lr_e(u), b * bp));
    }
}

// ---- β (tensor-border) ----

#[test]
fn beta_m_22_matches_eq_15() {
    let p = Bi5Reparams::for_valence(6);
    let c = p.c;
    let (x, y) = p.beta_m_22();
    // Diagonal point by symmetry.
    assert!(approx_eq(x, y));
    assert!(approx_eq(x, 1.0 / 3.0 + c / 216.0));
}

#[test]
fn beta_m_32_matches_eq_15() {
    let p = Bi5Reparams::for_valence(5);
    let c = p.c;
    let (x, y) = p.beta_m_32();
    assert!(approx_eq(x, 0.5 + 7.0 * c / 1296.0));
    assert!(approx_eq(y, 1.0 / 3.0 + c / 72.0));
}

#[test]
fn beta_r_32_matches_eq_15() {
    let p = Bi5Reparams::for_valence(7);
    let c = p.c;
    let (x, y) = p.beta_r_32();
    assert!(approx_eq(x, 1.0));
    assert!(approx_eq(y, 1.0 / 3.0 + c / 18.0));
}

// ---- inter-sector join weights ----

#[test]
fn join_weights_match_eq_12() {
    let p = Bi5Reparams::for_valence(6);
    let c = p.c;
    assert!(approx_eq(p.w1(), -c / 96.0));
    assert!(approx_eq(p.w2(), c / 96.0));
    assert!(approx_eq(p.w3(), c / 30.0));
    // w1 + w2 = 0 (antisymmetric pair).
    assert!(approx_eq(p.w1() + p.w2(), 0.0));
}

#[test]
fn validate_accepts_constructor_output() {
    for n in [3, 4, 5, 6, 7, 8, 9, 10, 12, 20] {
        let p = Bi5Reparams::for_valence(n);
        p.validate().unwrap();
    }
}

#[test]
fn validate_rejects_nonfinite_c() {
    let p = Bi5Reparams { n: 5, c: f64::NAN };
    assert!(p.validate().is_err());
}
