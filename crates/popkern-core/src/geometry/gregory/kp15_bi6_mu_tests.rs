//! Tests for the K-P 2015 bi-6 µ coefficient tables.
//!
//! Strategy: the µ vectors are pure rational functions of
//! `(c, τ)`. At `n = 6` we have `c = cos(π/3) = 0.5` (exact), and
//! `τ = 0.85` (paper-supplied), so every entry has a closed-form
//! value computable from the formula. Tests spot-check each
//! vector at `n = 6` by evaluating the paper's formula
//! independently and comparing.

#![allow(clippy::approx_constant)]

use super::*;
use crate::geometry::gregory::kp15_bi6_reparam::Bi6Reparams;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

// ---- µ̇ ----

#[test]
fn mu_dot_explicit_formula_at_n6() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    assert!(approx_eq(c, 0.5)); // cos(π/3) ≈ 0.5 (FP rounding).
    let m = mu_dot(&p);
    assert!(approx_eq(m[0], 2.0 * c / 5.0));
    assert!(approx_eq(m[1], 2.0 - c));
    assert!(approx_eq(m[2], c));
    assert!(approx_eq(m[3], -2.0 * c / 5.0));
    assert!(approx_eq(m[4], c / 15.0));
    // Concrete numbers for the record.
    assert!(approx_eq(m[0], 0.2));
    assert!(approx_eq(m[1], 1.5));
    assert!(approx_eq(m[2], 0.5));
    assert!(approx_eq(m[3], -0.2));
    assert!(approx_eq(m[4], 0.5 / 15.0));
}

#[test]
fn mu_dot_at_n5() {
    let p = Bi6Reparams::for_valence(5);
    let c = p.c;
    let m = mu_dot(&p);
    assert!(approx_eq(m[0], 2.0 * c / 5.0));
    assert!(approx_eq(m[1], 2.0 - c));
    assert!(approx_eq(m[2], c));
}

// ---- µ̈ ----

#[test]
fn mu_ddot_explicit_formula_at_n6() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let tau = p.tau;
    let t2 = tau * tau;
    let m = mu_ddot(&p);
    assert!(approx_eq(m[0], -4.0 * c * t2 / 5.0));
    assert!(approx_eq(
        m[1],
        4.0 * (1.0 - c + c * c / 5.0) + 24.0 * c * (c - 1.0) / 5.0 * tau + 4.0 * c * c * t2 / 5.0
    ));
    assert!(approx_eq(m[2], 24.0 * c / 5.0 * tau - 4.0));
    assert!(approx_eq(
        m[3],
        2.0 * c / 3.0 * (11.0 - 3.0 * c - 6.0 * c * tau)
    ));
    assert!(approx_eq(m[4], -4.0 * c));
    assert!(approx_eq(m[5], 2.0 * c * c));
    assert!(approx_eq(m[6], -4.0 * c * c / 5.0));
    assert!(approx_eq(m[7], 2.0 * c * c / 15.0));
}

// ---- ⃛µ ----

#[test]
fn mu_tdot_explicit_formula_at_n6() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let tau = p.tau;
    let t2 = tau * tau;
    let t3 = t2 * tau;
    let m = mu_tdot(&p);
    assert!(approx_eq(m[0], 4.0 * c * t3 / 25.0));
    assert!(approx_eq(
        m[1],
        4.0 * c / 25.0 * (5.0 - 3.0 * c * tau + 6.0 * (1.0 - c) * t2 - c * t3)
    ));
    assert!(approx_eq(m[2], -24.0 * c * t2 / 25.0));
    assert!(approx_eq(
        m[3],
        4.0 - 2.0 * c + 6.0 * c * (c - 2.0) / 5.0 * tau + 4.0 * c * c / 5.0 * t2
    ));
    assert!(approx_eq(m[4], 12.0 * c / 5.0 * tau - 4.0));
    assert!(approx_eq(m[5], 2.0 * c / 5.0 * (5.0 - 3.0 * c * tau)));
    assert!(approx_eq(m[6], 4.0 * c / 25.0 * (3.0 * c * tau - 5.0)));
    assert!(approx_eq(m[7], 2.0 * c / 75.0 * (5.0 - 3.0 * c * tau)));
    assert!(approx_eq(m[8], -3.0 * c / (5.0 * tau)));
    assert!(approx_eq(m[9], 6.0 * c * (5.0 - c * tau) / (125.0 * t2)));
    assert!(approx_eq(m[10], c * (6.0 * c * tau - 5.0) / (125.0 * t3)));
}

// ---- µ̃ ----

#[test]
fn mu_tilde_explicit_formula_at_n6() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let tau = p.tau;
    let t2 = tau * tau;
    let t3 = t2 * tau;
    let m = mu_tilde(&p);
    // µ̃_1 is a paper-pinned constant.
    assert_eq!(m[0], -11.0 / 20.0);
    assert!(approx_eq(m[1], (11.0 * c - 10.0) / (8.0 * c)));
    assert!(approx_eq(m[2], 5.0 / (4.0 * c)));
    // µ̃_4 at c = 0.5: 1/8 − 0.5/4 = 0.125 − 0.125 = 0 exactly.
    assert!(approx_eq(m[3], 0.0));
    assert!(approx_eq(m[4], 3.0 * c / 10.0 - 1.0 / 20.0));
    assert!(approx_eq(m[5], 1.0 / 120.0 - c / 20.0));
    assert!(approx_eq(m[6], 3.0 / (8.0 * tau)));
    assert!(approx_eq(m[7], (3.0 * c * tau - 15.0) / (100.0 * t2)));
    assert!(approx_eq(m[8], (5.0 - 6.0 * c * tau) / (200.0 * t3)));
    assert!(approx_eq(m[9], -5.0 / (16.0 * c)));
    // Concrete: µ̃_10 at c = 0.5 = −5/8.
    assert!(approx_eq(m[9], -5.0 / 8.0));
}

// ---- µ̄ ----

#[test]
fn mu_bar_explicit_formula_at_n6() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let tau = p.tau;
    let t2 = tau * tau;
    let t3 = t2 * tau;
    let m = mu_bar(&p);
    assert_eq!(m[0], 1.0);
    assert!(approx_eq(m[1], 5.0 * (2.0 - c) / (2.0 * c)));
    assert!(approx_eq(m[2], -5.0 / c));
    assert_eq!(m[3], 5.0 / 2.0);
    assert_eq!(m[4], -1.0);
    assert!(approx_eq(m[5], 1.0 / 6.0));
    assert!(approx_eq(m[6], 3.0 / (10.0 * tau) * (c * tau - 5.0)));
    assert!(approx_eq(
        m[7],
        3.0 / (125.0 * t2) * (25.0 - 20.0 * c * tau + c * c * t2)
    ));
    assert!(approx_eq(
        m[8],
        -1.0 / (250.0 * t3) * (25.0 - 45.0 * c * tau + 6.0 * c * c * t2)
    ));
    assert!(approx_eq(m[9], 5.0 / (4.0 * c)));
    // Concrete at n=6: µ̄_2 = 5·1.5/1 = 7.5, µ̄_3 = −10, µ̄_10 = 2.5.
    assert!(approx_eq(m[1], 7.5));
    assert!(approx_eq(m[2], -10.0));
    assert!(approx_eq(m[9], 2.5));
}

// ---- w scalars ----

#[test]
fn w_scalars_explicit_formula_at_n6() {
    let p = Bi6Reparams::for_valence(6);
    let c = p.c;
    let tau = p.tau;
    let (w1, w2, w3) = w_scalars(&p);
    assert!(approx_eq(w1, -c / (50.0 * tau)));
    assert!(approx_eq(w2, c / 50.0));
    assert!(approx_eq(w3, c / 15.0));
    // Concrete: w₂ at c=0.5 = 0.01.
    assert!(approx_eq(w2, 0.01));
    // w₃ same as µ̇_5 at n=6.
    assert!(approx_eq(w3, mu_dot(&p)[4]));
}

// ---- cross-valence smoke ----

#[test]
fn all_mu_entries_finite_for_standard_valences() {
    for n in [3usize, 5, 6, 7, 8, 9, 10, 11, 12, 20] {
        let p = Bi6Reparams::for_valence(n);
        for (label, entry) in mu_dot(&p).iter().enumerate() {
            assert!(entry.is_finite(), "n={n} µ̇_{} not finite", label + 1);
        }
        for (label, entry) in mu_ddot(&p).iter().enumerate() {
            assert!(entry.is_finite(), "n={n} µ̈_{} not finite", label + 1);
        }
        for (label, entry) in mu_tdot(&p).iter().enumerate() {
            assert!(entry.is_finite(), "n={n} ⃛µ_{} not finite", label + 1);
        }
        for (label, entry) in mu_tilde(&p).iter().enumerate() {
            assert!(entry.is_finite(), "n={n} µ̃_{} not finite", label + 1);
        }
        for (label, entry) in mu_bar(&p).iter().enumerate() {
            assert!(entry.is_finite(), "n={n} µ̄_{} not finite", label + 1);
        }
        let (w1, w2, w3) = w_scalars(&p);
        assert!(w1.is_finite() && w2.is_finite() && w3.is_finite());
    }
}

#[test]
fn mu_bar_first_entry_is_always_one() {
    // The paper's µ̄_1 = 1 is a constant independent of (c, τ).
    for n in [3usize, 5, 6, 7, 8, 12, 20] {
        let p = Bi6Reparams::for_valence(n);
        assert_eq!(mu_bar(&p)[0], 1.0, "n={n}");
        assert_eq!(mu_bar(&p)[3], 2.5, "n={n}");
        assert_eq!(mu_bar(&p)[4], -1.0, "n={n}");
    }
}

#[test]
fn mu_tilde_first_entry_is_paper_constant() {
    // µ̃_1 = −11/20 — independent of (c, τ).
    for n in [3usize, 5, 6, 7, 8, 12, 20] {
        let p = Bi6Reparams::for_valence(n);
        assert_eq!(mu_tilde(&p)[0], -11.0 / 20.0);
    }
}

#[test]
fn mu_dot_scales_linearly_in_c() {
    // µ̇_1 = 2c/5, µ̇_3 = c, µ̇_4 = −2c/5, µ̇_5 = c/15 — all linear
    // in c. For two valences with different c, ratios match.
    let p5 = Bi6Reparams::for_valence(5);
    let p6 = Bi6Reparams::for_valence(6);
    let m5 = mu_dot(&p5);
    let m6 = mu_dot(&p6);
    let ratio = p5.c / p6.c;
    assert!(approx_eq(m5[0] / m6[0], ratio));
    assert!(approx_eq(m5[2] / m6[2], ratio));
    assert!(approx_eq(m5[3] / m6[3], ratio));
    assert!(approx_eq(m5[4] / m6[4], ratio));
}
