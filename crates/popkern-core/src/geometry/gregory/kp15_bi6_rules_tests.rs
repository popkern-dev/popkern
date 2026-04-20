//! Tests for K-P 2015 bi-6 BB-coefficient assignment rules.
//!
//! Strategy: verify paper-cited identities (ṕ_{00} = p̀_{00}
//! under (S1); constant-input preservation; edge-degenerate
//! rule (5) reduction). Spot-check numerical cases where the
//! formulas simplify.

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::kp15_bi6_mu::{mu_bar, mu_ddot, mu_dot, mu_tdot, mu_tilde};

const EPS: f64 = 1e-10;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_eq(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// Catmull-Clark star-point tests live in `kp_polar_shared_tests.rs`
// since that utility is shared between K-P 2015 and K-P 2013.

// ---- Rule (5) ----

#[test]
fn rule_5_constant_edge() {
    // All six CPs equal P → degree-reduced p̀_{30} = P (affine
    // convex-combination weights sum to 0.05·2 − 0.3·2 + 0.75·2 = 1).
    let p = pt(2.0, 3.0, -1.0);
    let r = rule_5_reduce_edge(p, p, p, p, p, p);
    assert!(approx_eq(r, p));
    // Also check the weights sum to 1.
    let w: f64 = 0.05 + 0.05 + (-0.3) + (-0.3) + 0.75 + 0.75;
    assert!((w - 1.0).abs() < 1e-15);
}

#[test]
fn rule_5_linear_edge_recovers_linear() {
    // If the edge CPs lie on a line at uniform parameter values
    // 0/6, 1/6, 2/6, 3/6, 4/6, 5/6, 6/6 — but since we only have
    // 6 inputs (no p̀_{30}), test with CPs at t_i = i/6 for
    // i ∈ {0, 1, 2, 4, 5, 6}. Then p̀_{30} should be at t=3/6=0.5.
    let t = |i: usize| pt((i as f64) / 6.0, 0.0, 0.0);
    let r = rule_5_reduce_edge(t(0), t(1), t(2), t(4), t(5), t(6));
    // Linear interp at t=0.5 gives x = 0.5.
    assert!((r.x - 0.5).abs() < 1e-12);
}

// ---- Rule (S1) ----

#[test]
fn rule_s1_identity_on_first_three_cps() {
    let p = Bi6Reparams::for_valence(6);
    let p00 = pt(1.0, 0.0, 0.0);
    let p01 = pt(1.1, 0.1, 0.0);
    let p02 = pt(1.2, 0.2, 0.0);
    let p10 = pt(2.0, 0.0, 0.0);
    let p11 = pt(2.1, 0.1, 0.0);
    let p20 = pt(3.0, 0.0, 0.0);
    let out = rule_s1_star_propagation(&p, p00, p01, p02, p10, p11, p20);
    assert_eq!(out.p00, p00);
    assert_eq!(out.p10, p10);
    assert_eq!(out.p20, p20);
}

#[test]
fn rule_s1_p01_formula_spot_check() {
    // ṕ_{01} = 2(1−c)·p̀_{00} + 2c·p̀_{10} − p̀_{01}
    // With p̀_{00} = 0, p̀_{10} = (1, 0, 0), p̀_{01} = 0,
    // ṕ_{01} = 2c·(1, 0, 0).
    let params = Bi6Reparams::for_valence(6);
    let c = params.c;
    let out = rule_s1_star_propagation(
        &params,
        pt(0.0, 0.0, 0.0),
        pt(0.0, 0.0, 0.0),
        pt(0.0, 0.0, 0.0),
        pt(1.0, 0.0, 0.0),
        pt(0.0, 0.0, 0.0),
        pt(0.0, 0.0, 0.0),
    );
    assert!(approx_eq(out.p01, pt(2.0 * c, 0.0, 0.0)));
}

#[test]
fn rule_s1_constant_input_gives_consistent_affine_behaviour() {
    // If all inputs equal P, ṕ_{01} = 2(1−c)P + 2cP − P = P.
    let params = Bi6Reparams::for_valence(6);
    let p = pt(5.0, -3.0, 1.0);
    let out = rule_s1_star_propagation(&params, p, p, p, p, p, p);
    // For ṕ_{01}: weights sum to 2(1−c) + 2c + (−1) = 1 ✓.
    assert!(approx_eq(out.p01, p));
    // For ṕ_{11}: weights sum to c/3 + 2(1−c) + 5c/3 + (−1)
    //                     = 2c + 2 − 2c + (−1) = 1 ✓.
    assert!(approx_eq(out.p11, p));
}

// ---- Rule (S2) ----

#[test]
fn rule_s2_p21_constant_input_gives_p00() {
    // Weight of p̀_{00} in ṕ_{21} = 2 − Σ µ̇_i.
    // Other weights sum to µ̇_1 + µ̇_2 + µ̇_3 + µ̇_4 + µ̇_5 + (−1)
    //   = Σ µ̇ − 1.
    // Total = (2 − Σµ̇) + (Σµ̇ − 1) = 1. Affine.
    // So constant P input gives P.
    let params = Bi6Reparams::for_valence(6);
    let m = mu_dot(&params);
    let p = pt(1.0, 2.0, 3.0);
    let r = rule_s2_p21(m, p, p, p, p, p, p, p);
    assert!(approx_eq(r, p));
}

#[test]
fn rule_s2_p12_constant_input_affine() {
    // ṕ_{12} weights: −Σµ̈ + µ̈_1 + … + µ̈_8 + 1 = (Σµ̈ − Σµ̈) + 1 = 1.
    let params = Bi6Reparams::for_valence(6);
    let m = mu_ddot(&params);
    let p = pt(4.0, -1.0, 2.5);
    let r = rule_s2_p12(m, p, p, p, p, p, p, p, p, p, p);
    assert!(approx_eq(r, p));
}

// ---- Rule (S3) ----

#[test]
fn rule_s3_p22_constant_input_affine_when_diffs_zero() {
    // ṕ_{22} weights: −Σ⃛µ + ⃛µ_1..⃛µ_8 + 1 + (⃛µ_9,10,11)·0 = 1
    // (since diff_42/52/62 are zero when p̀_{i2} = ṕ_{i2}).
    let params = Bi6Reparams::for_valence(6);
    let m = mu_tdot(&params);
    let p = pt(1.0, 2.0, 3.0);
    let zero = pt(0.0, 0.0, 0.0);
    let r = rule_s3_p22(m, p, p, p, p, p, p, p, p, p, p, zero, zero, zero);
    assert!(approx_eq(r, p));
}

// ---- Rule (S4) — i=5,6 simple formulas ----

#[test]
fn rule_s4_i5_simple_formula() {
    // p̀_{51} = p̀_{50} + ¼(p̀_{52} − ṕ_{52})
    // ṕ_{51} = p̀_{50} − ¼(p̀_{52} − ṕ_{52})
    let params = Bi6Reparams::for_valence(6);
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let inp = S4Inputs {
        p_prev_00: pt(0.0, 0.0, 0.0),
        p_prev_10: pt(0.0, 0.0, 0.0),
        p_prev_20: pt(0.0, 0.0, 0.0),
        p_prev_21: pt(0.0, 0.0, 0.0),
        p_prev_40: pt(0.0, 0.0, 0.0),
        p_prev_50: pt(10.0, 0.0, 0.0),
        p_prev_60: pt(0.0, 0.0, 0.0),
        p_cur_21: pt(0.0, 0.0, 0.0),
        p_prev_50_row0: pt(0.0, 0.0, 0.0),
        p_prev_42: pt(0.0, 0.0, 0.0),
        p_prev_52: pt(2.0, 0.0, 0.0),
        p_prev_62: pt(0.0, 0.0, 0.0),
        p_cur_42: pt(0.0, 0.0, 0.0),
        p_cur_52: pt(0.0, 0.0, 0.0),
        p_cur_62: pt(0.0, 0.0, 0.0),
        p_prev_32: pt(0.0, 0.0, 0.0),
    };
    let out = rule_s4(&params, mt, mb, &inp);
    // diff_52 = (2,0,0), ¼·diff_52 = (0.5,0,0)
    // p̀_{51} = (10,0,0) + (0.5,0,0) = (10.5,0,0)
    assert!(approx_eq(out.p_prev_51, pt(10.5, 0.0, 0.0)));
    assert!(approx_eq(out.p_cur_51, pt(9.5, 0.0, 0.0)));
}

#[test]
fn rule_s4_i6_simple_formula() {
    let params = Bi6Reparams::for_valence(6);
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let inp = S4Inputs {
        p_prev_00: pt(0.0, 0.0, 0.0),
        p_prev_10: pt(0.0, 0.0, 0.0),
        p_prev_20: pt(0.0, 0.0, 0.0),
        p_prev_21: pt(0.0, 0.0, 0.0),
        p_prev_40: pt(0.0, 0.0, 0.0),
        p_prev_50: pt(0.0, 0.0, 0.0),
        p_prev_60: pt(20.0, 0.0, 0.0),
        p_cur_21: pt(0.0, 0.0, 0.0),
        p_prev_50_row0: pt(0.0, 0.0, 0.0),
        p_prev_42: pt(0.0, 0.0, 0.0),
        p_prev_52: pt(0.0, 0.0, 0.0),
        p_prev_62: pt(4.0, 0.0, 0.0),
        p_cur_42: pt(0.0, 0.0, 0.0),
        p_cur_52: pt(0.0, 0.0, 0.0),
        p_cur_62: pt(0.0, 0.0, 0.0),
        p_prev_32: pt(0.0, 0.0, 0.0),
    };
    let out = rule_s4(&params, mt, mb, &inp);
    // diff_62 = (4,0,0); ¼·diff_62 = (1,0,0)
    // p̀_{61} = (20,0,0) + (1,0,0) = (21,0,0)
    assert!(approx_eq(out.p_prev_61, pt(21.0, 0.0, 0.0)));
    assert!(approx_eq(out.p_cur_61, pt(19.0, 0.0, 0.0)));
}

#[test]
fn rule_s4_p41_zero_diffs_reduces_to_p40_plus_w3_term() {
    // With all diff_i2 = 0 and p̀_{60} = p̀_{50} + (1, 0, 0):
    //   p̀_{41} = p̀_{40} + 0 + 0 + 0 + w₃·(1, 0, 0)
    //   ṕ_{41} = p̀_{40} + w₃·(1, 0, 0)  (same, since all diff terms
    //             cancel with sign flips)
    let params = Bi6Reparams::for_valence(6);
    let c = params.c;
    let w3 = c / 15.0;
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let inp = S4Inputs {
        p_prev_00: pt(0.0, 0.0, 0.0),
        p_prev_10: pt(0.0, 0.0, 0.0),
        p_prev_20: pt(0.0, 0.0, 0.0),
        p_prev_21: pt(0.0, 0.0, 0.0),
        p_prev_40: pt(5.0, 0.0, 0.0),
        p_prev_50: pt(6.0, 0.0, 0.0),
        p_prev_60: pt(7.0, 0.0, 0.0),
        p_cur_21: pt(0.0, 0.0, 0.0),
        p_prev_50_row0: pt(0.0, 0.0, 0.0),
        p_prev_42: pt(0.0, 0.0, 0.0),
        p_prev_52: pt(0.0, 0.0, 0.0),
        p_prev_62: pt(0.0, 0.0, 0.0),
        p_cur_42: pt(0.0, 0.0, 0.0),
        p_cur_52: pt(0.0, 0.0, 0.0),
        p_cur_62: pt(0.0, 0.0, 0.0),
        p_prev_32: pt(0.0, 0.0, 0.0),
    };
    let out = rule_s4(&params, mt, mb, &inp);
    // p̀_{60} − p̀_{50} = (1, 0, 0), w₃·(1,0,0) = (w3,0,0).
    assert!(approx_eq(out.p_prev_41, pt(5.0 + w3, 0.0, 0.0)));
    assert!(approx_eq(out.p_cur_41, pt(5.0 + w3, 0.0, 0.0)));
}
