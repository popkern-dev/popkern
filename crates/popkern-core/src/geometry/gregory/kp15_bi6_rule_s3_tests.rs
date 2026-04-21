//! Tests for K-P 2015 bi-6 Rule (S3) — paper eq (9).
//!
//! Verifies the full µ̃_{10}/µ̄_{10} wiring (the `(p̀_{22} −
//! ṕ_{22})` term that was missing in the pre-audit revision).

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::kp15_bi6_mu::{mu_bar, mu_tilde};
use crate::geometry::gregory::kp15_bi6_reparam::Bi6Reparams;

const EPS: f64 = 1e-10;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_eq(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

fn s3_inputs_constant(p: Point3) -> S3Inputs {
    S3Inputs {
        p_prev_00: p,
        p_prev_10: p,
        p_prev_20: p,
        p_prev_40: p,
        p_prev_50: p,
        p_prev_60: p,
        p_prev_21: p,
        p_cur_21: p,
        p_prev_22: p,
        p_cur_22: p,
        p_prev_42: p,
        p_cur_42: p,
        p_prev_52: p,
        p_cur_52: p,
        p_prev_62: p,
        p_cur_62: p,
        p_prev_32: p,
    }
}

#[test]
fn constant_input_with_zero_diffs_is_fixed_point() {
    // Weights for p̀_{31}: (1 − Σµ̃_{1..6}) + µ̃_{1..6} + 0·diffs = 1
    // Weights for ṕ_{31}: same, diffs zero so sign flip irrelevant
    // Weights for ṕ_{32}: 1·p̀_{32} + (−Σµ̄_{1..6}) + µ̄_{1..6} + 0·diffs = 1
    let params = Bi6Reparams::for_valence(6);
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let p = pt(1.0, 2.0, -1.0);
    let inp = s3_inputs_constant(p);
    let out = rule_s3(mt, mb, &inp);
    assert!(approx_eq(out.p_prev_31, p));
    assert!(approx_eq(out.p_cur_31, p));
    assert!(approx_eq(out.p_cur_32, p));
}

#[test]
fn mu_tilde_10_term_present_in_p_prev_31_minus_p_cur_31() {
    // With zero CPs except p̀_{22} = (1,0,0), diff_22 = (1,0,0).
    //   p̀_{31}: + µ̃_{10}·(1,0,0)
    //   ṕ_{31}: − µ̃_{10}·(1,0,0)
    //   delta = 2·µ̃_{10}·(1,0,0)
    // Other diffs are zero so µ̃_{7,8,9} contribute nothing.
    let params = Bi6Reparams::for_valence(6);
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let zero = pt(0.0, 0.0, 0.0);
    let mut inp = s3_inputs_constant(zero);
    inp.p_prev_22 = pt(1.0, 0.0, 0.0);
    let out = rule_s3(mt, mb, &inp);
    let delta_x = out.p_prev_31.x - out.p_cur_31.x;
    let expected = 2.0 * mt[9];
    assert!(
        (delta_x - expected).abs() < 1e-12,
        "delta_x = {delta_x}, expected = {expected}"
    );
}

#[test]
fn mu_bar_10_term_present_in_p_cur_32() {
    // Zero everywhere except p̀_{22} = (1,0,0) → diff_22 = (1,0,0).
    // ṕ_{32} picks up µ̄_{10}·diff_22.
    let params = Bi6Reparams::for_valence(6);
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let zero = pt(0.0, 0.0, 0.0);
    let mut inp = s3_inputs_constant(zero);
    inp.p_prev_22 = pt(1.0, 0.0, 0.0);
    let out = rule_s3(mt, mb, &inp);
    assert!(
        (out.p_cur_32.x - mb[9]).abs() < 1e-12,
        "p_cur_32.x = {}, expected µ̄_10 = {}",
        out.p_cur_32.x,
        mb[9]
    );
}

#[test]
fn mu_tilde_7_term_present_in_diff_42() {
    // Zero everywhere except p̀_{42} = (1,0,0) → diff_42 = (1,0,0).
    // p̀_{31} = +µ̃_{7}·(1,0,0); ṕ_{31} = −µ̃_{7}·(1,0,0).
    let params = Bi6Reparams::for_valence(6);
    let mt = mu_tilde(&params);
    let mb = mu_bar(&params);
    let zero = pt(0.0, 0.0, 0.0);
    let mut inp = s3_inputs_constant(zero);
    inp.p_prev_42 = pt(1.0, 0.0, 0.0);
    let out = rule_s3(mt, mb, &inp);
    assert!((out.p_prev_31.x - mt[6]).abs() < 1e-12);
    assert!((out.p_cur_31.x - (-mt[6])).abs() < 1e-12);
}
