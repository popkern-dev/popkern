//! Tests for KP09a § 6 low-valence machinery.

use super::super::kp13_bi5_boundary::TensorBorder;
use super::super::kp13_bi5_circulant::{circulant_residual, solve_g21_circulant};
use super::super::kp13_bi5_free_params::{Bi5FreeParams, apply_free_params};
use super::super::kp13_bi5_gtri5::Gtri5;
use super::super::kp13_bi5_lsq::solve_g_tri_lsq;
use super::super::kp13_bi5_lsq_cc::solve_g_tri_lsq_cc;
use super::*;
use crate::geometry::Point3;

// ---- singular_modes lookup ----

#[test]
fn singular_modes_n3() {
    assert_eq!(singular_modes(3), vec![0]);
    assert_eq!(null_space_dim(3), 1);
}

#[test]
fn singular_modes_n4() {
    assert_eq!(singular_modes(4), vec![1, 3]);
    assert_eq!(null_space_dim(4), 2);
}

#[test]
fn singular_modes_n5_through_n8() {
    for n in [5, 7, 8] {
        assert!(singular_modes(n).is_empty(), "n={n} should be full-rank");
        assert_eq!(null_space_dim(n), 0);
    }
    assert_eq!(singular_modes(6), vec![3]);
    assert_eq!(null_space_dim(6), 1);
}

#[test]
fn singular_modes_monotone_pattern() {
    // Beyond n = 8 the circulant is always full-rank for the
    // KP13 denominator: 2c + cos(jα) = 0 only when cos(2π/n) is
    // a specific small set of values.
    for n in 9..=16 {
        assert_eq!(null_space_dim(n), 0, "n={n} should be full-rank");
    }
}

// ---- End-to-end: the full pipeline works for n = 3 and n = 4 ----

#[test]
fn apply_free_params_works_for_n3() {
    let params = Bi5FreeParams::zero(3);
    let star = Point3::new(0.5, -0.2, 1.3);
    let g = apply_free_params(&params, star).expect("n=3 supported");
    g.validate().expect("n=3 g△ valid after pseudo-inverse");
}

#[test]
fn apply_free_params_works_for_n4() {
    let params = Bi5FreeParams::zero(4);
    let star = Point3::new(0.0, 1.0, 0.0);
    let g = apply_free_params(&params, star).expect("n=4 supported");
    g.validate().expect("n=4 g△ valid after pseudo-inverse");
}

#[test]
fn solve_bilinear_works_for_n3_n4() {
    for n in [3, 4] {
        let target = TensorBorder::zero(n);
        let g = solve_g_tri_lsq(&target, Point3::new(0.0, 0.0, 0.0))
            .expect("bilinear LSQ works at low valence");
        assert_eq!(g.n, n);
    }
}

#[test]
fn solve_cc_works_for_n3_n4() {
    for n in [3, 4] {
        let target = TensorBorder::zero(n);
        let g = solve_g_tri_lsq_cc(&target, Point3::new(0.0, 0.0, 0.0))
            .expect("CC LSQ works at low valence");
        assert_eq!(g.n, n);
    }
}

// ---- Pseudo-inverse residual properties ----

#[test]
fn n3_pseudo_inverse_residual_in_null_space() {
    // For n = 3, the null mode is j = 0 (constant). A deliberately
    // constant-mode-rich RHS leaves residual proportional to
    // (1, 1, 1) — the constant pattern.
    let mut g = Gtri5::new(3);
    // Fill all sectors identically → compute_r_i returns the same
    // value for all i → R ∈ null space (constant mode).
    use super::super::kp13_bi5_connect::paper_flat;
    let idx_30 = paper_flat(3, 0);
    let idx_03 = paper_flat(0, 3);
    for k in 0..3 {
        g.sector_mut(k).cps[idx_30] = Point3::new(1.0, 0.0, 0.0);
        g.sector_mut(k).cps[idx_03] = Point3::new(1.0, 0.0, 0.0);
    }
    solve_g21_circulant(&mut g).unwrap();
    // All residuals equal: confirms the residual is in the constant
    // (null) mode.
    let r0 = circulant_residual(&g, 0);
    for i in 1..3 {
        let ri = circulant_residual(&g, i);
        assert!(
            (ri.x - r0.x).abs() < 1e-10 && (ri.y - r0.y).abs() < 1e-10,
            "n=3 residual at sector {i} = {ri:?} differs from sector 0 = {r0:?}"
        );
    }
}

#[test]
fn n4_pseudo_inverse_residual_in_null_space() {
    // For n = 4, null modes are j ∈ {1, 3} (alternating). A
    // constant RHS is orthogonal to these, so pseudo-inverse
    // gives an exact solution.
    let mut g = Gtri5::new(4);
    let p = Point3::new(1.0, 2.0, 3.0);
    for k in 0..4 {
        g.sector_mut(k).cps = [p; 21];
    }
    solve_g21_circulant(&mut g).unwrap();
    // All residuals should be zero (RHS has no null-mode content).
    for i in 0..4 {
        let r = circulant_residual(&g, i);
        assert!(
            r.x.abs() < 1e-9 && r.y.abs() < 1e-9 && r.z.abs() < 1e-9,
            "n=4 constant RHS: residual at sector {i} = {r:?}"
        );
    }
}

// ---- Determinism ----

#[test]
fn low_valence_solve_deterministic() {
    for n in [3, 4] {
        let params = Bi5FreeParams::zero(n);
        let star = Point3::new(0.1, 0.2, 0.3);
        let g_a = apply_free_params(&params, star).unwrap();
        let g_b = apply_free_params(&params, star).unwrap();
        for k in 0..n {
            for (cp_a, cp_b) in g_a.sector(k).cps.iter().zip(g_b.sector(k).cps.iter()) {
                assert!(
                    (cp_a.x - cp_b.x).abs() < 1e-14
                        && (cp_a.y - cp_b.y).abs() < 1e-14
                        && (cp_a.z - cp_b.z).abs() < 1e-14,
                    "n={n}, k={k}: non-deterministic CP {cp_a:?} vs {cp_b:?}"
                );
            }
        }
    }
}
