//! Tests for the σ F_3 solver.

use super::super::kp13_bi5_sigma::{SigmaPatch, f3_functional};
use super::*;

const EPS: f64 = 1e-9;

// ---- Flat-index encoding ----

#[test]
fn decode_flat_zero_is_zero_patch() {
    let x = [0.0_f64; SIGMA_FLAT_LEN];
    let patch = decode_flat(&x);
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(patch.cps[i][j], (0.0, 0.0), "cp[{i}][{j}]");
        }
    }
}

#[test]
fn decode_flat_respects_symmetry() {
    // Set one off-diagonal scalar and verify the swap partner is also set.
    let mut x = [0.0_f64; SIGMA_FLAT_LEN];
    let (ix, iy) = off_diag_flat_index(2, 4);
    x[ix] = 1.0;
    x[iy] = 2.0;
    let patch = decode_flat(&x);
    assert_eq!(patch.cps[2][4], (1.0, 2.0));
    assert_eq!(patch.cps[4][2], (2.0, 1.0)); // swap
    assert!(patch.is_symmetric());
}

#[test]
fn off_diag_flat_index_covers_30_scalars() {
    // The 15 pairs give flat indices 6..36.
    let mut used = std::collections::HashSet::new();
    for i in 0..6 {
        for j in (i + 1)..6 {
            let (ix, iy) = off_diag_flat_index(i, j);
            assert!(ix < SIGMA_FLAT_LEN);
            assert!(iy < SIGMA_FLAT_LEN);
            assert_eq!(iy, ix + 1);
            used.insert(ix);
            used.insert(iy);
        }
    }
    assert_eq!(used.len(), 30);
}

// ---- Hessian properties ----

#[test]
fn hessian_is_symmetric() {
    let h = build_f3_hessian();
    #[allow(clippy::needless_range_loop)]
    for i in 0..SIGMA_FLAT_LEN {
        for j in (i + 1)..SIGMA_FLAT_LEN {
            assert!(
                (h[i][j] - h[j][i]).abs() < 1e-10,
                "H[{i}][{j}] = {} vs H[{j}][{i}] = {}",
                h[i][j],
                h[j][i]
            );
        }
    }
}

#[test]
fn hessian_reproduces_f3_via_quadratic_form() {
    // F_3(x) = ½ x^T H x. Pick a random-ish x and verify.
    let h = build_f3_hessian();
    let x = {
        let mut x = [0.0_f64; SIGMA_FLAT_LEN];
        for (i, xi) in x.iter_mut().enumerate() {
            // Deterministic pseudo-random: f64 hash of index.
            *xi = ((i as f64 * 37.0).sin() * 0.5 + 0.3) * 0.1;
        }
        x
    };
    let mut quad = 0.0_f64;
    for i in 0..SIGMA_FLAT_LEN {
        for j in 0..SIGMA_FLAT_LEN {
            quad += x[i] * h[i][j] * x[j];
        }
    }
    quad *= 0.5;

    let patch = decode_flat(&x);
    let direct = f3_functional(&patch);
    assert!(
        (quad - direct).abs() < 1e-6,
        "½ x^T H x = {quad}, F_3(x) direct = {direct}"
    );
}

#[test]
fn hessian_positive_semi_definite() {
    // F_3 ≥ 0 implies H is PSD. Sample a few vectors and check
    // x^T H x ≥ 0.
    let h = build_f3_hessian();
    for seed in 0..20 {
        let mut x = [0.0_f64; SIGMA_FLAT_LEN];
        for (i, xi) in x.iter_mut().enumerate() {
            *xi = ((i as f64 * 13.0 + seed as f64 * 7.0).sin() * 0.7).tanh();
        }
        let mut q = 0.0;
        for i in 0..SIGMA_FLAT_LEN {
            for j in 0..SIGMA_FLAT_LEN {
                q += x[i] * h[i][j] * x[j];
            }
        }
        assert!(q >= -1e-10, "seed {seed}: x^T H x = {q} < 0");
    }
}

// ---- Solver ----

#[test]
fn solve_returns_symmetric_patch() {
    let patch = solve_sigma_f3_minimal();
    assert!(patch.is_symmetric());
}

#[test]
fn solve_pins_corners() {
    let patch = solve_sigma_f3_minimal();
    // Corners pinned by construction.
    let c00 = patch.cps[0][0];
    let c50 = patch.cps[5][0];
    let c05 = patch.cps[0][5];
    let c55 = patch.cps[5][5];
    assert!((c00.0 - 0.0).abs() < EPS && (c00.1 - 0.0).abs() < EPS);
    assert!((c55.0 - 1.0).abs() < EPS && (c55.1 - 1.0).abs() < EPS);
    assert!((c05.0 - 0.0).abs() < EPS && (c05.1 - 1.0).abs() < EPS);
    assert!((c50.0 - 1.0).abs() < EPS && (c50.1 - 0.0).abs() < EPS);
}

#[test]
fn solve_beats_identity_on_f3() {
    // F_3(solved) ≤ F_3(identity). With the same corner constraints,
    // the F_3-minimiser strictly improves over identity (unless
    // identity is already optimal, which it isn't — identity has
    // non-zero (u-v)^k σ contributions).
    let solved = solve_sigma_f3_minimal();
    let identity = SigmaPatch::identity();
    let f_solved = f3_functional(&solved);
    let f_identity = f3_functional(&identity);
    assert!(
        f_solved <= f_identity + 1e-8,
        "solved F_3 = {f_solved}, identity F_3 = {f_identity}"
    );
}

#[test]
fn solve_eval_at_corners_matches_expected() {
    let patch = solve_sigma_f3_minimal();
    // Eval at corners should hit the pinned values.
    let v00 = patch.eval(0.0, 0.0);
    let v10 = patch.eval(1.0, 0.0);
    let v01 = patch.eval(0.0, 1.0);
    let v11 = patch.eval(1.0, 1.0);
    assert!((v00.0 - 0.0).abs() < EPS && (v00.1 - 0.0).abs() < EPS);
    assert!((v10.0 - 1.0).abs() < EPS && (v10.1 - 0.0).abs() < EPS);
    assert!((v01.0 - 0.0).abs() < EPS && (v01.1 - 1.0).abs() < EPS);
    assert!((v11.0 - 1.0).abs() < EPS && (v11.1 - 1.0).abs() < EPS);
}

#[test]
fn solve_is_deterministic() {
    let a = solve_sigma_f3_minimal();
    let b = solve_sigma_f3_minimal();
    for i in 0..6 {
        for j in 0..6 {
            assert!(
                (a.cps[i][j].0 - b.cps[i][j].0).abs() < 1e-12
                    && (a.cps[i][j].1 - b.cps[i][j].1).abs() < 1e-12
            );
        }
    }
}
