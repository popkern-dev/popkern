//! Tests for the triangular total-degree-5 Bernstein basis.
//!
//! Cross-check: partition of unity, endpoint (barycentric-corner)
//! deltas, finite-difference derivative parity, and flat-index
//! correctness.

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

/// Small deterministic sample inside the open simplex.
fn interior_samples() -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    // Keep `u + v < 1` strictly so w > 0.
    for ui in 1..9 {
        for vi in 1..(9 - ui) {
            out.push((ui as f64 / 10.0, vi as f64 / 10.0));
        }
    }
    out
}

#[test]
fn cp_count_helpers() {
    assert_eq!(tri_cp_count(3), 10);
    assert_eq!(tri_cp_count(4), 15);
    assert_eq!(tri_cp_count(5), 21);
    assert_eq!(TRI5_CP_COUNT, 21);
    assert_eq!(TRI4_CP_COUNT, 15);
    assert_eq!(TRI3_CP_COUNT, 10);
}

#[test]
fn flat_index_bijection_degree_5() {
    let mut seen = [false; 21];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let idx = flat_index(5, i, j);
            assert!(idx < 21);
            assert!(!seen[idx], "duplicate flat index for (i={i}, j={j})");
            seen[idx] = true;
        }
    }
    assert!(seen.iter().all(|&x| x));
}

#[test]
fn flat_index_specific_values_degree_5() {
    // Sanity against the docstring formula.
    assert_eq!(flat_index(5, 0, 0), 0);
    assert_eq!(flat_index(5, 0, 5), 5);
    assert_eq!(flat_index(5, 1, 0), 6);
    assert_eq!(flat_index(5, 1, 4), 10);
    assert_eq!(flat_index(5, 2, 0), 11);
    assert_eq!(flat_index(5, 3, 0), 15);
    assert_eq!(flat_index(5, 4, 0), 18);
    assert_eq!(flat_index(5, 5, 0), 20);
}

#[test]
fn partition_of_unity() {
    for &(u, v) in &interior_samples() {
        let b = bernstein_tri5_all(u, v);
        let sum: f64 = b.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "(u,v)=({u},{v}), sum={sum}");
    }
    // Also at the three corners and on the edges.
    for (u, v) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)] {
        let b = bernstein_tri5_all(u, v);
        let sum: f64 = b.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }
}

#[test]
fn corner_delta_basis() {
    // At u=1, v=0 (corner (5, 0, 0)): only B^5_{5,0,0} = 1.
    let b = bernstein_tri5_all(1.0, 0.0);
    let target = flat_index(5, 5, 0);
    for (i, &v) in b.iter().enumerate() {
        let expected = if i == target { 1.0 } else { 0.0 };
        assert!(
            approx_eq(v, expected),
            "idx {i}: got {v}, expected {expected}"
        );
    }
    // At u=0, v=1 (corner (0, 5, 0)): only B^5_{0,5,0} = 1.
    let b = bernstein_tri5_all(0.0, 1.0);
    let target = flat_index(5, 0, 5);
    for (i, &v) in b.iter().enumerate() {
        let expected = if i == target { 1.0 } else { 0.0 };
        assert!(approx_eq(v, expected));
    }
    // At u=0, v=0 (corner (0, 0, 5)): only B^5_{0,0,5} = 1.
    let b = bernstein_tri5_all(0.0, 0.0);
    let target = flat_index(5, 0, 0);
    for (i, &v) in b.iter().enumerate() {
        let expected = if i == target { 1.0 } else { 0.0 };
        assert!(approx_eq(v, expected));
    }
}

#[test]
fn partial_u_matches_finite_difference() {
    let h = 1e-6;
    for &(u, v) in &interior_samples() {
        let analytic = bernstein_tri5_partial_u(u, v);
        let plus = bernstein_tri5_all(u + h, v);
        let minus = bernstein_tri5_all(u - h, v);
        for i in 0..TRI5_CP_COUNT {
            let fd = (plus[i] - minus[i]) / (2.0 * h);
            assert!(
                (analytic[i] - fd).abs() < 1e-5,
                "idx {i}, (u,v)=({u},{v}): analytic = {}, fd = {fd}",
                analytic[i]
            );
        }
    }
}

#[test]
fn partial_v_matches_finite_difference() {
    let h = 1e-6;
    for &(u, v) in &interior_samples() {
        let analytic = bernstein_tri5_partial_v(u, v);
        let plus = bernstein_tri5_all(u, v + h);
        let minus = bernstein_tri5_all(u, v - h);
        for i in 0..TRI5_CP_COUNT {
            let fd = (plus[i] - minus[i]) / (2.0 * h);
            assert!(
                (analytic[i] - fd).abs() < 1e-5,
                "idx {i}, (u,v)=({u},{v}): analytic = {}, fd = {fd}",
                analytic[i]
            );
        }
    }
}

#[test]
fn second_partial_uu_matches_finite_difference() {
    let h = 1e-4;
    for &(u, v) in &interior_samples() {
        let analytic = bernstein_tri5_second_partial_uu(u, v);
        let plus = bernstein_tri5_all(u + h, v);
        let zero = bernstein_tri5_all(u, v);
        let minus = bernstein_tri5_all(u - h, v);
        for i in 0..TRI5_CP_COUNT {
            let fd = (plus[i] - 2.0 * zero[i] + minus[i]) / (h * h);
            assert!(
                (analytic[i] - fd).abs() < 1e-3,
                "idx {i}, (u,v)=({u},{v}): analytic = {}, fd = {fd}",
                analytic[i]
            );
        }
    }
}

#[test]
fn second_partial_vv_matches_finite_difference() {
    let h = 1e-4;
    for &(u, v) in &interior_samples() {
        let analytic = bernstein_tri5_second_partial_vv(u, v);
        let plus = bernstein_tri5_all(u, v + h);
        let zero = bernstein_tri5_all(u, v);
        let minus = bernstein_tri5_all(u, v - h);
        for i in 0..TRI5_CP_COUNT {
            let fd = (plus[i] - 2.0 * zero[i] + minus[i]) / (h * h);
            assert!(
                (analytic[i] - fd).abs() < 1e-3,
                "idx {i}, (u,v)=({u},{v}): analytic = {}, fd = {fd}",
                analytic[i]
            );
        }
    }
}

#[test]
fn second_partial_uv_matches_finite_difference() {
    // ∂²B/∂u∂v ≈ (B(u+h,v+h) - B(u+h,v-h) - B(u-h,v+h) + B(u-h,v-h)) / (4 h²).
    let h = 1e-4;
    for &(u, v) in &interior_samples() {
        let analytic = bernstein_tri5_second_partial_uv(u, v);
        let pp = bernstein_tri5_all(u + h, v + h);
        let pm = bernstein_tri5_all(u + h, v - h);
        let mp = bernstein_tri5_all(u - h, v + h);
        let mm = bernstein_tri5_all(u - h, v - h);
        for i in 0..TRI5_CP_COUNT {
            let fd = (pp[i] - pm[i] - mp[i] + mm[i]) / (4.0 * h * h);
            assert!(
                (analytic[i] - fd).abs() < 1e-3,
                "idx {i}, (u,v)=({u},{v}): analytic = {}, fd = {fd}",
                analytic[i]
            );
        }
    }
}

#[test]
fn partial_sums_are_zero() {
    // Σ_{i,j,k} ∂B^5_{i,j,k}/∂u = 0 (derivative of the constant 1).
    for &(u, v) in &interior_samples() {
        let du: f64 = bernstein_tri5_partial_u(u, v).iter().sum();
        let dv: f64 = bernstein_tri5_partial_v(u, v).iter().sum();
        assert!(du.abs() < 1e-10, "∂u sum = {du}");
        assert!(dv.abs() < 1e-10, "∂v sum = {dv}");
    }
}

#[test]
fn degree_4_and_3_partition_of_unity() {
    for &(u, v) in &interior_samples() {
        let b4: f64 = bernstein_tri4_all(u, v).iter().sum();
        let b3: f64 = bernstein_tri3_all(u, v).iter().sum();
        assert!((b4 - 1.0).abs() < 1e-12);
        assert!((b3 - 1.0).abs() < 1e-12);
    }
}
