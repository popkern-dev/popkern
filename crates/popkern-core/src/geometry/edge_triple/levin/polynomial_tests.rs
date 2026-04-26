//! Sanity tests for the polynomial root solvers used by the Levin
//! pipeline. Tolerances reflect the conservative TOL constants in
//! `polynomial.rs`.

use super::polynomial::{solve_cubic, solve_quadratic, solve_quartic};

const TOL: f64 = 1e-9;

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < TOL
}

#[test]
fn quadratic_two_real_roots() {
    // (x − 2)(x + 3) = x² + x − 6
    let r = solve_quadratic([1.0, 1.0, -6.0]);
    assert_eq!(r.len(), 2);
    assert!(approx_eq(r[0], -3.0));
    assert!(approx_eq(r[1], 2.0));
}

#[test]
fn quadratic_double_root() {
    // (x − 5)² = x² − 10x + 25
    let r = solve_quadratic([1.0, -10.0, 25.0]);
    assert!(!r.is_empty());
    for v in &r {
        assert!(approx_eq(*v, 5.0));
    }
}

#[test]
fn quadratic_no_real_roots() {
    // x² + 1
    let r = solve_quadratic([1.0, 0.0, 1.0]);
    assert!(r.is_empty());
}

#[test]
fn cubic_three_distinct_real_roots() {
    // (x − 1)(x − 2)(x − 3) = x³ − 6x² + 11x − 6
    let r = solve_cubic([1.0, -6.0, 11.0, -6.0]);
    assert_eq!(r.len(), 3);
    assert!(approx_eq(r[0], 1.0));
    assert!(approx_eq(r[1], 2.0));
    assert!(approx_eq(r[2], 3.0));
}

#[test]
fn cubic_one_real_two_complex() {
    // x³ + x + 1 — one real root near −0.6824
    let r = solve_cubic([1.0, 0.0, 1.0, 1.0]);
    assert_eq!(r.len(), 1);
    let val = r[0].powi(3) + r[0] + 1.0;
    assert!(val.abs() < TOL, "p({}) = {}", r[0], val);
}

#[test]
fn cubic_with_repeated_root() {
    // (x − 1)²(x − 2) = x³ − 4x² + 5x − 2
    let r = solve_cubic([1.0, -4.0, 5.0, -2.0]);
    // Should report at least both distinct values (1 and 2), with 1 possibly twice.
    let has_one = r.iter().any(|&v| approx_eq(v, 1.0));
    let has_two = r.iter().any(|&v| approx_eq(v, 2.0));
    assert!(has_one && has_two, "got {r:?}");
}

#[test]
fn quartic_four_distinct_real_roots() {
    // (x − 1)(x − 2)(x − 3)(x − 4) = x⁴ − 10x³ + 35x² − 50x + 24
    let r = solve_quartic([1.0, -10.0, 35.0, -50.0, 24.0]);
    assert_eq!(r.len(), 4);
    for (i, target) in [1.0, 2.0, 3.0, 4.0].iter().enumerate() {
        assert!(approx_eq(r[i], *target), "root {i}: {} vs {target}", r[i]);
    }
}

#[test]
fn quartic_two_real_two_complex() {
    // (x² + 1)(x − 1)(x + 2) = x⁴ + x³ − x² + x − 2
    let r = solve_quartic([1.0, 1.0, -1.0, 1.0, -2.0]);
    assert_eq!(r.len(), 2);
    assert!(approx_eq(r[0], -2.0));
    assert!(approx_eq(r[1], 1.0));
}

#[test]
fn quartic_no_real_roots() {
    // (x² + 1)(x² + 2) = x⁴ + 3x² + 2
    let r = solve_quartic([1.0, 0.0, 3.0, 0.0, 2.0]);
    assert!(r.is_empty());
}

#[test]
fn quartic_biquadratic() {
    // x⁴ − 5x² + 4 = (x² − 1)(x² − 4) → roots ±1, ±2
    let r = solve_quartic([1.0, 0.0, -5.0, 0.0, 4.0]);
    assert_eq!(r.len(), 4);
    assert!(approx_eq(r[0], -2.0));
    assert!(approx_eq(r[1], -1.0));
    assert!(approx_eq(r[2], 1.0));
    assert!(approx_eq(r[3], 2.0));
}

#[test]
fn quartic_drops_to_cubic() {
    // 0·x⁴ + (x − 1)(x − 2)(x − 3)
    let r = solve_quartic([0.0, 1.0, -6.0, 11.0, -6.0]);
    assert_eq!(r.len(), 3);
    assert!(approx_eq(r[0], 1.0));
    assert!(approx_eq(r[1], 2.0));
    assert!(approx_eq(r[2], 3.0));
}
