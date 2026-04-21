//! Tests for the KP09a C¹/C² inter-sector join operators.
//!
//! # Strategy
//!
//! 1. **Weights**: `ν₀ + ν₁ + ν₂ = 1` (affine property) — verified
//!    across the valence range.
//! 2. **Constant-patch invariance**: applying `c1_join` / `c2_join`
//!    to a sector whose CPs are all a single point `P` returns `P`
//!    (follows from affine property).
//! 3. **Global-polynomial match**: build a globally-defined
//!    bivariate polynomial `P(x, y)` and sample its BB-coefs on two
//!    adjacent wedges. Applying KP09a eq (3)/(4) to sector `i`'s
//!    coefs should reproduce sector `i+1`'s observed depth-1 /
//!    depth-2 CPs.

use super::super::kp13_bi5_gtri5::Gtri5Sector;
use super::*;

const EPS: f64 = 1e-10;

fn approx_point(a: Point3, b: Point3) -> bool {
    (a.x - b.x).abs() < EPS && (a.y - b.y).abs() < EPS && (a.z - b.z).abs() < EPS
}

// ---- Weights ----

#[test]
fn weights_sum_to_one() {
    for n in 3..=20 {
        let w = ConnectWeights::for_valence(n);
        assert!((w.sum() - 1.0).abs() < EPS, "n = {n}: sum {} != 1", w.sum());
    }
}

#[test]
fn weights_nu2_is_minus_one() {
    for n in 3..=15 {
        let w = ConnectWeights::for_valence(n);
        assert_eq!(w.nu2, -1.0);
    }
}

#[test]
fn weights_n6_match_regular() {
    // For n = 6 (regular): c = cos(π/3) = 0.5, c̄ = 0.5.
    // ν₀ = 2·0.5 = 1,  ν₁ = 2·0.5 = 1,  ν₂ = −1. Sum = 1.
    let w = ConnectWeights::for_valence(6);
    assert!((w.nu0 - 1.0).abs() < EPS);
    assert!((w.nu1 - 1.0).abs() < EPS);
    assert_eq!(w.nu2, -1.0);
}

// ---- Constant-patch invariance ----

#[test]
fn c1_join_constant_patch() {
    let p = Point3::new(1.5, -2.0, 3.25);
    let sector = Gtri5Sector { cps: [p; 21] };
    let w = ConnectWeights::for_valence(5);
    for j in 0..=4 {
        let out = c1_join(&w, &sector, j);
        assert!(approx_point(out, p), "c1_join at j = {j}: got {out:?}");
    }
}

#[test]
fn c2_join_constant_patch() {
    let p = Point3::new(-0.3, 0.7, 1.1);
    let sector = Gtri5Sector { cps: [p; 21] };
    let w = ConnectWeights::for_valence(7);
    for j in 0..=3 {
        let out = c2_join(&w, &sector, j);
        assert!(approx_point(out, p), "c2_join at j = {j}: got {out:?}");
    }
}

// ---- Global-polynomial match ----

/// Sample the scalar field `f(x, y) = a*x + b*y + c` over sector
/// `k`'s wedge, producing a Gtri5Sector (z-coord = 0).
fn sample_linear(n: usize, k: usize, a: f64, b: f64, c: f64) -> Gtri5Sector {
    let (xi, yi) = super::super::kp13_bi5_gtri5::sector_corner(n, k);
    let (xi1, yi1) = super::super::kp13_bi5_gtri5::sector_corner(n, k + 1);
    let mut s = Gtri5Sector::zero();
    // Iterate (i, j) with i + j ≤ 5 → CP at barycentric (i/5, j/5, (5-i-j)/5)
    // → geometric at i/5·v_k + j/5·v_{k+1} (star at origin).
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let u = i as f64 / 5.0;
            let v = j as f64 / 5.0;
            let x = u * xi + v * xi1;
            let y = u * yi + v * yi1;
            s.set_cp(i, j, Point3::new(a * x + b * y + c, 0.0, 0.0));
        }
    }
    s
}

#[test]
fn c1_join_reproduces_linear_across_boundary() {
    // For a linear (degree-1) polynomial, every BB-coef is just the
    // polynomial value at that CP's geometric position — so the
    // sampled sectors should be automatically C² connected (they
    // come from a single global polynomial). Therefore applying
    // eq (3) to sector i's CPs should give exactly sector i+1's
    // paper-(j, 1) CPs.
    let (a, b, c) = (0.3, -0.7, 1.2);
    let n = 6;
    let w = ConnectWeights::for_valence(n);

    for i in 0..n {
        let s_i = sample_linear(n, i, a, b, c);
        let s_ip1 = sample_linear(n, i + 1, a, b, c);
        for j in 0..=4 {
            let predicted = c1_join(&w, &s_i, j);
            let actual = s_ip1.cp(j, 1);
            assert!(
                (predicted.x - actual.x).abs() < 1e-10,
                "n = {n}, i = {i}, j = {j}: predicted = {predicted:?}, actual = {actual:?}"
            );
        }
    }
}

#[test]
fn c2_join_reproduces_linear_across_boundary() {
    // Same as above for C² — a linear polynomial has trivially
    // vanishing higher derivatives, so every layer matches.
    let (a, b, c) = (-0.5, 0.4, 0.1);
    let n = 5;
    let w = ConnectWeights::for_valence(n);

    for i in 0..n {
        let s_i = sample_linear(n, i, a, b, c);
        let s_ip1 = sample_linear(n, i + 1, a, b, c);
        for j in 0..=3 {
            let predicted = c2_join(&w, &s_i, j);
            let actual = s_ip1.cp(j, 2);
            assert!(
                (predicted.x - actual.x).abs() < 1e-10,
                "n = {n}, i = {i}, j = {j}: predicted = {predicted:?}, actual = {actual:?}"
            );
        }
    }
}

// ---- paper_flat ----

#[test]
fn paper_flat_star_corner_ids() {
    // Paper (0, 0) — star — maps to flat 0.
    assert_eq!(paper_flat(0, 0), 0);
    // Paper (5, 0) — v_i — flat 20.
    assert_eq!(paper_flat(5, 0), 20);
    // Paper (0, 5) — v_{i+1} — flat 5.
    assert_eq!(paper_flat(0, 5), 5);
}

#[test]
#[should_panic]
fn c1_join_panics_on_oob_j() {
    let sector = Gtri5Sector::zero();
    let w = ConnectWeights::for_valence(6);
    let _ = c1_join(&w, &sector, 5); // j + 1 > 5
}

#[test]
#[should_panic]
fn c2_join_panics_on_oob_j() {
    let sector = Gtri5Sector::zero();
    let w = ConnectWeights::for_valence(6);
    let _ = c2_join(&w, &sector, 4); // j + 2 > 5
}
