//! Tests for [`pns_f_kappa`] — F_κ quality functional on bi-4 patches.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;

#[test]
fn constant_patch_has_zero_energy() {
    // A constant patch has all partial derivatives zero → F_κ = 0.
    let p = Point3::new(1.0, 2.0, 3.0);
    let grid = [[p; 5]; 5];
    for kappa in 1..=4 {
        let e = f_kappa_bi4(&grid, kappa);
        assert!(e.abs() < 1e-20, "F_{kappa} on constant patch = {e}");
    }
}

#[test]
fn linear_patch_has_zero_f3_and_f4() {
    // A linear patch (affine in s, t) has all ≥ 2nd derivatives
    // zero → F_2 = F_3 = F_4 = 0.
    //
    // Set grid[i][j] = A + B·(i/4) + C·(j/4) for some A, B, C.
    let a = Point3::new(1.0, 0.0, 0.0);
    let b = Point3::new(0.5, 1.0, 0.0);
    let c = Point3::new(0.0, 0.3, 1.0);
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            let u = (i as f64) / 4.0;
            let v = (j as f64) / 4.0;
            grid[i][j] = Point3::new(
                a.x + b.x * u + c.x * v,
                a.y + b.y * u + c.y * v,
                a.z + b.z * u + c.z * v,
            );
        }
    }
    // F_2, F_3, F_4 should all be zero (≥2nd partial derivs vanish).
    for kappa in 2..=4 {
        let e = f_kappa_bi4(&grid, kappa);
        assert!(
            e.abs() < 1e-12,
            "F_{kappa} on linear patch should be 0, got {e}"
        );
    }
}

#[test]
fn f_kappa_is_non_negative() {
    // F_κ is a sum of squares integrated with positive weights.
    // Must be ≥ 0 for any patch.
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            grid[i][j] = Point3::new(
                (i as f64 - 2.0) * (j as f64 - 2.0),
                (i as f64) * 0.3 + (j as f64) * 0.2,
                (i as f64).powi(2) - (j as f64),
            );
        }
    }
    for kappa in 1..=4 {
        let e = f_kappa_bi4(&grid, kappa);
        assert!(e >= -1e-12, "F_{kappa} must be ≥ 0, got {e}");
    }
}

#[test]
fn f_kappa_is_translation_invariant() {
    // Translating all CPs by the same Point3 doesn't change any
    // partial derivative → F_κ is invariant.
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            grid[i][j] = Point3::new(
                0.1 * (i as f64) * (j as f64),
                0.2 * (i as f64).powi(2),
                0.3 * (j as f64).powi(2),
            );
        }
    }
    let delta = Point3::new(100.0, -50.0, 25.5);
    let mut shifted = grid;
    for i in 0..5 {
        for j in 0..5 {
            shifted[i][j] = Point3::new(
                grid[i][j].x + delta.x,
                grid[i][j].y + delta.y,
                grid[i][j].z + delta.z,
            );
        }
    }
    for kappa in 1..=4 {
        let orig = f_kappa_bi4(&grid, kappa);
        let shifted_val = f_kappa_bi4(&shifted, kappa);
        assert!(
            (orig - shifted_val).abs() < 1e-8,
            "F_{kappa}: orig={orig}, shifted={shifted_val}"
        );
    }
}

#[test]
fn f_kappa_is_quadratic_in_cps() {
    // Scaling all CPs by α scales (∂p)² by α² and thus F_κ by α².
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    for i in 0..5 {
        for j in 0..5 {
            grid[i][j] = Point3::new(
                ((i as f64) * 0.3).sin(),
                ((j as f64) * 0.5).cos(),
                ((i + j) as f64) * 0.1,
            );
        }
    }
    let alpha = 2.5_f64;
    let mut scaled = grid;
    for i in 0..5 {
        for j in 0..5 {
            scaled[i][j] = Point3::new(
                alpha * grid[i][j].x,
                alpha * grid[i][j].y,
                alpha * grid[i][j].z,
            );
        }
    }
    for kappa in 1..=4 {
        let orig = f_kappa_bi4(&grid, kappa);
        let scaled_val = f_kappa_bi4(&scaled, kappa);
        assert!(
            (scaled_val - alpha * alpha * orig).abs() < 1e-8,
            "F_{kappa} quadratic-scaling: orig={orig}, scaled={scaled_val}, alpha²·orig={}",
            alpha * alpha * orig
        );
    }
}

#[test]
fn partial_at_0_0_equals_bb_corner_formula() {
    // ∂^i_s ∂^j_t p(0, 0) = (4!/(4−i)!)·(4!/(4−j)!)·Δ^i_s Δ^j_t p_{0, 0}
    // For bi-4 with i=j=1: 4·4 = 16, acting on first mixed difference
    // (p_{11} − p_{01} − p_{10} + p_{00}).
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 5]; 5];
    grid[0][0] = Point3::new(1.0, 0.0, 0.0);
    grid[1][0] = Point3::new(2.0, 0.0, 0.0);
    grid[0][1] = Point3::new(0.0, 1.0, 0.0);
    grid[1][1] = Point3::new(3.0, 1.5, 0.0);
    // ∂_s ∂_t p(0, 0) = 16 · (p_{11} − p_{01} − p_{10} + p_{00})
    //                  = 16 · ((3, 1.5, 0) − (0, 1, 0) − (2, 0, 0) + (1, 0, 0))
    //                  = 16 · (2, 0.5, 0)
    //                  = (32, 8, 0)
    let d = super::partial_bi4(&grid, 1, 1, 0.0, 0.0);
    assert!((d.x - 32.0).abs() < 1e-10, "∂_s∂_t p_x = {}", d.x);
    assert!((d.y - 8.0).abs() < 1e-10, "∂_s∂_t p_y = {}", d.y);
    assert!(d.z.abs() < 1e-10, "∂_s∂_t p_z = {}", d.z);
}
