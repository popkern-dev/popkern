//! Tests for [`super::KpGuidePatchBiNonic`]. Mirrors the
//! bi-quintic test suite at degree 9 — FD parity on all
//! partials, derivative-CP identities, elevation identities.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::{Point3, Vec3};

const EPS: f64 = 1e-10;

fn flat_guide() -> KpGuidePatchBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for i in 0..10 {
        for j in 0..10 {
            grid[i][j] = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.0);
        }
    }
    KpGuidePatchBiNonic { grid }
}

fn bumpy_guide() -> KpGuidePatchBiNonic {
    let mut g = flat_guide();
    // z-bump on the interior. Asymmetric so all partials are nonzero.
    for i in 1..9 {
        for j in 1..9 {
            let s = i as f64 / 9.0;
            let t = j as f64 / 9.0;
            g.grid[i][j].z = 0.4 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.15 * (s - 0.5) * (t - 0.5);
        }
    }
    g
}

fn eval_binonic_curve(cps: &[Vec3; 10], t: f64) -> Vec3 {
    use crate::geometry::gregory::bernstein_binonic::bernstein_binonic_all;
    let b = bernstein_binonic_all(t);
    let mut out = Vec3::zero();
    for i in 0..10 {
        out.x += cps[i].x * b[i];
        out.y += cps[i].y * b[i];
        out.z += cps[i].z * b[i];
    }
    out
}

fn vec3_max(a: Vec3, b: Vec3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

// ---- Basic evaluation ----

#[test]
fn zeros_has_all_zero_cps() {
    let g = KpGuidePatchBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            assert_eq!(g.grid[i][j].x, 0.0);
            assert_eq!(g.grid[i][j].y, 0.0);
            assert_eq!(g.grid[i][j].z, 0.0);
        }
    }
}

#[test]
fn flat_guide_evaluates_in_plane() {
    let g = flat_guide();
    for &u in &[0.0, 0.2, 0.5, 0.7, 1.0] {
        for &v in &[0.0, 0.1, 0.5, 0.9, 1.0] {
            let pt = g.point_at(u, v);
            assert!(pt.z.abs() < EPS);
            assert!((pt.x - u).abs() < EPS);
            assert!((pt.y - v).abs() < EPS);
        }
    }
}

#[test]
fn corner_evaluation_equals_corner_cp() {
    let g = bumpy_guide();
    assert!(point_max(g.point_at(0.0, 0.0), g.grid[0][0]) < EPS);
    assert!(point_max(g.point_at(1.0, 0.0), g.grid[9][0]) < EPS);
    assert!(point_max(g.point_at(0.0, 1.0), g.grid[0][9]) < EPS);
    assert!(point_max(g.point_at(1.0, 1.0), g.grid[9][9]) < EPS);
}

// ---- Partial derivatives vs finite difference ----

#[test]
fn partial_u_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-5;
    for &u in &[0.15, 0.4, 0.6, 0.85] {
        for &v in &[0.15, 0.4, 0.6, 0.85] {
            let analytic = g.partial_u(u, v);
            let plus = g.point_at(u + h, v);
            let minus = g.point_at(u - h, v);
            let fd = Vec3::new(
                (plus.x - minus.x) / (2.0 * h),
                (plus.y - minus.y) / (2.0 * h),
                (plus.z - minus.z) / (2.0 * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-6);
        }
    }
}

#[test]
fn partial_v_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-5;
    for &u in &[0.15, 0.4, 0.6, 0.85] {
        for &v in &[0.15, 0.4, 0.6, 0.85] {
            let analytic = g.partial_v(u, v);
            let plus = g.point_at(u, v + h);
            let minus = g.point_at(u, v - h);
            let fd = Vec3::new(
                (plus.x - minus.x) / (2.0 * h),
                (plus.y - minus.y) / (2.0 * h),
                (plus.z - minus.z) / (2.0 * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-6);
        }
    }
}

#[test]
fn second_partial_uu_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-4;
    for &u in &[0.2, 0.5, 0.8] {
        for &v in &[0.2, 0.5, 0.8] {
            let analytic = g.second_partial_uu(u, v);
            let plus = g.point_at(u + h, v);
            let center = g.point_at(u, v);
            let minus = g.point_at(u - h, v);
            let fd = Vec3::new(
                (plus.x - 2.0 * center.x + minus.x) / (h * h),
                (plus.y - 2.0 * center.y + minus.y) / (h * h),
                (plus.z - 2.0 * center.z + minus.z) / (h * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-3);
        }
    }
}

#[test]
fn second_partial_uv_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-4;
    for &u in &[0.2, 0.5, 0.8] {
        for &v in &[0.2, 0.5, 0.8] {
            let analytic = g.second_partial_uv(u, v);
            let pp = g.point_at(u + h, v + h);
            let pm = g.point_at(u + h, v - h);
            let mp = g.point_at(u - h, v + h);
            let mm = g.point_at(u - h, v - h);
            let fd = Vec3::new(
                (pp.x - pm.x - mp.x + mm.x) / (4.0 * h * h),
                (pp.y - pm.y - mp.y + mm.y) / (4.0 * h * h),
                (pp.z - pm.z - mp.z + mm.z) / (4.0 * h * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-3);
        }
    }
}

#[test]
fn second_partial_vv_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-4;
    for &u in &[0.2, 0.5, 0.8] {
        for &v in &[0.2, 0.5, 0.8] {
            let analytic = g.second_partial_vv(u, v);
            let plus = g.point_at(u, v + h);
            let center = g.point_at(u, v);
            let minus = g.point_at(u, v - h);
            let fd = Vec3::new(
                (plus.x - 2.0 * center.x + minus.x) / (h * h),
                (plus.y - 2.0 * center.y + minus.y) / (h * h),
                (plus.z - 2.0 * center.z + minus.z) / (h * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-3);
        }
    }
}

#[test]
fn flat_guide_has_zero_second_partials() {
    let g = flat_guide();
    for &u in &[0.2, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.8] {
            assert!(g.second_partial_uu(u, v).length() < EPS);
            assert!(g.second_partial_uv(u, v).length() < EPS);
            assert!(g.second_partial_vv(u, v).length() < EPS);
        }
    }
}

// ---- Derivative-CP identities ----

#[test]
fn cross_deriv_cps_u0_reconstructs_partial_u_at_u_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_u0();
    for &v in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_binonic_curve(&cps, v);
        let via_partial = g.partial_u(0.0, v);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_deriv_cps_u1_reconstructs_partial_u_at_u_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_u1();
    for &v in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_binonic_curve(&cps, v);
        let via_partial = g.partial_u(1.0, v);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_deriv_cps_v0_reconstructs_partial_v_at_v_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_v0();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_binonic_curve(&cps, u);
        let via_partial = g.partial_v(u, 0.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_deriv_cps_v1_reconstructs_partial_v_at_v_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_v1();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_binonic_curve(&cps, u);
        let via_partial = g.partial_v(u, 1.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_second_deriv_cps_uu_u0_reconstructs_second_partial_uu_at_u_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_uu_u0();
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_binonic_curve(&cps, v);
        let via_partial = g.second_partial_uu(0.0, v);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_second_deriv_cps_vv_v1_reconstructs_second_partial_vv_at_v_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_vv_v1();
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_binonic_curve(&cps, u);
        let via_partial = g.second_partial_vv(u, 1.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

// ---- Degree elevation identities ----

#[test]
fn from_bicubic_elevated_is_pointwise_identical() {
    let mut bi = crate::geometry::gregory::KpGuidePatch::zeros();
    for i in 0..4 {
        for j in 0..4 {
            let s = i as f64 / 3.0;
            let t = j as f64 / 3.0;
            bi.grid[i][j] = Point3::new(
                s,
                t,
                0.3 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.2 * (s - 0.5) * (t - 0.5),
            );
        }
    }
    let bn = KpGuidePatchBiNonic::from_bicubic_elevated(&bi);
    for &(u, v) in &[(0.0, 0.0), (0.3, 0.4), (0.5, 0.5), (0.8, 0.7), (1.0, 1.0)] {
        assert!(point_max(bn.point_at(u, v), bi.point_at(u, v)) < 1e-11);
        assert!(vec3_max(bn.partial_u(u, v), bi.partial_u(u, v)) < 1e-11);
        assert!(vec3_max(bn.partial_v(u, v), bi.partial_v(u, v)) < 1e-11);
        assert!(vec3_max(bn.second_partial_uu(u, v), bi.second_partial_uu(u, v)) < 1e-10);
        assert!(vec3_max(bn.second_partial_uv(u, v), bi.second_partial_uv(u, v)) < 1e-10);
        assert!(vec3_max(bn.second_partial_vv(u, v), bi.second_partial_vv(u, v)) < 1e-10);
    }
}

#[test]
fn from_biquintic_elevated_is_pointwise_identical() {
    use crate::geometry::gregory::{KpGuidePatch, KpGuidePatchBiQuintic};
    // Build a bi-quintic via its from_bicubic_elevated.
    let mut bi = KpGuidePatch::zeros();
    for i in 0..4 {
        for j in 0..4 {
            let s = i as f64 / 3.0;
            let t = j as f64 / 3.0;
            bi.grid[i][j] = Point3::new(s, t, 0.35 * (s - 0.5) * (t - 0.3));
        }
    }
    let bq = KpGuidePatchBiQuintic::from_bicubic_elevated(&bi);
    let bn = KpGuidePatchBiNonic::from_biquintic_elevated(&bq);
    for &(u, v) in &[(0.0, 0.0), (0.3, 0.4), (0.5, 0.5), (0.8, 0.7), (1.0, 1.0)] {
        assert!(point_max(bn.point_at(u, v), bq.point_at(u, v)) < 1e-11);
        assert!(vec3_max(bn.partial_u(u, v), bq.partial_u(u, v)) < 1e-11);
        assert!(vec3_max(bn.partial_v(u, v), bq.partial_v(u, v)) < 1e-11);
    }
}

// ---- Higher partials (for G2 corner reconciliation) ----

#[test]
fn third_partial_uuv_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-3;
    for &u in &[0.3, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.7] {
            let analytic = g.third_partial_uuv(u, v);
            let suu_plus = g.second_partial_uu(u, v + h);
            let suu_minus = g.second_partial_uu(u, v - h);
            let fd = Vec3::new(
                (suu_plus.x - suu_minus.x) / (2.0 * h),
                (suu_plus.y - suu_minus.y) / (2.0 * h),
                (suu_plus.z - suu_minus.z) / (2.0 * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-2);
        }
    }
}

#[test]
fn fourth_partial_uuvv_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 2e-3;
    for &u in &[0.3, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.7] {
            let analytic = g.fourth_partial_uuvv(u, v);
            let suv_pp = g.second_partial_uv(u + h, v + h);
            let suv_pm = g.second_partial_uv(u + h, v - h);
            let suv_mp = g.second_partial_uv(u - h, v + h);
            let suv_mm = g.second_partial_uv(u - h, v - h);
            let fd = Vec3::new(
                (suv_pp.x - suv_pm.x - suv_mp.x + suv_mm.x) / (4.0 * h * h),
                (suv_pp.y - suv_pm.y - suv_mp.y + suv_mm.y) / (4.0 * h * h),
                (suv_pp.z - suv_pm.z - suv_mp.z + suv_mm.z) / (4.0 * h * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-1);
        }
    }
}
