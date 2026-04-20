//! Tests for [`super::KpGuidePatch`].
//!
//! Two layers of verification:
//!
//! 1. **Point / partial evaluation** — the bicubic must agree with
//!    a finite-difference reference on every first + second
//!    partial, and match classical Bezier identities at corners
//!    (corner CP = corner value, boundary depends only on boundary
//!    CPs, etc.).
//! 2. **Derivative-CP identities** — the `cross_deriv_cps_*`
//!    and `cross_second_deriv_cps_*` helpers return 4 Bernstein
//!    CPs of a cubic curve. Those CPs must reconstruct the
//!    corresponding boundary derivative exactly — this is the
//!    load-bearing contract for the H.4/H.5 row-1/row-2 formulas
//!    (the K-P construction uses the CPs as input, NOT point
//!    samples of the derivative curve).

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-10;

/// Flat unit-square guide on `z = 0`. Linear in both `u` and `v`,
/// so `point_at(u, v) = (u, v, 0)` exactly.
fn flat_guide() -> KpGuidePatch {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            grid[i][j] = Point3::new(i as f64 / 3.0, j as f64 / 3.0, 0.0);
        }
    }
    KpGuidePatch { grid }
}

/// Non-trivial guide with curvature. Base grid is the unit square
/// on `z = 0`; interior CPs get a smooth z-bump so every first
/// and second partial is nonzero somewhere.
fn bumpy_guide() -> KpGuidePatch {
    let mut g = flat_guide();
    for i in 1..3 {
        for j in 1..3 {
            let s = i as f64 / 3.0;
            let t = j as f64 / 3.0;
            g.grid[i][j].z = 0.4 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.15 * (s - 0.5) * (t - 0.5);
        }
    }
    g
}

/// Evaluate a cubic Bezier curve (4 CPs) at parameter `t`.
fn eval_cubic_curve(cps: &[Vec3; 4], t: f64) -> Vec3 {
    let omt = 1.0 - t;
    let b = [
        omt * omt * omt,
        3.0 * t * omt * omt,
        3.0 * t * t * omt,
        t * t * t,
    ];
    let mut out = Vec3::zero();
    for i in 0..4 {
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
    let g = KpGuidePatch::zeros();
    for i in 0..4 {
        for j in 0..4 {
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
            assert!(pt.z.abs() < EPS, "z = {} at ({u}, {v})", pt.z);
            assert!((pt.x - u).abs() < EPS);
            assert!((pt.y - v).abs() < EPS);
        }
    }
}

#[test]
fn corner_evaluation_equals_corner_cp() {
    let g = bumpy_guide();
    assert!(point_max(g.point_at(0.0, 0.0), g.grid[0][0]) < EPS);
    assert!(point_max(g.point_at(1.0, 0.0), g.grid[3][0]) < EPS);
    assert!(point_max(g.point_at(0.0, 1.0), g.grid[0][3]) < EPS);
    assert!(point_max(g.point_at(1.0, 1.0), g.grid[3][3]) < EPS);
}

#[test]
fn boundary_accessors_return_correct_edges() {
    let g = bumpy_guide();
    let b_u0 = g.boundary_u0();
    let b_u1 = g.boundary_u1();
    let b_v0 = g.boundary_v0();
    let b_v1 = g.boundary_v1();
    for j in 0..4 {
        assert_eq!(b_u0[j].x, g.grid[0][j].x);
        assert_eq!(b_u1[j].x, g.grid[3][j].x);
    }
    for i in 0..4 {
        assert_eq!(b_v0[i].x, g.grid[i][0].x);
        assert_eq!(b_v1[i].x, g.grid[i][3].x);
    }
}

// ---- Partial derivatives match finite difference ----

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
            assert!(
                vec3_max(analytic, fd) < 1e-6,
                "∂u at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
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
            assert!(
                vec3_max(analytic, fd) < 1e-6,
                "∂v at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
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
            assert!(
                vec3_max(analytic, fd) < 1e-3,
                "∂²uu at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
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
            assert!(
                vec3_max(analytic, fd) < 1e-3,
                "∂²vv at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
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
            assert!(
                vec3_max(analytic, fd) < 1e-3,
                "∂²uv at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
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

// ---- Derivative-CP identities (the load-bearing H.4 contract) ----

#[test]
fn cross_deriv_cps_u0_reconstructs_partial_u_at_u_equals_zero() {
    // The 4 CPs returned by `cross_deriv_cps_u0` are the Bernstein
    // CPs of the cubic curve ∂G/∂u(0, v) in v. Evaluating that
    // cubic at any v must equal `partial_u(0, v)`.
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_u0();
    for &v in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_cubic_curve(&cps, v);
        let via_partial = g.partial_u(0.0, v);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂u(0, {v}): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_deriv_cps_u1_reconstructs_partial_u_at_u_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_u1();
    for &v in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_cubic_curve(&cps, v);
        let via_partial = g.partial_u(1.0, v);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂u(1, {v}): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_deriv_cps_v0_reconstructs_partial_v_at_v_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_v0();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_cubic_curve(&cps, u);
        let via_partial = g.partial_v(u, 0.0);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂v({u}, 0): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_deriv_cps_v1_reconstructs_partial_v_at_v_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_v1();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_cubic_curve(&cps, u);
        let via_partial = g.partial_v(u, 1.0);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂v({u}, 1): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_second_deriv_cps_uu_u0_reconstructs_second_partial_uu_at_u_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_uu_u0();
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_cubic_curve(&cps, v);
        let via_partial = g.second_partial_uu(0.0, v);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂²uu(0, {v}): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_second_deriv_cps_uu_u1_reconstructs_second_partial_uu_at_u_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_uu_u1();
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_cubic_curve(&cps, v);
        let via_partial = g.second_partial_uu(1.0, v);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂²uu(1, {v}): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_second_deriv_cps_vv_v0_reconstructs_second_partial_vv_at_v_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_vv_v0();
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_cubic_curve(&cps, u);
        let via_partial = g.second_partial_vv(u, 0.0);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂²vv({u}, 0): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn cross_second_deriv_cps_vv_v1_reconstructs_second_partial_vv_at_v_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_vv_v1();
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_cubic_curve(&cps, u);
        let via_partial = g.second_partial_vv(u, 1.0);
        assert!(
            vec3_max(via_cps, via_partial) < EPS,
            "∂²vv({u}, 1): via CPs {via_cps:?}, via partial {via_partial:?}"
        );
    }
}

#[test]
fn endpoint_derivative_cps_agree_at_corners() {
    // At v = 0: cps[0] of cross_deriv_cps_u0 equals partial_u(0, 0).
    // At v = 1: cps[3] equals partial_u(0, 1).
    // (Classical Bezier-curve-endpoint identity.)
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_u0();
    let at_0 = g.partial_u(0.0, 0.0);
    let at_1 = g.partial_u(0.0, 1.0);
    assert!(vec3_max(cps[0], at_0) < EPS, "cps[0] vs ∂u(0, 0)");
    assert!(vec3_max(cps[3], at_1) < EPS, "cps[3] vs ∂u(0, 1)");
}
