//! Tests for [`super::KpGuidePatchBiQuintic`]. Mirrors
//! `kp_guide_patch_tests.rs` at degree 5 (6 CPs per boundary,
//! quartic cross-derivatives).

#![allow(clippy::needless_range_loop)]

use super::*;

const EPS: f64 = 1e-10;

fn flat_guide() -> KpGuidePatchBiQuintic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            grid[i][j] = Point3::new(i as f64 / 5.0, j as f64 / 5.0, 0.0);
        }
    }
    KpGuidePatchBiQuintic { grid }
}

fn bumpy_guide() -> KpGuidePatchBiQuintic {
    let mut g = flat_guide();
    for i in 1..5 {
        for j in 1..5 {
            let s = i as f64 / 5.0;
            let t = j as f64 / 5.0;
            g.grid[i][j].z = 0.4 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.15 * (s - 0.5) * (t - 0.5);
        }
    }
    g
}

fn eval_quartic_curve(cps: &[Vec3; 6], t: f64) -> Vec3 {
    let omt = 1.0 - t;
    let omt2 = omt * omt;
    let omt3 = omt2 * omt;
    let omt4 = omt3 * omt;
    let omt5 = omt4 * omt;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    let b = [
        omt5,
        5.0 * t * omt4,
        10.0 * t2 * omt3,
        10.0 * t3 * omt2,
        5.0 * t4 * omt,
        t5,
    ];
    let mut out = Vec3::zero();
    for i in 0..6 {
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
    let g = KpGuidePatchBiQuintic::zeros();
    for i in 0..6 {
        for j in 0..6 {
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
    assert!(point_max(g.point_at(1.0, 0.0), g.grid[5][0]) < EPS);
    assert!(point_max(g.point_at(0.0, 1.0), g.grid[0][5]) < EPS);
    assert!(point_max(g.point_at(1.0, 1.0), g.grid[5][5]) < EPS);
}

#[test]
fn boundary_accessors_return_correct_edges() {
    let g = bumpy_guide();
    let b_u0 = g.boundary_u0();
    let b_u1 = g.boundary_u1();
    let b_v0 = g.boundary_v0();
    let b_v1 = g.boundary_v1();
    for j in 0..6 {
        assert_eq!(b_u0[j].x, g.grid[0][j].x);
        assert_eq!(b_u1[j].x, g.grid[5][j].x);
    }
    for i in 0..6 {
        assert_eq!(b_v0[i].x, g.grid[i][0].x);
        assert_eq!(b_v1[i].x, g.grid[i][5].x);
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
        let via_cps = eval_quartic_curve(&cps, v);
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
        let via_cps = eval_quartic_curve(&cps, v);
        let via_partial = g.partial_u(1.0, v);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_deriv_cps_v0_reconstructs_partial_v_at_v_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_v0();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_quartic_curve(&cps, u);
        let via_partial = g.partial_v(u, 0.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_deriv_cps_v1_reconstructs_partial_v_at_v_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_v1();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let via_cps = eval_quartic_curve(&cps, u);
        let via_partial = g.partial_v(u, 1.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_second_deriv_cps_uu_u0_reconstructs_second_partial_uu_at_u_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_uu_u0();
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_quartic_curve(&cps, v);
        let via_partial = g.second_partial_uu(0.0, v);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_second_deriv_cps_uu_u1_reconstructs_second_partial_uu_at_u_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_uu_u1();
    for &v in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_quartic_curve(&cps, v);
        let via_partial = g.second_partial_uu(1.0, v);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_second_deriv_cps_vv_v0_reconstructs_second_partial_vv_at_v_equals_zero() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_vv_v0();
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_quartic_curve(&cps, u);
        let via_partial = g.second_partial_vv(u, 0.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn cross_second_deriv_cps_vv_v1_reconstructs_second_partial_vv_at_v_equals_one() {
    let g = bumpy_guide();
    let cps = g.cross_second_deriv_cps_vv_v1();
    for &u in &[0.0, 0.15, 0.4, 0.6, 0.85, 1.0] {
        let via_cps = eval_quartic_curve(&cps, u);
        let via_partial = g.second_partial_vv(u, 1.0);
        assert!(vec3_max(via_cps, via_partial) < EPS);
    }
}

#[test]
fn endpoint_derivative_cps_agree_at_corners() {
    let g = bumpy_guide();
    let cps = g.cross_deriv_cps_u0();
    let at_0 = g.partial_u(0.0, 0.0);
    let at_1 = g.partial_u(0.0, 1.0);
    assert!(vec3_max(cps[0], at_0) < EPS);
    assert!(vec3_max(cps[5], at_1) < EPS);
}

// ---- Degree-elevation identity ----

#[test]
fn from_bicubic_elevated_is_pointwise_identical() {
    // Build a bicubic guide with curvature and degree-elevate it.
    // At every (u, v) sample the bi-quintic must match the bicubic
    // exactly (pointwise, first + second partials). This is the
    // load-bearing contract for Option C: the optimizer starts
    // from an exact replica of the current bicubic result and only
    // diverges via the extra DoFs.
    let mut bicubic = crate::geometry::gregory::KpGuidePatch::zeros();
    for i in 0..4 {
        for j in 0..4 {
            let s = i as f64 / 3.0;
            let t = j as f64 / 3.0;
            bicubic.grid[i][j] = Point3::new(
                s,
                t,
                0.4 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.15 * (s - 0.5) * (t - 0.5),
            );
        }
    }
    let bq = KpGuidePatchBiQuintic::from_bicubic_elevated(&bicubic);
    let sample_uvs = [
        (0.0, 0.0),
        (0.25, 0.4),
        (0.5, 0.5),
        (0.7, 0.3),
        (0.9, 0.95),
        (1.0, 1.0),
    ];
    for &(u, v) in &sample_uvs {
        assert!(
            point_max(bq.point_at(u, v), bicubic.point_at(u, v)) < 1e-12,
            "point_at mismatch at ({u}, {v})"
        );
        assert!(
            vec3_max(bq.partial_u(u, v), bicubic.partial_u(u, v)) < 1e-12,
            "partial_u mismatch at ({u}, {v})"
        );
        assert!(
            vec3_max(bq.partial_v(u, v), bicubic.partial_v(u, v)) < 1e-12,
            "partial_v mismatch at ({u}, {v})"
        );
        assert!(
            vec3_max(bq.second_partial_uu(u, v), bicubic.second_partial_uu(u, v)) < 1e-11,
            "second_partial_uu at ({u}, {v})"
        );
        assert!(
            vec3_max(bq.second_partial_uv(u, v), bicubic.second_partial_uv(u, v)) < 1e-11,
            "second_partial_uv at ({u}, {v})"
        );
        assert!(
            vec3_max(bq.second_partial_vv(u, v), bicubic.second_partial_vv(u, v)) < 1e-11,
            "second_partial_vv at ({u}, {v})"
        );
    }
}

// ---- Third and fourth partials (corner G2 reconciliation contract) ----

#[test]
fn third_partial_uuv_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-3;
    for &u in &[0.3, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.7] {
            let analytic = g.third_partial_uuv(u, v);
            // ∂³/∂u²∂v via central diff: ∂_v of ∂_uu
            let suu_plus = g.second_partial_uu(u, v + h);
            let suu_minus = g.second_partial_uu(u, v - h);
            let fd = Vec3::new(
                (suu_plus.x - suu_minus.x) / (2.0 * h),
                (suu_plus.y - suu_minus.y) / (2.0 * h),
                (suu_plus.z - suu_minus.z) / (2.0 * h),
            );
            assert!(
                vec3_max(analytic, fd) < 1e-3,
                "∂³uuv at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
        }
    }
}

#[test]
fn third_partial_uvv_matches_finite_difference() {
    let g = bumpy_guide();
    let h = 1e-3;
    for &u in &[0.3, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.7] {
            let analytic = g.third_partial_uvv(u, v);
            let svv_plus = g.second_partial_vv(u + h, v);
            let svv_minus = g.second_partial_vv(u - h, v);
            let fd = Vec3::new(
                (svv_plus.x - svv_minus.x) / (2.0 * h),
                (svv_plus.y - svv_minus.y) / (2.0 * h),
                (svv_plus.z - svv_minus.z) / (2.0 * h),
            );
            assert!(vec3_max(analytic, fd) < 1e-3);
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
            // ∂⁴/∂u²∂v² via central diff on ∂²uv twice
            let suv_pp = g.second_partial_uv(u + h, v + h);
            let suv_pm = g.second_partial_uv(u + h, v - h);
            let suv_mp = g.second_partial_uv(u - h, v + h);
            let suv_mm = g.second_partial_uv(u - h, v - h);
            let fd = Vec3::new(
                (suv_pp.x - suv_pm.x - suv_mp.x + suv_mm.x) / (4.0 * h * h),
                (suv_pp.y - suv_pm.y - suv_mp.y + suv_mm.y) / (4.0 * h * h),
                (suv_pp.z - suv_pm.z - suv_mp.z + suv_mm.z) / (4.0 * h * h),
            );
            assert!(
                vec3_max(analytic, fd) < 1e-2,
                "∂⁴uuvv at ({u}, {v}): analytic {analytic:?}, FD {fd:?}"
            );
        }
    }
}
