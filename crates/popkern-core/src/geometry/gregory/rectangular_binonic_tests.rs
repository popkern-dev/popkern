//! Tests for [`super::RectangularGregoryBiNonic`].
//!
//! Philosophy: verify the tensor-product Bezier evaluation
//! matches finite-difference references for every partial
//! derivative, and that corner / boundary evaluations follow
//! the classical Bezier identities (corner CP = corner value,
//! boundary depends only on boundary CPs).

#![allow(clippy::needless_range_loop)]

use super::*;

/// Flat unit-square patch on the `z = 0` plane with
/// configurable origin + basis. Used as a sanity-check
/// baseline — the surface is exactly the input plane.
fn flat_patch(
    origin: (f64, f64),
    u_dir: (f64, f64),
    v_dir: (f64, f64),
) -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for i in 0..10 {
        for j in 0..10 {
            let s = (i as f64) / 9.0;
            let t = (j as f64) / 9.0;
            grid[i][j] = Point3::new(
                origin.0 + s * u_dir.0 + t * v_dir.0,
                origin.1 + s * u_dir.1 + t * v_dir.1,
                0.0,
            );
        }
    }
    RectangularGregoryBiNonic { grid }
}

/// A non-flat test patch with variation in z — for non-trivial
/// partial-derivative tests. `grid[i][j].z = sin(i/9 · π) · cos(j/9 · π)`
/// gives a smoothly varying surface.
fn bumpy_patch() -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for i in 0..10 {
        for j in 0..10 {
            let s = (i as f64) / 9.0;
            let t = (j as f64) / 9.0;
            grid[i][j] = Point3::new(
                s,
                t,
                (std::f64::consts::PI * s).sin() * (std::f64::consts::PI * t).cos(),
            );
        }
    }
    RectangularGregoryBiNonic { grid }
}

#[test]
fn zeros_has_all_zero_cps() {
    let p = RectangularGregoryBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            assert_eq!(p.grid[i][j].x, 0.0);
            assert_eq!(p.grid[i][j].y, 0.0);
            assert_eq!(p.grid[i][j].z, 0.0);
        }
    }
}

#[test]
fn flat_patch_evaluates_in_plane() {
    let p = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    for &u in &[0.0, 0.2, 0.5, 0.7, 1.0] {
        for &v in &[0.0, 0.1, 0.5, 0.9, 1.0] {
            let pt = p.point_at(u, v);
            assert!(pt.z.abs() < 1e-12, "z = {} at ({u}, {v})", pt.z);
            assert!((pt.x - u).abs() < 1e-10);
            assert!((pt.y - v).abs() < 1e-10);
        }
    }
}

#[test]
fn corner_evaluation_equals_corner_cp() {
    let p = bumpy_patch();
    assert_point_eq(p.point_at(0.0, 0.0), p.grid[0][0]);
    assert_point_eq(p.point_at(1.0, 0.0), p.grid[9][0]);
    assert_point_eq(p.point_at(0.0, 1.0), p.grid[0][9]);
    assert_point_eq(p.point_at(1.0, 1.0), p.grid[9][9]);
}

#[test]
fn boundary_evaluation_only_depends_on_boundary_cps() {
    // On the u=0 boundary, S(0, v) depends only on grid[0][j] for j=0..=9.
    // We perturb a non-boundary CP and confirm boundary values don't change.
    let mut a = bumpy_patch();
    let b = a.clone();
    a.grid[5][5].z += 100.0;
    for k in 0..=9 {
        let v = (k as f64) / 9.0;
        // u = 0 boundary.
        let p_a = a.point_at(0.0, v);
        let p_b = b.point_at(0.0, v);
        assert_point_eq(p_a, p_b);
        // u = 1 boundary.
        let p_a = a.point_at(1.0, v);
        let p_b = b.point_at(1.0, v);
        assert_point_eq(p_a, p_b);
        // v = 0 boundary.
        let u = (k as f64) / 9.0;
        let p_a = a.point_at(u, 0.0);
        let p_b = b.point_at(u, 0.0);
        assert_point_eq(p_a, p_b);
        // v = 1 boundary.
        let p_a = a.point_at(u, 1.0);
        let p_b = b.point_at(u, 1.0);
        assert_point_eq(p_a, p_b);
    }
}

#[test]
fn partial_u_matches_finite_difference() {
    let p = bumpy_patch();
    let h = 1e-5;
    for &u in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        for &v in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let analytic = p.partial_u(u, v);
            let plus = p.point_at(u + h, v);
            let minus = p.point_at(u - h, v);
            let fd = Vec3::new(
                (plus.x - minus.x) / (2.0 * h),
                (plus.y - minus.y) / (2.0 * h),
                (plus.z - minus.z) / (2.0 * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(
                diff < 1e-6,
                "({u}, {v}): analytic {analytic:?} vs FD {fd:?}, diff {diff}"
            );
        }
    }
}

#[test]
fn partial_v_matches_finite_difference() {
    let p = bumpy_patch();
    let h = 1e-5;
    for &u in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        for &v in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let analytic = p.partial_v(u, v);
            let plus = p.point_at(u, v + h);
            let minus = p.point_at(u, v - h);
            let fd = Vec3::new(
                (plus.x - minus.x) / (2.0 * h),
                (plus.y - minus.y) / (2.0 * h),
                (plus.z - minus.z) / (2.0 * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(diff < 1e-6, "({u}, {v}): diff {diff}");
        }
    }
}

#[test]
fn second_partial_uu_matches_finite_difference() {
    let p = bumpy_patch();
    let h = 1e-4;
    for &u in &[0.2, 0.4, 0.6, 0.8] {
        for &v in &[0.2, 0.4, 0.6, 0.8] {
            let analytic = p.second_partial_uu(u, v);
            let plus = p.point_at(u + h, v);
            let center = p.point_at(u, v);
            let minus = p.point_at(u - h, v);
            let fd = Vec3::new(
                (plus.x - 2.0 * center.x + minus.x) / (h * h),
                (plus.y - 2.0 * center.y + minus.y) / (h * h),
                (plus.z - 2.0 * center.z + minus.z) / (h * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(diff < 1e-3, "({u}, {v}): diff {diff}");
        }
    }
}

#[test]
fn second_partial_vv_matches_finite_difference() {
    let p = bumpy_patch();
    let h = 1e-4;
    for &u in &[0.2, 0.4, 0.6, 0.8] {
        for &v in &[0.2, 0.4, 0.6, 0.8] {
            let analytic = p.second_partial_vv(u, v);
            let plus = p.point_at(u, v + h);
            let center = p.point_at(u, v);
            let minus = p.point_at(u, v - h);
            let fd = Vec3::new(
                (plus.x - 2.0 * center.x + minus.x) / (h * h),
                (plus.y - 2.0 * center.y + minus.y) / (h * h),
                (plus.z - 2.0 * center.z + minus.z) / (h * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(diff < 1e-3, "({u}, {v}): diff {diff}");
        }
    }
}

#[test]
fn second_partial_uv_matches_finite_difference() {
    let p = bumpy_patch();
    let h = 1e-4;
    for &u in &[0.2, 0.4, 0.6, 0.8] {
        for &v in &[0.2, 0.4, 0.6, 0.8] {
            let analytic = p.second_partial_uv(u, v);
            // Mixed-partial FD: (f(u+h, v+h) - f(u+h, v-h) - f(u-h, v+h) + f(u-h, v-h)) / (4h²).
            let pp = p.point_at(u + h, v + h);
            let pm = p.point_at(u + h, v - h);
            let mp = p.point_at(u - h, v + h);
            let mm = p.point_at(u - h, v - h);
            let fd = Vec3::new(
                (pp.x - pm.x - mp.x + mm.x) / (4.0 * h * h),
                (pp.y - pm.y - mp.y + mm.y) / (4.0 * h * h),
                (pp.z - pm.z - mp.z + mm.z) / (4.0 * h * h),
            );
            let diff = Vec3::new(analytic.x - fd.x, analytic.y - fd.y, analytic.z - fd.z).length();
            assert!(
                diff < 1e-3,
                "({u}, {v}): analytic {analytic:?}, FD {fd:?}, diff {diff}"
            );
        }
    }
}

#[test]
fn flat_patch_second_partials_are_zero() {
    // A linearly-parameterised flat patch has zero second
    // derivatives in every direction.
    let p = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    for &u in &[0.2, 0.5, 0.7] {
        for &v in &[0.3, 0.5, 0.8] {
            let uu = p.second_partial_uu(u, v);
            let uv = p.second_partial_uv(u, v);
            let vv = p.second_partial_vv(u, v);
            assert!(uu.length() < 1e-10, "uu = {uu:?}");
            assert!(uv.length() < 1e-10, "uv = {uv:?}");
            assert!(vv.length() < 1e-10, "vv = {vv:?}");
        }
    }
}

#[test]
fn endpoint_partials_match_classical_bezier() {
    // At (u=0, v=0): ∂S/∂u(0, 0) = 9·(P[1][0] − P[0][0]).
    // At (u=1, v=0): ∂S/∂u(1, 0) = 9·(P[9][0] − P[8][0]).
    let p = bumpy_patch();
    let du_00_analytic = p.partial_u(0.0, 0.0);
    let du_00_expected = Vec3::new(
        9.0 * (p.grid[1][0].x - p.grid[0][0].x),
        9.0 * (p.grid[1][0].y - p.grid[0][0].y),
        9.0 * (p.grid[1][0].z - p.grid[0][0].z),
    );
    assert!((du_00_analytic.x - du_00_expected.x).abs() < 1e-10);
    assert!((du_00_analytic.y - du_00_expected.y).abs() < 1e-10);
    assert!((du_00_analytic.z - du_00_expected.z).abs() < 1e-10);

    let du_10_analytic = p.partial_u(1.0, 0.0);
    let du_10_expected = Vec3::new(
        9.0 * (p.grid[9][0].x - p.grid[8][0].x),
        9.0 * (p.grid[9][0].y - p.grid[8][0].y),
        9.0 * (p.grid[9][0].z - p.grid[8][0].z),
    );
    assert!((du_10_analytic.x - du_10_expected.x).abs() < 1e-10);
    assert!((du_10_analytic.y - du_10_expected.y).abs() < 1e-10);
    assert!((du_10_analytic.z - du_10_expected.z).abs() < 1e-10);
}

#[test]
fn serde_round_trip_preserves_grid() {
    let p = bumpy_patch();
    let serialized = serde_json::to_string(&p).expect("serialize");
    let deserialized: RectangularGregoryBiNonic =
        serde_json::from_str(&serialized).expect("deserialize");
    for i in 0..10 {
        for j in 0..10 {
            let a = p.grid[i][j];
            let b = deserialized.grid[i][j];
            assert!((a.x - b.x).abs() < 1e-12);
            assert!((a.y - b.y).abs() < 1e-12);
            assert!((a.z - b.z).abs() < 1e-12);
        }
    }
}

#[test]
fn clone_is_independent() {
    let mut a = bumpy_patch();
    let b = a.clone();
    a.grid[5][5].z += 1000.0;
    assert!((a.grid[5][5].z - b.grid[5][5].z).abs() > 100.0);
}

/// Assert two Point3s are equal within a reasonable tolerance.
/// Kept out of the test body for readability.
fn assert_point_eq(a: Point3, b: Point3) {
    assert!(
        (a.x - b.x).abs() < 1e-10 && (a.y - b.y).abs() < 1e-10 && (a.z - b.z).abs() < 1e-10,
        "expected {b:?}, got {a:?}"
    );
}
