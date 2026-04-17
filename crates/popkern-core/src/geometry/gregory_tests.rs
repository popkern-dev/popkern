//! Tests for the rational Gregory patch evaluator.
//!
//! Focus: correctness of boundary evaluation, derivative
//! agreement with finite differences, and the twist-compatibility
//! property at corners. The fillet-side construction (extracting
//! `inner_cps` from incident blend strips) is exercised in the
//! Phase 3 commits; these tests use canned CPs representing a
//! mildly-warped bicubic patch so we can independently verify the
//! maths.

use super::*;

/// A flat planar rectangular Gregory patch — all 16 boundary CPs
/// lie on the z=0 plane on a regular grid, inner CPs are chosen
/// to match a standard bicubic Bezier (both candidates equal).
/// Evaluating must produce points on z=0 regardless of (u, v).
fn flat_rect_patch() -> RectangularGregory {
    let mut boundary = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in boundary.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = Point3::new(i as f64 / 3.0, j as f64 / 3.0, 0.0);
        }
    }
    // Inner CPs: both candidates equal the ordinary bicubic
    // interior CP — degenerate twist-compat case.
    let inner = (
        Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.0),
        Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.0),
    );
    let inner12 = (
        Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
        Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
    );
    let inner22 = (
        Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.0),
        Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.0),
    );
    let inner21 = (
        Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
        Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
    );
    RectangularGregory {
        boundary_cps: boundary,
        inner_cps: [inner, inner12, inner22, inner21],
    }
}

#[test]
fn flat_rect_patch_evaluates_on_plane() {
    let patch = flat_rect_patch();
    for u in [0.0, 0.2, 0.5, 0.8, 1.0] {
        for v in [0.0, 0.2, 0.5, 0.8, 1.0] {
            let p = patch.point_at(u, v);
            assert!(p.z.abs() < 1e-12, "z should be ~0, got {}", p.z);
            // Bilinear on a flat patch: x = u, y = v (scaled 0..1).
            assert!((p.x - u).abs() < 1e-10, "x = {} expected {}", p.x, u);
            assert!((p.y - v).abs() < 1e-10, "y = {} expected {}", p.y, v);
        }
    }
}

#[test]
fn flat_rect_patch_corners_exact() {
    let patch = flat_rect_patch();
    let p00 = patch.point_at(0.0, 0.0);
    assert!(p00.x.abs() < 1e-12 && p00.y.abs() < 1e-12);
    let p10 = patch.point_at(1.0, 0.0);
    assert!((p10.x - 1.0).abs() < 1e-12 && p10.y.abs() < 1e-12);
    let p11 = patch.point_at(1.0, 1.0);
    assert!((p11.x - 1.0).abs() < 1e-12 && (p11.y - 1.0).abs() < 1e-12);
    let p01 = patch.point_at(0.0, 1.0);
    assert!(p01.x.abs() < 1e-12 && (p01.y - 1.0).abs() < 1e-12);
}

/// Build a patch where the inner CP pair at corner (0,0) is
/// deliberately asymmetric: `P_u` and `P_v` at different z
/// heights. On the two boundaries through the corner, only one
/// candidate should contribute (the other blended to zero) so
/// the patch must hit exactly the expected boundary heights.
fn asymmetric_rect_patch() -> RectangularGregory {
    let mut boundary = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in boundary.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = Point3::new(i as f64 / 3.0, j as f64 / 3.0, 0.0);
        }
    }
    // Keep corners (0,1), (1,0), (2,1), (1,2), (2,2), (1,1)
    // other than the one at (0,0) symmetric. Asymmetric pair at
    // corner k=0 only: P_u.z = 0.1, P_v.z = 0.3.
    let asym = (
        Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.1),
        Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.3),
    );
    let sym = (
        Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
        Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
    );
    let sym22 = (
        Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.0),
        Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.0),
    );
    let sym21 = (
        Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
        Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
    );
    RectangularGregory {
        boundary_cps: boundary,
        inner_cps: [asym, sym, sym22, sym21],
    }
}

/// Analytic partial derivatives match finite differences on an
/// asymmetric Gregory patch to within 1e-5 — this is the core
/// correctness check for the quotient-rule `rect_inner_cp_partial_*`
/// implementations.
#[test]
fn rect_partials_agree_with_finite_differences() {
    let patch = asymmetric_rect_patch();
    let h = 1e-5;
    for u in [0.2, 0.4, 0.5, 0.7] {
        for v in [0.2, 0.4, 0.5, 0.7] {
            let du_fd = {
                let pp = patch.point_at(u + h, v);
                let pm = patch.point_at(u - h, v);
                Vec3::new(
                    (pp.x - pm.x) / (2.0 * h),
                    (pp.y - pm.y) / (2.0 * h),
                    (pp.z - pm.z) / (2.0 * h),
                )
            };
            let du_an = patch.partial_u(u, v);
            let dx = (du_an.x - du_fd.x).abs();
            let dy = (du_an.y - du_fd.y).abs();
            let dz = (du_an.z - du_fd.z).abs();
            assert!(
                dx < 1e-5 && dy < 1e-5 && dz < 1e-5,
                "∂/∂u mismatch at ({u}, {v}): analytic ({:.6}, {:.6}, {:.6}) vs FD ({:.6}, {:.6}, {:.6})",
                du_an.x,
                du_an.y,
                du_an.z,
                du_fd.x,
                du_fd.y,
                du_fd.z,
            );

            let dv_fd = {
                let pp = patch.point_at(u, v + h);
                let pm = patch.point_at(u, v - h);
                Vec3::new(
                    (pp.x - pm.x) / (2.0 * h),
                    (pp.y - pm.y) / (2.0 * h),
                    (pp.z - pm.z) / (2.0 * h),
                )
            };
            let dv_an = patch.partial_v(u, v);
            let dx = (dv_an.x - dv_fd.x).abs();
            let dy = (dv_an.y - dv_fd.y).abs();
            let dz = (dv_an.z - dv_fd.z).abs();
            assert!(
                dx < 1e-5 && dy < 1e-5 && dz < 1e-5,
                "∂/∂v mismatch at ({u}, {v}): analytic ({:.6}, {:.6}, {:.6}) vs FD ({:.6}, {:.6}, {:.6})",
                dv_an.x,
                dv_an.y,
                dv_an.z,
                dv_fd.x,
                dv_fd.y,
                dv_fd.z,
            );
        }
    }
}

/// Along the u=0 boundary (v varies), the patch value must be a
/// pure cubic Bezier in `v` through `boundary_cps[0][0..4]` —
/// independent of the inner-CP asymmetry. This is the defining
/// property that makes Gregory stitch G1 across boundaries.
#[test]
fn u_zero_boundary_is_cubic_bezier_in_v() {
    let patch = asymmetric_rect_patch();
    for v in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let p = patch.point_at(0.0, v);
        // Expected: cubic Bezier through (0,0,0) (0,1/3,0) (0,2/3,0) (0,1,0).
        // On the flat boundary, z = 0; y = v; x = 0.
        assert!(p.x.abs() < 1e-12, "x should be 0 at u=0, got {}", p.x);
        assert!(p.z.abs() < 1e-12, "z should be 0 at u=0, got {}", p.z);
        assert!((p.y - v).abs() < 1e-10, "y mismatch at v={v}: {}", p.y);
    }
}

/// Along the v=0 boundary (u varies), same invariant in the u
/// direction.
#[test]
fn v_zero_boundary_is_cubic_bezier_in_u() {
    let patch = asymmetric_rect_patch();
    for u in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        let p = patch.point_at(u, 0.0);
        assert!(p.y.abs() < 1e-12, "y should be 0 at v=0, got {}", p.y);
        assert!(p.z.abs() < 1e-12, "z should be 0 at v=0, got {}", p.z);
        assert!((p.x - u).abs() < 1e-10, "x mismatch at u={u}: {}", p.x);
    }
}
