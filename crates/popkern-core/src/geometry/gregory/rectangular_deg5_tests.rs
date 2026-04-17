//! Tests for the degree-5 rectangular Gregory patch evaluator.
//!
//! Focus: flat-plane identity, boundary interpolation, derivative
//! agreement with finite differences, and corner-exact position.

use super::rectangular_deg5::RectangularGregoryDeg5;
use super::*;

fn flat_deg5_boundaries() -> [[Point3; 6]; 4] {
    // A unit square in the z=0 plane with 5-segment spacing.
    let corner_00 = Point3::new(0.0, 0.0, 0.0);
    let corner_05 = Point3::new(0.0, 1.0, 0.0);
    let corner_55 = Point3::new(1.0, 1.0, 0.0);
    let corner_50 = Point3::new(1.0, 0.0, 0.0);

    let lin = |a: Point3, b: Point3, t: f64| -> Point3 {
        Point3::new(
            (1.0 - t) * a.x + t * b.x,
            (1.0 - t) * a.y + t * b.y,
            (1.0 - t) * a.z + t * b.z,
        )
    };

    // boundary 0: corner 0 → corner 1, 6 CPs at t = 0, 1/5, ..., 1.
    let b0: [Point3; 6] = core::array::from_fn(|i| lin(corner_00, corner_05, i as f64 / 5.0));
    // boundary 1: corner 1 → corner 2
    let b1: [Point3; 6] = core::array::from_fn(|i| lin(corner_05, corner_55, i as f64 / 5.0));
    // boundary 2: corner 2 → corner 3
    let b2: [Point3; 6] = core::array::from_fn(|i| lin(corner_55, corner_50, i as f64 / 5.0));
    // boundary 3: corner 3 → corner 0
    let b3: [Point3; 6] = core::array::from_fn(|i| lin(corner_50, corner_00, i as f64 / 5.0));

    [b0, b1, b2, b3]
}

#[test]
fn deg5_no_twist_builds_on_unit_square() {
    let boundaries = flat_deg5_boundaries();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9)
        .expect("no-twist builder must succeed on closed CCW square");
    // Corner positions must match.
    let p00 = patch.point_at(0.0, 0.0);
    assert!((p00.x - 0.0).abs() < 1e-12);
    assert!((p00.y - 0.0).abs() < 1e-12);
    let p11 = patch.point_at(1.0, 1.0);
    assert!((p11.x - 1.0).abs() < 1e-12);
    assert!((p11.y - 1.0).abs() < 1e-12);
}

#[test]
fn deg5_no_twist_stays_planar_on_flat_boundaries() {
    // With all 4 boundaries in z=0 and no-twist inner CPs
    // (bilinear interp of corners — also in z=0), every
    // evaluation must have z = 0 to within f64 precision.
    let boundaries = flat_deg5_boundaries();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    for &u in &[0.1, 0.3, 0.5, 0.7, 0.9] {
        for &v in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let p = patch.point_at(u, v);
            assert!(
                p.z.abs() < 1e-12,
                "deg5 no-twist on flat boundaries must stay in z=0 plane; \
                 got z = {} at (u, v) = ({u}, {v})",
                p.z,
            );
        }
    }
}

#[test]
fn deg5_corner_mismatch_rejected() {
    let mut boundaries = flat_deg5_boundaries();
    // Perturb the shared corner between boundary 0 and boundary 1
    // by more than tol_linear.
    boundaries[1][0].x += 0.1;
    let res = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-6);
    assert!(
        matches!(res, Err(GregoryBuildError::CornerMismatch { .. })),
        "corner mismatch must surface as CornerMismatch error"
    );
}

#[test]
fn deg5_partial_u_matches_finite_difference() {
    // Mildly-warped patch: shift the central row of boundary CPs in z
    // so the patch has non-zero curvature, then verify `partial_u`
    // matches a central-difference `point_at` approximation.
    let mut boundaries = flat_deg5_boundaries();
    // Bump the middle CPs of boundaries 0 and 2 up in z — creates a
    // ridge along v = 0.5 in the patch.
    boundaries[0][3].z = 0.2;
    boundaries[2][2].z = 0.2; // boundary 2 is reversed; index 2 in storage = v close to middle
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    let h = 1e-5;
    let test_points = [(0.3, 0.4), (0.5, 0.5), (0.7, 0.2), (0.2, 0.8)];
    for &(u, v) in &test_points {
        let analytic = patch.partial_u(u, v);
        let fwd = patch.point_at(u + h, v);
        let bwd = patch.point_at(u - h, v);
        let fd_x = (fwd.x - bwd.x) / (2.0 * h);
        let fd_y = (fwd.y - bwd.y) / (2.0 * h);
        let fd_z = (fwd.z - bwd.z) / (2.0 * h);
        assert!(
            (analytic.x - fd_x).abs() < 1e-4,
            "partial_u.x at ({u},{v}): analytic {}, FD {}",
            analytic.x,
            fd_x
        );
        assert!(
            (analytic.y - fd_y).abs() < 1e-4,
            "partial_u.y at ({u},{v}): analytic {}, FD {}",
            analytic.y,
            fd_y
        );
        assert!(
            (analytic.z - fd_z).abs() < 1e-4,
            "partial_u.z at ({u},{v}): analytic {}, FD {}",
            analytic.z,
            fd_z
        );
    }
}

#[test]
fn deg5_partial_v_matches_finite_difference() {
    let mut boundaries = flat_deg5_boundaries();
    boundaries[1][3].z = 0.15;
    boundaries[3][2].z = 0.15;
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    let h = 1e-5;
    for &(u, v) in &[(0.3, 0.4), (0.5, 0.5), (0.7, 0.2), (0.2, 0.8)] {
        let analytic = patch.partial_v(u, v);
        let fwd = patch.point_at(u, v + h);
        let bwd = patch.point_at(u, v - h);
        let fd_x = (fwd.x - bwd.x) / (2.0 * h);
        let fd_y = (fwd.y - bwd.y) / (2.0 * h);
        let fd_z = (fwd.z - bwd.z) / (2.0 * h);
        assert!((analytic.x - fd_x).abs() < 1e-4);
        assert!((analytic.y - fd_y).abs() < 1e-4);
        assert!((analytic.z - fd_z).abs() < 1e-4);
    }
}

/// Boundary interpolation: evaluating on `u = 0` must trace
/// the boundary 0 quintic-Bezier exactly (no contribution from
/// interior CPs on the boundary line).
#[test]
fn deg5_boundary_interpolation() {
    let boundaries = flat_deg5_boundaries();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    // Boundary 0 is corner 0 → corner 1, which in this flat
    // fixture is the line x=0, 0 ≤ y ≤ 1. Sample and verify.
    for &v in &[0.0, 0.1, 0.3, 0.5, 0.7, 1.0] {
        let p = patch.point_at(0.0, v);
        assert!(p.x.abs() < 1e-12, "u=0 boundary: expected x=0, got {}", p.x);
        assert!(
            (p.y - v).abs() < 1e-12,
            "u=0 boundary: expected y={v}, got {}",
            p.y
        );
        assert!(p.z.abs() < 1e-12);
    }
}
