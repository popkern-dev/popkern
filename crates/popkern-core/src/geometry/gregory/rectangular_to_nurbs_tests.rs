//! Tests for `RectangularGregory::to_nurbs_exact`.
//!
//! The conversion is polynomial-identity exact in the interior and
//! ε-regularized at the four corners. Tests verify: (1) corner points
//! match analytically, (2) boundary curves match to f64 precision,
//! (3) interior evaluation matches `RectangularGregory::point_at` to
//! 1e-10 at a grid of test points, (4) a flat (planar) Gregory
//! produces an exactly-planar NURBS.

use super::rectangular::RectangularGregory;
use super::*;
use crate::geometry::surface::SurfaceEval;

/// A simple flat Gregory patch: unit square on z=0 plane, all inner CPs
/// set to the bilinear-interior values so there's no rational twist.
fn flat_unit_square_gregory() -> RectangularGregory {
    let lin = |u: f64, v: f64| Point3::new(u, v, 0.0);
    let mut boundary_cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in boundary_cps.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = lin(i as f64 / 3.0, j as f64 / 3.0);
        }
    }
    // Inner CPs: bilinear interior values at the 4 inner Bezier indices.
    // Both P_u and P_v at each corner equal the bilinear point so the
    // rational blend collapses to that value.
    let inner_cps = [
        (lin(1.0 / 3.0, 1.0 / 3.0), lin(1.0 / 3.0, 1.0 / 3.0)),
        (lin(1.0 / 3.0, 2.0 / 3.0), lin(1.0 / 3.0, 2.0 / 3.0)),
        (lin(2.0 / 3.0, 2.0 / 3.0), lin(2.0 / 3.0, 2.0 / 3.0)),
        (lin(2.0 / 3.0, 1.0 / 3.0), lin(2.0 / 3.0, 1.0 / 3.0)),
    ];
    RectangularGregory {
        boundary_cps,
        inner_cps,
    }
}

/// A mildly-twisted Gregory patch: flat boundaries but non-bilinear
/// inner CP pairs so the rational blend is non-trivial.
fn twisted_gregory() -> RectangularGregory {
    let mut g = flat_unit_square_gregory();
    // Offset the inner_cps so P_u and P_v differ — this exercises
    // the rational blend.
    g.inner_cps[0].0 = Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.15);
    g.inner_cps[0].1 = Point3::new(1.0 / 3.0, 1.0 / 3.0, -0.05);
    g.inner_cps[1].0 = Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.05);
    g.inner_cps[1].1 = Point3::new(1.0 / 3.0, 2.0 / 3.0, 0.10);
    g.inner_cps[2].0 = Point3::new(2.0 / 3.0, 2.0 / 3.0, -0.08);
    g.inner_cps[2].1 = Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.12);
    g.inner_cps[3].0 = Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.03);
    g.inner_cps[3].1 = Point3::new(2.0 / 3.0, 1.0 / 3.0, -0.11);
    g
}

fn max_abs(a: Point3, b: Point3) -> f64 {
    (a.x - b.x)
        .abs()
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

#[test]
fn to_nurbs_exact_corners_flat_match() {
    let g = flat_unit_square_gregory();
    let nurbs = g.to_nurbs_exact();
    for &(u, v, expect_idx_i, expect_idx_j) in &[
        (0.0_f64, 0.0_f64, 0usize, 0usize),
        (0.0, 1.0, 0, 3),
        (1.0, 0.0, 3, 0),
        (1.0, 1.0, 3, 3),
    ] {
        let expected = g.boundary_cps[expect_idx_i][expect_idx_j];
        let got = nurbs.point_at(u, v);
        assert!(
            max_abs(got, expected) < 1e-10,
            "corner (u,v)=({u}, {v}): expected {expected:?}, got {got:?}"
        );
    }
}

#[test]
fn to_nurbs_exact_flat_stays_planar() {
    let g = flat_unit_square_gregory();
    let nurbs = g.to_nurbs_exact();
    for &u in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
        for &v in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let p = nurbs.point_at(u, v);
            assert!(
                p.z.abs() < 1e-10,
                "flat Gregory's NURBS image should be planar; got z = {} at ({u}, {v})",
                p.z
            );
        }
    }
}

#[test]
fn to_nurbs_exact_flat_matches_gregory_eval() {
    let g = flat_unit_square_gregory();
    let nurbs = g.to_nurbs_exact();
    for &u in &[0.0, 0.1, 0.37, 0.5, 0.63, 0.9, 1.0] {
        for &v in &[0.0, 0.1, 0.37, 0.5, 0.63, 0.9, 1.0] {
            let p_g = g.point_at(u, v);
            let p_n = nurbs.point_at(u, v);
            assert!(
                max_abs(p_g, p_n) < 1e-10,
                "flat Gregory vs NURBS disagree at ({u}, {v}): {p_g:?} vs {p_n:?}"
            );
        }
    }
}

#[test]
fn to_nurbs_exact_twisted_matches_gregory_eval() {
    let g = twisted_gregory();
    let nurbs = g.to_nurbs_exact();
    // Sample interior points; corners are handled by the ε-regularized path
    // which is covered separately in `corners_flat_match`.
    let params = [0.0_f64, 0.1, 0.25, 0.37, 0.5, 0.63, 0.75, 0.9, 1.0];
    let mut max_err = 0.0_f64;
    for &u in &params {
        for &v in &params {
            let p_g = g.point_at(u, v);
            let p_n = nurbs.point_at(u, v);
            let err = max_abs(p_g, p_n);
            if err > max_err {
                max_err = err;
            }
        }
    }
    assert!(
        max_err < 1e-10,
        "twisted Gregory vs NURBS max error = {max_err:e}, expected < 1e-10"
    );
}

#[test]
fn to_nurbs_exact_shape_metadata() {
    let g = twisted_gregory();
    let nurbs = g.to_nurbs_exact();
    assert_eq!(nurbs.degree_u, 7);
    assert_eq!(nurbs.degree_v, 7);
    assert_eq!(nurbs.num_u, 8);
    assert_eq!(nurbs.num_v, 8);
    assert_eq!(nurbs.control_points.len(), 64);
    assert_eq!(nurbs.weights.len(), 64);
    assert_eq!(nurbs.knots_u.len(), 16);
    assert_eq!(nurbs.knots_v.len(), 16);
    // All weights must be positive (ε-regularization guarantee).
    for &w in &nurbs.weights {
        assert!(w > 0.0, "all weights must be positive post-regularization");
    }
}

// ---- RectangularGregoryDeg5::to_nurbs_exact ----

use super::build::GregoryBuildError;
use super::rectangular_deg5::RectangularGregoryDeg5;

fn flat_deg5_boundaries_unit_square() -> [[Point3; 6]; 4] {
    let c00 = Point3::new(0.0, 0.0, 0.0);
    let c05 = Point3::new(0.0, 1.0, 0.0);
    let c55 = Point3::new(1.0, 1.0, 0.0);
    let c50 = Point3::new(1.0, 0.0, 0.0);
    let lin = |a: Point3, b: Point3, t: f64| -> Point3 {
        Point3::new(
            (1.0 - t) * a.x + t * b.x,
            (1.0 - t) * a.y + t * b.y,
            (1.0 - t) * a.z + t * b.z,
        )
    };
    let b0: [Point3; 6] = core::array::from_fn(|i| lin(c00, c05, i as f64 / 5.0));
    let b1: [Point3; 6] = core::array::from_fn(|i| lin(c05, c55, i as f64 / 5.0));
    let b2: [Point3; 6] = core::array::from_fn(|i| lin(c55, c50, i as f64 / 5.0));
    let b3: [Point3; 6] = core::array::from_fn(|i| lin(c50, c00, i as f64 / 5.0));
    [b0, b1, b2, b3]
}

#[test]
fn deg5_to_nurbs_exact_corners_match() {
    let boundaries = flat_deg5_boundaries_unit_square();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    let nurbs = patch.to_nurbs_exact();
    for &(u, v, exp) in &[
        (0.0, 0.0, Point3::new(0.0, 0.0, 0.0)),
        (0.0, 1.0, Point3::new(0.0, 1.0, 0.0)),
        (1.0, 0.0, Point3::new(1.0, 0.0, 0.0)),
        (1.0, 1.0, Point3::new(1.0, 1.0, 0.0)),
    ] {
        let got = nurbs.point_at(u, v);
        assert!(
            max_abs(got, exp) < 1e-10,
            "deg5 corner ({u}, {v}): expected {exp:?}, got {got:?}"
        );
    }
}

#[test]
fn deg5_to_nurbs_exact_flat_stays_planar() {
    let boundaries = flat_deg5_boundaries_unit_square();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    let nurbs = patch.to_nurbs_exact();
    for &u in &[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
        for &v in &[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            let p = nurbs.point_at(u, v);
            assert!(
                p.z.abs() < 1e-10,
                "flat deg5 Gregory's NURBS image should be planar; got z = {} at ({u}, {v})",
                p.z
            );
        }
    }
}

#[test]
fn deg5_to_nurbs_exact_flat_matches_gregory() {
    let boundaries = flat_deg5_boundaries_unit_square();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    let nurbs = patch.to_nurbs_exact();
    for &u in &[0.0, 0.17, 0.4, 0.5, 0.63, 0.9, 1.0] {
        for &v in &[0.0, 0.17, 0.4, 0.5, 0.63, 0.9, 1.0] {
            let p_g = patch.point_at(u, v);
            let p_n = nurbs.point_at(u, v);
            assert!(
                max_abs(p_g, p_n) < 1e-10,
                "deg5 flat Gregory vs NURBS disagree at ({u}, {v}): {p_g:?} vs {p_n:?}"
            );
        }
    }
}

#[test]
fn deg5_to_nurbs_exact_twisted_matches_gregory() {
    // Twisted patch: non-zero inner CPs per corner with both P_u and
    // P_v differing so the rational blend is non-trivial.
    let boundaries = flat_deg5_boundaries_unit_square();
    let mut patch =
        RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    // Perturb inner CPs with different (P_u, P_v) offsets in z.
    patch.inner_cps[0] = (Point3::new(0.2, 0.2, 0.15), Point3::new(0.2, 0.2, -0.05));
    patch.inner_cps[1] = (Point3::new(0.2, 0.8, 0.07), Point3::new(0.2, 0.8, 0.12));
    patch.inner_cps[2] = (Point3::new(0.8, 0.8, -0.09), Point3::new(0.8, 0.8, 0.14));
    patch.inner_cps[3] = (Point3::new(0.8, 0.2, 0.04), Point3::new(0.8, 0.2, -0.13));

    let nurbs = patch.to_nurbs_exact();
    let params = [0.0_f64, 0.1, 0.25, 0.4, 0.5, 0.63, 0.75, 0.9, 1.0];
    let mut max_err = 0.0_f64;
    for &u in &params {
        for &v in &params {
            let p_g = patch.point_at(u, v);
            let p_n = nurbs.point_at(u, v);
            let err = max_abs(p_g, p_n);
            if err > max_err {
                max_err = err;
            }
        }
    }
    assert!(
        max_err < 1e-10,
        "deg5 twisted Gregory vs NURBS max error = {max_err:e}, expected < 1e-10"
    );
}

#[test]
fn deg5_to_nurbs_exact_shape_metadata() {
    let boundaries = flat_deg5_boundaries_unit_square();
    let patch = RectangularGregoryDeg5::from_boundary_beziers_no_twist(boundaries, 1e-9).unwrap();
    let nurbs = patch.to_nurbs_exact();
    assert_eq!(nurbs.degree_u, 9);
    assert_eq!(nurbs.degree_v, 9);
    assert_eq!(nurbs.num_u, 10);
    assert_eq!(nurbs.num_v, 10);
    assert_eq!(nurbs.control_points.len(), 100);
    assert_eq!(nurbs.weights.len(), 100);
    for &w in &nurbs.weights {
        assert!(
            w > 0.0,
            "all deg5 weights must be positive post-regularization"
        );
    }
}

// Silence unused-import warning for GregoryBuildError, pulled in for
// documentation parity with the cubic tests.
#[allow(dead_code)]
fn _touch_build_error(_: GregoryBuildError) {}
