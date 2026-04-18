//! Tests for `TriangularGregory::to_nurbs_exact`.
//!
//! Validates the bidegree-(5, 5) rational NURBS conversion matches
//! the analytic triangular evaluator on corners, boundaries, and a
//! grid of interior points. The `s = 1` edge of the rectangular
//! domain collapses onto triangular corner 0 (degenerate edge) —
//! tested separately.

use super::triangular::TriangularGregory;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;

/// Flat triangular Gregory on the unit barycentric simplex (corners
/// at standard R² positions). Inner CPs pinned to the centroid so
/// the rational blend collapses to the Bezier mean.
fn flat_triangular() -> TriangularGregory {
    let c0 = Point3::new(1.0, 0.0, 0.0);
    let c1 = Point3::new(0.0, 1.0, 0.0);
    let c2 = Point3::new(0.0, 0.0, 1.0);
    // Edge CPs at one-third positions on each edge (boundary order:
    // [300, 030, 003, 210, 120, 021, 012, 102, 201]).
    let lerp = |a: Point3, b: Point3, t: f64| -> Point3 {
        Point3::new(
            (1.0 - t) * a.x + t * b.x,
            (1.0 - t) * a.y + t * b.y,
            (1.0 - t) * a.z + t * b.z,
        )
    };
    let centroid = Point3::new(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    TriangularGregory {
        boundary_cps: [
            c0,                      // 300
            c1,                      // 030
            c2,                      // 003
            lerp(c0, c1, 1.0 / 3.0), // 210
            lerp(c0, c1, 2.0 / 3.0), // 120
            lerp(c1, c2, 1.0 / 3.0), // 021
            lerp(c1, c2, 2.0 / 3.0), // 012
            lerp(c2, c0, 1.0 / 3.0), // 102
            lerp(c2, c0, 2.0 / 3.0), // 201
        ],
        inner_cps: [centroid, centroid, centroid],
    }
}

fn max_abs(a: Point3, b: Point3) -> f64 {
    (a.x - b.x)
        .abs()
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

#[test]
fn triangular_to_nurbs_exact_corners_match() {
    let g = flat_triangular();
    let nurbs = g.to_nurbs_exact();
    // (s=1, t=*): degenerate edge → triangular corner 0.
    // (s=0, t=0) → corner 1; (s=0, t=1) → corner 2.
    let cases = [
        (1.0_f64, 0.0_f64, g.boundary_cps[0]), // corner 0
        (1.0, 0.5, g.boundary_cps[0]),         // middle of degenerate edge
        (1.0, 1.0, g.boundary_cps[0]),         // other end of degenerate edge
        (0.0, 0.0, g.boundary_cps[1]),         // corner 1
        (0.0, 1.0, g.boundary_cps[2]),         // corner 2
    ];
    for (s, t, exp) in cases {
        let got = nurbs.point_at(s, t);
        assert!(
            max_abs(got, exp) < 1e-10,
            "corner (s, t) = ({s}, {t}): expected {exp:?}, got {got:?}"
        );
    }
}

#[test]
fn triangular_to_nurbs_exact_matches_analytic_interior() {
    let g = flat_triangular();
    let nurbs = g.to_nurbs_exact();
    // Interior (s, t) grid. Map each to (u, v, w) = (s, (1-s)(1-t), (1-s)t).
    let samples = [0.1_f64, 0.25, 0.4, 0.6, 0.75, 0.9];
    let mut max_err = 0.0_f64;
    for &s in &samples {
        for &t in &samples {
            let u = s;
            let v = (1.0 - s) * (1.0 - t);
            let _w = (1.0 - s) * t;
            // Triangular Gregory is parameterised by (u, v) with
            // w = 1 - u - v implicit.
            let p_tri = g.point_at_uv(u, v);
            let p_nurbs = nurbs.point_at(s, t);
            let err = max_abs(p_tri, p_nurbs);
            if err > max_err {
                max_err = err;
            }
        }
    }
    assert!(
        max_err < 1e-10,
        "triangular Gregory vs NURBS max error = {max_err:e}, expected < 1e-10"
    );
}

#[test]
fn triangular_to_nurbs_exact_stays_planar_on_planar_patch() {
    let g = flat_triangular();
    let nurbs = g.to_nurbs_exact();
    // Every boundary + inner CP lies on the plane x + y + z = 1.
    // Both the triangular eval and the NURBS image should live on
    // that plane.
    let samples = [0.0_f64, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0];
    for &s in &samples {
        for &t in &samples {
            let p = nurbs.point_at(s, t);
            let plane_dev = (p.x + p.y + p.z - 1.0).abs();
            assert!(
                plane_dev < 1e-10,
                "flat triangular NURBS must lie on x+y+z=1; got deviation {plane_dev:e} at (s, t) = ({s}, {t})"
            );
        }
    }
}

#[test]
fn triangular_to_nurbs_exact_twisted_matches_analytic() {
    // Non-trivial inner-CP offsets so the rational blend is active.
    let mut g = flat_triangular();
    g.inner_cps[0] = Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.4);
    g.inner_cps[1] = Point3::new(0.35, 1.0 / 3.0 + 0.1, 1.0 / 3.0);
    g.inner_cps[2] = Point3::new(1.0 / 3.0 + 0.05, 0.3, 0.4);
    let nurbs = g.to_nurbs_exact();
    let samples = [0.1_f64, 0.25, 0.4, 0.55, 0.7, 0.85];
    let mut max_err = 0.0_f64;
    for &s in &samples {
        for &t in &samples {
            let u = s;
            let v = (1.0 - s) * (1.0 - t);
            let p_tri = g.point_at_uv(u, v);
            let p_nurbs = nurbs.point_at(s, t);
            let err = max_abs(p_tri, p_nurbs);
            if err > max_err {
                max_err = err;
            }
        }
    }
    assert!(
        max_err < 1e-10,
        "twisted triangular Gregory vs NURBS max error = {max_err:e}"
    );
}

#[test]
fn triangular_to_nurbs_exact_shape_metadata() {
    let g = flat_triangular();
    let nurbs = g.to_nurbs_exact();
    assert_eq!(nurbs.degree_u, 5);
    assert_eq!(nurbs.degree_v, 5);
    assert_eq!(nurbs.num_u, 6);
    assert_eq!(nurbs.num_v, 6);
    assert_eq!(nurbs.control_points.len(), 36);
    assert_eq!(nurbs.weights.len(), 36);
    assert_eq!(nurbs.knots_u.len(), 12);
    assert_eq!(nurbs.knots_v.len(), 12);
    for &w in &nurbs.weights {
        assert!(w > 0.0, "all weights must be positive post-regularization");
    }
}
