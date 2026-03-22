//! Tests for the Bezier clipping intersection algorithm.

#[cfg(test)]
use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;

fn p(x: f64, y: f64) -> Point3 {
    Point3::new(x, y, 0.0)
}

#[test]
fn clip_edge_cases() {
    let h = vec![p(0.0, 0.0), p(1.0, 0.0)];
    // Non-overlapping
    assert!(clip_against_fat_line(&h, &[p(0.0, 2.0), p(1.0, 2.0)], (0.0, 1.0)).is_none());
    // Empty/short inputs
    assert!(clip_against_fat_line(&[p(0.0, 0.0)], &h, (0.0, 1.0)).is_none());
    assert!(clip_against_fat_line(&h, &vec![], (0.0, 1.0)).is_none());
    // Degenerate A (coincident)
    assert!(clip_against_fat_line(&[p(0.5, 0.5), p(0.5, 0.5)], &h, (0.0, 1.0)).is_none());
}

#[test]
fn clip_crossing_lines() {
    let a = vec![p(0.0, 0.0), p(1.0, 0.0)];
    let b = vec![p(0.0, -1.0), p(1.0, 1.0)];
    let (t_min, t_max) = clip_against_fat_line(&a, &b, (0.0, 1.0)).unwrap();
    assert!(t_min <= 0.5 && t_max >= 0.5);
}

#[test]
fn subdivide_cases() {
    // Single point
    assert_eq!(subdivide_bezier(&[p(1.0, 2.0)], 0.0, 1.0).len(), 1);
    // Full range returns original
    let pts = vec![p(0.0, 0.0), p(0.5, 1.0), p(1.0, 0.0)];
    assert_eq!(subdivide_bezier(&pts, 0.0, 1.0).len(), 3);
    // t_hi = 0 degenerate
    assert_eq!(
        subdivide_bezier(&[p(0.0, 0.0), p(1.0, 0.0)], 0.0, 0.0).len(),
        2
    );
    // Quadratic left half
    let left = subdivide_bezier(&pts, 0.0, 0.5);
    assert_eq!(left.len(), 3);
    assert!((left[1].x - 0.25).abs() < 1e-10 && (left[1].y - 0.5).abs() < 1e-10);
}

fn run_intersect(
    a: &[Point3],
    aw: &[f64],
    b: &[Point3],
    bw: &[f64],
    depth: usize,
) -> Vec<(f64, f64)> {
    let mut r = Vec::new();
    bezier_intersect(
        a,
        aw,
        (0.0, 1.0),
        b,
        bw,
        (0.0, 1.0),
        &Tolerance::default(),
        depth,
        &mut r,
    );
    r
}

#[test]
fn bezier_intersect_cases() {
    let arc = vec![p(0.0, 0.0), p(0.5, 1.0), p(1.0, 0.0)];
    let line = vec![p(0.0, 0.5), p(1.0, 0.5)];
    assert!(!run_intersect(&arc, &[1.0; 3], &line, &[1.0; 2], 50).is_empty());
    // Non-intersecting
    assert!(
        run_intersect(
            &[p(0.0, 0.0), p(1.0, 0.0)],
            &[1.0; 2],
            &[p(0.0, 5.0), p(1.0, 5.0)],
            &[1.0; 2],
            50
        )
        .is_empty()
    );
    // Rational weights
    assert!(
        !run_intersect(
            &arc,
            &[1.0, 2.0, 1.0],
            &[p(0.0, 0.25), p(1.0, 0.25)],
            &[1.0; 2],
            50
        )
        .is_empty()
    );
    // max_depth = 0 fallback
    assert!(!run_intersect(&arc, &[1.0; 3], &line, &[1.0; 2], 0).is_empty());
}

#[test]
fn compute_fat_line_normal_collinear_fallback() {
    let dir = crate::geometry::Vec3::new(1.0, 0.0, 0.0);
    assert!(
        compute_fat_line_normal(
            &dir,
            &[p(0.0, 0.0), p(1.0, 0.0)],
            &[p(0.0, 0.0), p(2.0, 0.0)]
        )
        .is_some()
    );
}

#[test]
fn evaluate_bezier_midpoint_linear() {
    // Midpoint of a linear segment from (0,0) to (1,0) should be (0.5, 0)
    let pts = vec![p(0.0, 0.0), p(1.0, 0.0)];
    let mid = evaluate_bezier_midpoint(&pts);
    assert!((mid.x - 0.5).abs() < 1e-12);
    assert!(mid.y.abs() < 1e-12);
}

#[test]
fn evaluate_bezier_midpoint_quadratic() {
    // Midpoint of parabola arc: control points (0,0), (0.5,1), (1,0).
    // De Casteljau at t=0.5: result = (0.5, 0.5)
    let pts = vec![p(0.0, 0.0), p(0.5, 1.0), p(1.0, 0.0)];
    let mid = evaluate_bezier_midpoint(&pts);
    assert!((mid.x - 0.5).abs() < 1e-12);
    assert!((mid.y - 0.5).abs() < 1e-12);
}

#[test]
fn evaluate_bezier_midpoint_empty() {
    // Empty points should not panic, returns origin
    let mid = evaluate_bezier_midpoint(&[]);
    assert_eq!(mid.x, 0.0);
    assert_eq!(mid.y, 0.0);
    assert_eq!(mid.z, 0.0);
}

#[test]
fn near_tangent_curves_both_found() {
    // A quadratic Bezier that starts above y=0, dips just barely below y=0 near
    // its center, then rises back above y=0.  This creates two near-tangent
    // crossings of the y=0 line at interior parameter values (t≈0.42 and t≈0.578).
    //
    // Control points: (0, 0.1), (0.5, -0.105), (1, 0.1)
    // The minimum y is about -0.0025 at t=0.5, giving two close crossings.
    let curve_pts = vec![p(0.0, 0.1), p(0.5, -0.105), p(1.0, 0.1)];
    let line_pts = vec![p(0.0, 0.0), p(1.0, 0.0)];
    let weights_3 = vec![1.0, 1.0, 1.0];
    let weights_2 = vec![1.0, 1.0];

    let mut results = Vec::new();
    bezier_intersect(
        &curve_pts,
        &weights_3,
        (0.0, 1.0),
        &line_pts,
        &weights_2,
        (0.0, 1.0),
        &Tolerance::default(),
        60,
        &mut results,
    );

    assert!(
        results.len() >= 2,
        "Should find both near-tangent crossings, got {}",
        results.len()
    );

    // Verify two distinct crossings — one near t≈0.42 (left) and one near t≈0.578 (right).
    let has_left = results.iter().any(|(t1, _)| *t1 < 0.5);
    let has_right = results.iter().any(|(t1, _)| *t1 > 0.5);
    assert!(
        has_left,
        "Expected a crossing in left half (t<0.5), results: {:?}",
        results
    );
    assert!(
        has_right,
        "Expected a crossing in right half (t>0.5), results: {:?}",
        results
    );
}

#[test]
fn stall_bisection_triggered_on_symmetric_curves() {
    // Two identical quadratic arcs. The fat-line clip gives poor reduction
    // (since they overlap everywhere). The stall counter should trigger
    // double-bisection and still terminate without panic or runaway.
    let arc = vec![p(0.0, 0.0), p(0.5, 1.0), p(1.0, 0.0)];
    let weights = vec![1.0, 1.0, 1.0];
    let mut results = Vec::new();
    // This should not panic or infinite-loop; any outcome is acceptable.
    bezier_intersect(
        &arc,
        &weights,
        (0.0, 1.0),
        &arc,
        &weights,
        (0.0, 1.0),
        &Tolerance::default(),
        30,
        &mut results,
    );
    // At minimum, the algorithm must terminate. No assertion on count needed.
}
