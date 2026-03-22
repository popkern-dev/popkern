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

#[test]
fn symmetric_convergence_base_case() {
    // Lines 314-318: hit the base case where BOTH a_span and b_span are below
    // parametric tolerance on the SAME recursive call entry.
    //
    // Construction: A = parabola (0,0)→(0.5,0.5)→(1,0) with fat-line
    // half-thickness d_max=0.25.  B = vertical line (0.3,-1)→(0.3,1).
    // parametric = 0.5.
    //
    // On the first call (both spans=1.0 > 0.5) the tangent check fires but
    // finds large midpoint distance (A-mid=(0.5,0.5), B-mid=(0.3,0), dist≈0.54)
    // so it falls through to clipping.
    //
    // First clip narrows B to the slice crossing A's fat-line:
    //   clipped_b_span ≈ 0.25 < parametric=0.5.
    // Reverse clip (A against that tiny B slice's fat-line, normal≈(1,0)) finds
    // A's Bernstein polygon crossing x=0.3 at a single point:
    //   clipped_a_span ≈ 1e-14 < parametric=0.5.
    //
    // The recursive call therefore enters with BOTH a_span < 0.5 AND
    // b_span < 0.5, hitting the symmetric base case at lines 314-318.
    // A peak at y=0.25 gives d_max=0.25 → clipped_b_span=0.25, which is strictly
    // less than parametric=0.5.  (A peak at y=0.5 would give clipped_b_span=0.5
    // which is NOT strictly less than parametric, so line 314 would not fire.)
    let a = vec![p(0.0, 0.0), p(0.5, 0.25), p(1.0, 0.0)]; // parabola, d_max=0.25
    let b = vec![p(0.3, -1.0), p(0.3, 1.0)]; // vertical line, crosses A near x=0.3
    let tol = Tolerance::default().with_parametric(0.5);
    let mut results = Vec::new();
    bezier_intersect(
        &a,
        &[1.0, 1.0, 1.0],
        (0.0, 1.0),
        &b,
        &[1.0, 1.0],
        (0.0, 1.0),
        &tol,
        50,
        &mut results,
    );
    // The parabola A crosses the vertical line B at y=0, x=0.3.
    assert!(
        !results.is_empty(),
        "Should find crossing of parabola with vertical line"
    );
}

#[test]
fn asymmetric_b_convergence_base_case() {
    // Lines 331-335: b_span < parametric but a_span is not.
    //
    // Construction: A = long horizontal line (0,0)→(4,0), range (0,1).
    // B = short vertical line (3,-0.5)→(3,0.5), range (0,1).
    // Their midpoints A(0.5)=(2,0) and B(0.5)=(3,0) are 1.0 apart → distance
    // > linear=1e-6, so the tangent check at line 349 fires but fails.
    //
    // First clip narrows B to the crossing of A's thin fat-line:
    //   clipped_b_span ≈ 1e-14, which is < parametric=0.6.
    //
    // Reverse clip (A against tiny near-degenerate B slice) returns
    // Some(full a_range=(0,1)) because the degenerate B slice's fat line
    // degenerates to the Z-axis normal, giving signed distance ≈ 0 for all
    // XY-plane points — clipped_a_span = 1.0, which is > parametric=0.6.
    //
    // Recursive call therefore has a_span=1.0 (not < 0.6) and
    // b_span≈1e-14 (< 0.6), hitting the asymmetric-B base case at lines 331-335.
    let a = vec![p(0.0, 0.0), p(4.0, 0.0)]; // long horizontal, midpoint at (2,0)
    let b = vec![p(3.0, -0.5), p(3.0, 0.5)]; // vertical at x=3, midpoint at (3,0)
    let tol = Tolerance::default().with_parametric(0.6);
    let mut results = Vec::new();
    bezier_intersect(
        &a,
        &[1.0, 1.0],
        (0.0, 1.0),
        &b,
        &[1.0, 1.0],
        (0.0, 1.0),
        &tol,
        50,
        &mut results,
    );
    assert!(
        !results.is_empty(),
        "Should find intersection of long horizontal and vertical lines"
    );
}

#[test]
fn tangent_termination_base_case() {
    // Lines 349-359: both spans within 10×parametric but not within parametric,
    // and the midpoints of the two curves are within linear tolerance.
    //
    // Strategy: use two nearly-identical parabolas and very small initial
    // parameter ranges so the spans start in (parametric, 10*parametric).
    // The fat line clip for coincident curves gives b_reduction=1.0 (stall),
    // which after 3 consecutive stalls triggers the 4-way bisect. Each
    // bisected sub-call starts with span = a_span/2 still in
    // (parametric, 10*parametric), so the tangent check fires.
    //
    // Use parametric=1e-9 and initial ranges of size 5e-9 so each span
    // satisfies: 5e-9 > 1e-9 (parametric) and 5e-9 < 1e-8 (10*parametric).
    // Use linear=1.0 so any spatial distance < 1.0 counts as "close enough".
    let a = vec![p(0.0, 0.0), p(0.5, 0.01), p(1.0, 0.0)]; // parabola
    let b = vec![p(0.0, 0.0), p(0.5, 0.01), p(1.0, 0.0)]; // identical — midpoints coincide
    let tol = Tolerance::default().with_parametric(1e-9).with_linear(1.0);
    let mut results = Vec::new();
    bezier_intersect(
        &a,
        &[1.0, 1.0, 1.0],
        (0.0, 5e-9), // span = 5e-9, in (1e-9, 1e-8)
        &b,
        &[1.0, 1.0, 1.0],
        (0.0, 5e-9),
        &tol,
        50,
        &mut results,
    );
    // The coincident curves trigger the tangent-termination path at lines 349-359.
    assert!(
        !results.is_empty(),
        "Tangent termination should report a result for coincident curves"
    );
}

#[test]
fn reverse_clip_returns_none() {
    // Line 387: after clipping B against A's fat line succeeds (B passes through
    // A's fat line band), clipping A against the clipped B's fat line returns
    // None because A lies entirely to one side of B's fat line.
    //
    // Construction:
    //   A = parabola from (0,0) to (1,0) with peak at (0.5, 0.5) — width 1,
    //       fat-line half-thickness d_max = 0.25.
    //   B = vertical line from (5,-1) to (5,1) — at x=5, does NOT intersect A.
    //
    // First clip (B against A's fat line, which measures Y-distance):
    //   B crosses A's Y-band [-0.25, +0.25] near y=0, giving clipped_b_range
    //   of roughly (0.375, 0.625) — span = 0.25, well above parametric=1e-8.
    //
    // Reverse clip (A against clipped-B's fat line, which measures X-distance):
    //   b_sub2 ≈ vertical segment at x=5 near y=0; its fat-line direction is
    //   (0,1,0) and normal is (1,0,0). All of A's control points are at x≤1,
    //   giving signed distances ≈ -4 to -5 — entirely outside b_sub2's band
    //   (half-thickness ~ 0) → clip_against_fat_line returns None → line 387.
    let a = vec![p(0.0, 0.0), p(0.5, 0.5), p(1.0, 0.0)]; // parabola at x∈[0,1]
    let b = vec![p(5.0, -1.0), p(5.0, 1.0)]; // vertical line at x=5
    let tol = Tolerance::default();
    let mut results = Vec::new();
    bezier_intersect(
        &a,
        &[1.0, 1.0, 1.0],
        (0.0, 1.0),
        &b,
        &[1.0, 1.0],
        (0.0, 1.0),
        &tol,
        50,
        &mut results,
    );
    // These curves do not intersect — A is at x∈[0,1], B is at x=5.
    // The important thing is that line 387 (None return from reverse clip) is hit.
    assert!(
        results.is_empty(),
        "Non-intersecting curves should yield no results, got: {results:?}"
    );
}
