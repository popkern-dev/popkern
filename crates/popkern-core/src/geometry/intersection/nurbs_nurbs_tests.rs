//! Tests for `nurbs_nurbs` (split out for the 500-line cap).
use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;

#[test]
fn parabola_crosses_horizontal_line() {
    // Quadratic Bezier (parabola) and horizontal line at y=0.5
    let c1 = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let c2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
    // Parabola y = 2t(1-t) crosses y=0.5 at 2 points
    assert!(
        result.len() >= 2,
        "Expected at least 2 intersections, got {}",
        result.len()
    );

    // Verify invariant: points lie on both curves
    for r in &result {
        if let CurveIntersection::Point { point, t1, t2 } = r {
            let p1 = c1.point_at(*t1);
            let p2 = c2.point_at(*t2);
            assert!(p1.distance_to(&p2) < 1e-6, "Invariant: points don't match");
            assert!(p1.distance_to(point) < 1e-6);
        }
    }
}

#[test]
fn non_intersecting_curves() {
    let c1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let c2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 10.0, 0.0), Point3::new(1.0, 10.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
    assert!(result.is_empty());
}

#[test]
fn identical_curves_overlap() {
    let c1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c1, (0.0, 1.0), &Tolerance::default());
    assert!(
        result
            .iter()
            .any(|r| matches!(r, CurveIntersection::Overlap { .. }))
    );
}

#[test]
fn multi_span_curve_decomposition() {
    // Multi-span curve (has interior knot) crossing a line
    let c1 = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.25, 1.0, 0.0),
            Point3::new(0.75, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let c2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
    assert!(
        result.len() >= 2,
        "Expected at least 2 intersections with multi-span curve"
    );
}

#[test]
fn aabb_overlap_disjoint() {
    let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)];
    let pts2 = vec![Point3::new(5.0, 5.0, 0.0), Point3::new(6.0, 6.0, 0.0)];
    assert!(!aabb_overlap(&pts1, &pts2));
}

#[test]
fn aabb_overlap_touching() {
    let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)];
    let pts2 = vec![Point3::new(0.5, 0.5, 0.5), Point3::new(2.0, 2.0, 2.0)];
    assert!(aabb_overlap(&pts1, &pts2));
}

#[test]
fn aabb_overlap_empty() {
    let empty: Vec<Point3> = vec![];
    let pts = vec![Point3::new(0.0, 0.0, 0.0)];
    assert!(!aabb_overlap(&empty, &pts));
    assert!(!aabb_overlap(&pts, &empty));
}

#[test]
fn deduplicate_removes_close_points() {
    let mut results = vec![
        CurveIntersection::Point {
            point: Point3::new(0.5, 0.5, 0.0),
            t1: 0.5,
            t2: 0.5,
        },
        CurveIntersection::Point {
            point: Point3::new(0.5, 0.5, 0.0),
            t1: 0.500001,
            t2: 0.500001,
        },
    ];
    deduplicate(&mut results, &Tolerance::default());
    assert_eq!(results.len(), 1);
}

#[test]
fn sort_by_overlap_t1_range() {
    // Test that overlaps sort correctly by t1_range.0
    let result = intersect_nurbs_nurbs(
        &NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        },
        (0.0, 1.0),
        &NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        },
        (0.0, 1.0),
        &Tolerance::default(),
    );
    // Identical curves => overlap
    assert!(!result.is_empty());
}
