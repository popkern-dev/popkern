use super::*;
use crate::geometry::curve::CurveEval;

/// Build a rectangular boundary from four line segments.
fn rectangular_boundary(x0: f64, y0: f64, x1: f64, y1: f64) -> CoonsBoundary {
    let bottom = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x0, y0, 0.0), Point3::new(x1, y0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let top = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x0, y1, 0.0), Point3::new(x1, y1, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let left = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x0, y0, 0.0), Point3::new(x0, y1, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let right = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(x1, y0, 0.0), Point3::new(x1, y1, 0.0)],
        weights: vec![1.0, 1.0],
    };
    CoonsBoundary {
        bottom,
        top,
        left,
        right,
    }
}

#[test]
fn rectangular_coons_patch_has_small_boundary_deviation() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let result = construct_coons_patch(&boundary, 8, 8, 1e-6);
    assert!(result.is_some(), "construction should succeed");
    let r = result.as_ref();
    // The NURBS approximation uses sampled control points (not true
    // interpolation), so deviation is non-zero but should be small
    // for a planar bilinear case.
    assert!(
        r.map(|r| r.max_boundary_deviation).unwrap_or(f64::MAX) < 0.5,
        "boundary deviation {} is too large",
        r.map(|r| r.max_boundary_deviation).unwrap_or(f64::MAX),
    );
}

#[test]
fn rectangular_corners_are_correct() {
    let boundary = rectangular_boundary(0.0, 0.0, 2.0, 3.0);
    let result = construct_coons_patch(&boundary, 8, 8, 1e-6);
    assert!(result.is_some(), "construction should succeed");
    let r = result.as_ref();
    let corners = r.map(|r| r.corners);
    assert!(corners.is_some());
    let c = corners.as_ref();
    // P_00 = (0,0,0), P_10 = (2,0,0), P_01 = (0,3,0), P_11 = (2,3,0)
    assert!(
        c.map(|c| c[0].distance_to(&Point3::new(0.0, 0.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
    assert!(
        c.map(|c| c[1].distance_to(&Point3::new(2.0, 0.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
    assert!(
        c.map(|c| c[2].distance_to(&Point3::new(0.0, 3.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
    assert!(
        c.map(|c| c[3].distance_to(&Point3::new(2.0, 3.0, 0.0)))
            .unwrap_or(1.0)
            < 1e-10
    );
}

#[test]
fn coons_patch_from_four_line_segments() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let result = construct_coons_patch(&boundary, 8, 8, 1e-6);
    assert!(result.is_some());
    let r = result.as_ref();
    assert!(r.map(|r| r.surface.num_u).unwrap_or(0) == 8);
    assert!(r.map(|r| r.surface.num_v).unwrap_or(0) == 8);
}

#[test]
fn coons_patch_non_closing_boundary_returns_none() {
    // Left curve starts at a different point than bottom curve starts.
    let bottom = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let top = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let left = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(5.0, 5.0, 5.0), // wrong start
            Point3::new(0.0, 1.0, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };
    let right = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let boundary = CoonsBoundary {
        bottom,
        top,
        left,
        right,
    };
    let result = construct_coons_patch(&boundary, 8, 8, 1e-6);
    assert!(result.is_none(), "non-closing boundary should return None");
}

#[test]
fn coons_patch_with_curved_boundaries() {
    // Four quadratic arcs forming a lens-like shape
    let bottom = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, -0.2, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let top = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.5, 1.2, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let left = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(-0.2, 0.5, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let right = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.2, 0.5, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };

    let boundary = CoonsBoundary {
        bottom,
        top,
        left,
        right,
    };
    let result = construct_coons_patch(&boundary, 10, 10, 1e-6);
    assert!(
        result.is_some(),
        "curved boundary construction should succeed"
    );
    let r = result.as_ref();
    // Boundary deviation should be reasonably small
    assert!(
        r.map(|r| r.max_boundary_deviation).unwrap_or(f64::MAX) < 1.0,
        "curved boundary deviation {} is too large",
        r.map(|r| r.max_boundary_deviation).unwrap_or(f64::MAX),
    );
}

#[test]
fn coons_patch_too_small_grid_returns_none() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    // Grid too small for degree 3
    assert!(construct_coons_patch(&boundary, 3, 8, 1e-6).is_none());
    assert!(construct_coons_patch(&boundary, 8, 3, 1e-6).is_none());
    assert!(construct_coons_patch(&boundary, 2, 2, 1e-6).is_none());
}

#[test]
fn coons_eval_at_corners_matches_corner_points() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];

    let p00 = vec_to_point(coons_eval(
        &boundary.bottom,
        &boundary.top,
        &boundary.left,
        &boundary.right,
        &corners,
        0.0,
        0.0,
    ));
    let p10 = vec_to_point(coons_eval(
        &boundary.bottom,
        &boundary.top,
        &boundary.left,
        &boundary.right,
        &corners,
        1.0,
        0.0,
    ));
    let p01 = vec_to_point(coons_eval(
        &boundary.bottom,
        &boundary.top,
        &boundary.left,
        &boundary.right,
        &corners,
        0.0,
        1.0,
    ));
    let p11 = vec_to_point(coons_eval(
        &boundary.bottom,
        &boundary.top,
        &boundary.left,
        &boundary.right,
        &corners,
        1.0,
        1.0,
    ));

    assert!(p00.distance_to(&corners[0]) < 1e-10, "P_00 mismatch");
    assert!(p10.distance_to(&corners[1]) < 1e-10, "P_10 mismatch");
    assert!(p01.distance_to(&corners[2]) < 1e-10, "P_01 mismatch");
    assert!(p11.distance_to(&corners[3]) < 1e-10, "P_11 mismatch");
}

#[test]
fn coons_eval_along_bottom_boundary_matches_curve() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];

    // Along v=0, the Coons patch should match the bottom curve
    for i in 0..=10 {
        let u = i as f64 / 10.0;
        let patch_pt = vec_to_point(coons_eval(
            &boundary.bottom,
            &boundary.top,
            &boundary.left,
            &boundary.right,
            &corners,
            u,
            0.0,
        ));
        let (t0, t1) = boundary.bottom.parameter_range();
        let curve_pt = boundary.bottom.point_at(t0 + u * (t1 - t0));
        assert!(
            patch_pt.distance_to(&curve_pt) < 1e-10,
            "bottom boundary mismatch at u={u}: patch={patch_pt:?} curve={curve_pt:?}",
        );
    }
}

#[test]
fn coons_eval_along_left_boundary_matches_curve() {
    let boundary = rectangular_boundary(0.0, 0.0, 1.0, 1.0);
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];

    // Along u=0, the Coons patch should match the left curve
    for i in 0..=10 {
        let v = i as f64 / 10.0;
        let patch_pt = vec_to_point(coons_eval(
            &boundary.bottom,
            &boundary.top,
            &boundary.left,
            &boundary.right,
            &corners,
            0.0,
            v,
        ));
        let (t0, t1) = boundary.left.parameter_range();
        let curve_pt = boundary.left.point_at(t0 + v * (t1 - t0));
        assert!(
            patch_pt.distance_to(&curve_pt) < 1e-10,
            "left boundary mismatch at v={v}: patch={patch_pt:?} curve={curve_pt:?}",
        );
    }
}

#[test]
fn fill_four_sided_hole_rejects_non_four_sided() {
    let brep = BRep::new(crate::tolerance::Tolerance::default());
    let hole = DetectedHole {
        boundary_edges: vec![],
        num_sides: 3,
    };
    assert!(fill_four_sided_hole(&brep, &hole, 1e-6).is_none());
}

#[test]
fn build_clamped_knots_degree_3_basic() {
    let knots = build_clamped_knots(8, 3, 0.0, 1.0);
    // Should have 8 + 3 + 1 = 12 knots
    assert_eq!(knots.len(), 12);
    // First 4 should be 0.0
    for &k in &knots[..4] {
        assert!((k - 0.0).abs() < 1e-14);
    }
    // Last 4 should be 1.0
    for &k in &knots[8..] {
        assert!((k - 1.0).abs() < 1e-14);
    }
    // Knots should be non-decreasing
    for w in knots.windows(2) {
        assert!(w[1] >= w[0], "knots not non-decreasing: {:?}", knots);
    }
}

#[test]
fn build_clamped_knots_minimum_control_points() {
    // Minimum for degree 3: 4 control points, no interior knots
    let knots = build_clamped_knots(4, 3, 0.0, 1.0);
    assert_eq!(knots.len(), 8); // 4 + 3 + 1
    assert_eq!(knots, vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
}
