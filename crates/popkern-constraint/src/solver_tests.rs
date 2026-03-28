//! Convergence tests for the Newton–Raphson constraint solver.

use crate::constraint::Constraint;
use crate::error::SolveStatus;
use crate::solver::SolverOptions;
use crate::Sketch;

/// Helper: distance between two parameter-space points.
fn dist(sketch: &Sketch, a: crate::PointHandle, b: crate::PointHandle) -> f64 {
    let (ax, ay) = sketch.point_position(a).unwrap();
    let (bx, by) = sketch.point_position(b).unwrap();
    ((bx - ax).powi(2) + (by - ay).powi(2)).sqrt()
}

#[test]
fn horizontal_line_converges() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 2.0);
    let l = sketch.add_line(p1, p2);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(l));

    let result = sketch
        .solve_with_options(&SolverOptions::default())
        .unwrap();
    assert!(
        result.residual_norm < 1e-10,
        "residual too large: {}",
        result.residual_norm
    );

    let (_, y2) = sketch.point_position(p2).unwrap();
    assert!(
        y2.abs() < 1e-8,
        "p2.y should be 0 for horizontal line, got {}",
        y2
    );

    // p1 should remain at origin.
    let (x1, y1) = sketch.point_position(p1).unwrap();
    assert!(x1.abs() < 1e-8, "p1.x moved: {}", x1);
    assert!(y1.abs() < 1e-8, "p1.y moved: {}", y1);
}

#[test]
fn fixed_distance_converges() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(1.0, 1.0);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Distance(p1, p2, 2.0));

    let result = sketch
        .solve_with_options(&SolverOptions::default())
        .unwrap();
    assert!(
        result.residual_norm < 1e-10,
        "residual too large: {}",
        result.residual_norm
    );

    let d = dist(&sketch, p1, p2);
    assert!((d - 2.0).abs() < 1e-8, "distance should be 2.0, got {}", d);
}

#[test]
fn coincident_points_converge() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(1.0, 2.0);
    let p2 = sketch.add_point(4.0, 6.0);

    sketch.add_constraint(Constraint::Coincident(p1, p2));

    let result = sketch
        .solve_with_options(&SolverOptions::default())
        .unwrap();
    assert!(
        result.residual_norm < 1e-10,
        "residual too large: {}",
        result.residual_norm
    );

    let (x1, y1) = sketch.point_position(p1).unwrap();
    let (x2, y2) = sketch.point_position(p2).unwrap();
    assert!(
        (x1 - x2).abs() < 1e-8 && (y1 - y2).abs() < 1e-8,
        "points should coincide: ({}, {}) vs ({}, {})",
        x1,
        y1,
        x2,
        y2
    );
}

#[test]
fn rectangle_converges() {
    // Build a rectangle: 4 points, 4 lines.
    // Constraints: fix p1 at origin, horizontal bottom and top,
    // vertical left and right, width=3, height=2.
    let mut sketch = Sketch::new();

    // Initial positions (slightly off to test convergence).
    let p1 = sketch.add_point(0.0, 0.0); // bottom-left
    let p2 = sketch.add_point(2.5, 0.3); // bottom-right
    let p3 = sketch.add_point(2.8, 1.7); // top-right
    let p4 = sketch.add_point(0.2, 2.1); // top-left

    let bottom = sketch.add_line(p1, p2);
    let right = sketch.add_line(p2, p3);
    let top = sketch.add_line(p3, p4);
    let left = sketch.add_line(p4, p1);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(bottom));
    sketch.add_constraint(Constraint::Horizontal(top));
    sketch.add_constraint(Constraint::Vertical(right));
    sketch.add_constraint(Constraint::Vertical(left));
    sketch.add_constraint(Constraint::Distance(p1, p2, 3.0)); // width
    sketch.add_constraint(Constraint::Distance(p1, p4, 2.0)); // height

    let result = sketch
        .solve_with_options(&SolverOptions::default())
        .unwrap();
    assert!(
        result.residual_norm < 1e-10,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify the rectangle geometry.
    let (x1, y1) = sketch.point_position(p1).unwrap();
    let (x2, y2) = sketch.point_position(p2).unwrap();
    let (x3, y3) = sketch.point_position(p3).unwrap();
    let (x4, y4) = sketch.point_position(p4).unwrap();

    // p1 at origin.
    assert!(x1.abs() < 1e-8 && y1.abs() < 1e-8, "p1 not at origin");

    // Bottom horizontal: y1 == y2.
    assert!((y1 - y2).abs() < 1e-8, "bottom not horizontal");

    // Right vertical: x2 == x3.
    assert!((x2 - x3).abs() < 1e-8, "right not vertical");

    // Top horizontal: y3 == y4.
    assert!((y3 - y4).abs() < 1e-8, "top not horizontal");

    // Left vertical: x4 == x1.
    assert!((x4 - x1).abs() < 1e-8, "left not vertical");

    // Width = 3.
    assert!(
        (dist(&sketch, p1, p2) - 3.0).abs() < 1e-8,
        "width wrong: {}",
        dist(&sketch, p1, p2)
    );

    // Height = 2.
    assert!(
        (dist(&sketch, p1, p4) - 2.0).abs() < 1e-8,
        "height wrong: {}",
        dist(&sketch, p1, p4)
    );
}

#[test]
fn empty_sketch_solves_immediately() {
    let mut sketch = Sketch::new();
    let result = sketch.solve().unwrap();
    assert_eq!(result.iterations, 0);
    assert_eq!(result.residual_norm, 0.0);
}

#[test]
fn already_satisfied_returns_immediately() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.0);
    let l = sketch.add_line(p1, p2);
    sketch.add_constraint(Constraint::Horizontal(l));

    let result = sketch.solve().unwrap();
    assert!(result.residual_norm < 1e-10);
    assert_eq!(result.iterations, 0);
}

#[test]
fn perpendicular_lines_converge() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(2.0, 1.0); // slightly off-axis
    let p3 = sketch.add_point(1.0, 3.0);
    let l1 = sketch.add_line(p1, p2);
    let l2 = sketch.add_line(p1, p3);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Perpendicular(l1, l2));

    let result = sketch.solve().unwrap();
    assert!(
        result.residual_norm < 1e-10,
        "residual: {}",
        result.residual_norm
    );

    // Verify dot product of direction vectors is zero.
    let (x1, y1) = sketch.point_position(p1).unwrap();
    let (x2, y2) = sketch.point_position(p2).unwrap();
    let (x3, y3) = sketch.point_position(p3).unwrap();
    let dot = (x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1);
    assert!(dot.abs() < 1e-8, "not perpendicular: dot = {}", dot);
}

#[test]
fn solve_status_reports_fully_constrained() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.5);
    let l = sketch.add_line(p1, p2);

    // 4 params, 4 equations (Fixed=2, Horizontal=1, Distance=1) — fully constrained.
    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(l));
    sketch.add_constraint(Constraint::Distance(p1, p2, 3.0));

    let result = sketch.solve().unwrap();
    assert!(result.residual_norm < 1e-10);
    // 2 (Fixed) + 1 (Horizontal) + 1 (Distance) = 4 equations, 4 params.
    assert!(matches!(result.status, SolveStatus::FullyConstrained));
}

#[test]
fn under_constrained_status() {
    let mut sketch = Sketch::new();
    let _p1 = sketch.add_point(0.0, 0.0);
    let _p2 = sketch.add_point(1.0, 0.0);
    // No constraints at all — under-constrained.
    let result = sketch.solve().unwrap();
    assert!(matches!(result.status, SolveStatus::UnderConstrained(_)));
}
