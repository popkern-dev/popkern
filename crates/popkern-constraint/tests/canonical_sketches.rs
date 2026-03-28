//! Integration tests for canonical sketch scenarios: fully constrained
//! rectangles, tangent arcs, under/over-constrained systems, drag
//! solving, and parallel lines.

use popkern_constraint::{Constraint, Sketch, SolveStatus};
use std::f64::consts::FRAC_PI_2;

/// Tolerance for floating-point comparisons.
const TOL: f64 = 1e-6;

/// Assert that two floats are within [`TOL`] of each other.
fn approx(a: f64, b: f64, label: &str) {
    assert!(
        (a - b).abs() < TOL,
        "{}: expected {}, got {} (delta {})",
        label,
        b,
        a,
        (a - b).abs()
    );
}

// ── 1. Fully constrained rectangle ─────────────────────────────────────

#[test]
fn fully_constrained_rectangle() {
    let mut sketch = Sketch::new();

    // Four corner points (initial guess: rough rectangle).
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(4.0, 0.1);
    let p2 = sketch.add_point(4.1, 2.9);
    let p3 = sketch.add_point(-0.1, 3.1);

    // Four sides.
    let bottom = sketch.add_line(p0, p1);
    let right = sketch.add_line(p1, p2);
    let top = sketch.add_line(p2, p3);
    let left = sketch.add_line(p3, p0);

    // Fix origin.
    sketch.add_constraint(Constraint::Fixed(p0));

    // Horizontal bottom and top, vertical sides.
    sketch.add_constraint(Constraint::Horizontal(bottom));
    sketch.add_constraint(Constraint::Horizontal(top));
    sketch.add_constraint(Constraint::Vertical(right));
    sketch.add_constraint(Constraint::Vertical(left));

    // Width = 5, Height = 3.
    sketch.add_constraint(Constraint::Distance(p0, p1, 5.0));
    sketch.add_constraint(Constraint::Distance(p0, p3, 3.0));

    let result = sketch.solve().expect("solve");
    assert!(
        matches!(result.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained, got {:?}",
        result.status
    );

    // Verify corner positions.
    let (x0, y0) = sketch.point_position(p0).unwrap();
    let (x1, y1) = sketch.point_position(p1).unwrap();
    let (x2, y2) = sketch.point_position(p2).unwrap();
    let (x3, y3) = sketch.point_position(p3).unwrap();

    approx(x0, 0.0, "p0.x");
    approx(y0, 0.0, "p0.y");
    approx(x1, 5.0, "p1.x");
    approx(y1, 0.0, "p1.y");
    approx(x2, 5.0, "p2.x");
    approx(y2, 3.0, "p2.y");
    approx(x3, 0.0, "p3.x");
    approx(y3, 3.0, "p3.y");
}

// ── 2. Tangent arc between lines ────────────────────────────────────────

#[test]
fn tangent_arc_between_lines() {
    let mut sketch = Sketch::new();

    // Horizontal line along x-axis.
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(5.0, 0.0);
    let h_line = sketch.add_line(a, b);

    // Vertical line along y-axis.
    let c = sketch.add_point(0.0, 0.0);
    let d = sketch.add_point(0.0, 5.0);
    let v_line = sketch.add_line(c, d);

    // Arc in the first quadrant, tangent to both lines.
    let arc = sketch.add_arc(2.0, 2.0, 2.0, 0.0, FRAC_PI_2);

    // Fix the lines in place.
    sketch.add_constraint(Constraint::Fixed(a));
    sketch.add_constraint(Constraint::Fixed(b));
    sketch.add_constraint(Constraint::Fixed(c));
    sketch.add_constraint(Constraint::Fixed(d));

    // Tangent constraints.
    sketch.add_constraint(Constraint::Tangent(h_line, arc));
    sketch.add_constraint(Constraint::Tangent(v_line, arc));

    let result = sketch.solve().expect("solve");

    // Verify the arc centre is at (r, r) where r is the radius.
    // Use arc_params() accessor.
    let (cx, cy, r, _, _) = sketch.arc_params(arc).expect("arc params");

    // A circle tangent to x-axis and y-axis in the first quadrant has
    // centre at (r, r).
    approx(cx, r, "arc cx == r");
    approx(cy, r, "arc cy == r");
    assert!(
        r > 0.5,
        "radius should be positive and reasonable, got {}",
        r
    );
    assert!(
        result.residual_norm < 1e-8,
        "residual: {}",
        result.residual_norm
    );
}

// ── 3. Under-constrained sketch ─────────────────────────────────────────

#[test]
fn under_constrained_sketch() {
    let mut sketch = Sketch::new();

    // Two free points, one line, horizontal constraint only.
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 1.0);
    let l = sketch.add_line(p1, p2);
    sketch.add_constraint(Constraint::Horizontal(l));

    let result = sketch.solve().expect("solve");

    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            // 4 DOF (2 points * 2 coords) minus 1 (horizontal) = 3 free.
            // The DOF analyser may merge co-directional freedoms, so we
            // just check that there is at least one free element.
            assert!(!freedoms.is_empty(), "expected at least one free element");
            let total_dof: usize = freedoms.iter().map(|f| f.directions.len()).sum();
            assert!(total_dof >= 2, "expected at least 2 DOF, got {}", total_dof);
        }
        other => panic!("expected UnderConstrained, got {:?}", other),
    }
}

// ── 4. Over-constrained sketch ──────────────────────────────────────────

#[test]
fn over_constrained_sketch() {
    let mut sketch = Sketch::new();

    // A line with horizontal + vertical constraints (conflicting unless
    // the line is a point). The solver should detect redundancy.
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.0);
    let l = sketch.add_line(p1, p2);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(l));
    sketch.add_constraint(Constraint::Vertical(l));
    sketch.add_constraint(Constraint::Distance(p1, p2, 3.0));

    let result = sketch.solve().expect("solve");

    // The system is over-constrained because horizontal + vertical +
    // non-zero distance conflict. We accept either OverConstrained or
    // Failed as valid -- both indicate the conflict was detected.
    match &result.status {
        SolveStatus::OverConstrained(handles) => {
            assert!(
                !handles.is_empty(),
                "expected conflicting handles, got empty"
            );
        }
        SolveStatus::Failed(_) => {
            // Solver could not satisfy all constraints -- this is also a
            // valid detection of over-constraint.
        }
        other => panic!("expected OverConstrained or Failed, got {:?}", other),
    }
}

// ── 5. Drag test ────────────────────────────────────────────────────────

#[test]
fn drag_maintains_horizontal_constraint() {
    let mut sketch = Sketch::new();

    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.0);
    let l = sketch.add_line(p1, p2);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(l));

    // First solve to establish consistent state.
    sketch.solve().expect("initial solve");

    // Drag p2 to (7.0, 3.0).
    sketch.begin_drag(p2).expect("begin_drag");
    sketch.drag_to(7.0, 3.0).expect("drag_to");
    sketch.end_drag();

    // Re-solve to snap to constraints.
    sketch.solve().expect("re-solve after drag");

    let (_x2, y2) = sketch.point_position(p2).unwrap();
    // Horizontal constraint means p2.y must equal p1.y (0.0).
    assert!(
        y2.abs() < TOL,
        "horizontal constraint violated: p2.y = {}",
        y2
    );
}

// ── 6. Parallel lines ───────────────────────────────────────────────────

#[test]
fn parallel_lines_with_distance() {
    let mut sketch = Sketch::new();

    // Line 1: roughly along x-axis.
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(4.0, 0.1);
    let l1 = sketch.add_line(a, b);

    // Line 2: offset above, roughly parallel.
    let c = sketch.add_point(0.0, 3.1);
    let d = sketch.add_point(4.0, 2.9);
    let l2 = sketch.add_line(c, d);

    // Fix line 1 endpoints.
    sketch.add_constraint(Constraint::Fixed(a));
    sketch.add_constraint(Constraint::Horizontal(l1));
    sketch.add_constraint(Constraint::Distance(a, b, 4.0));

    // Line 2 is parallel and 3 units away.
    sketch.add_constraint(Constraint::Parallel(l1, l2));
    sketch.add_constraint(Constraint::PointLineDistance(c, l1, 3.0));
    sketch.add_constraint(Constraint::PointLineDistance(d, l1, 3.0));

    let result = sketch.solve().expect("solve");
    assert!(
        result.residual_norm < 1e-8,
        "residual: {}",
        result.residual_norm
    );

    // Check parallelism: both lines should have the same direction.
    let (xa, ya) = sketch.point_position(a).unwrap();
    let (xb, yb) = sketch.point_position(b).unwrap();
    let (xc, yc) = sketch.point_position(c).unwrap();
    let (xd, yd) = sketch.point_position(d).unwrap();

    let dx1 = xb - xa;
    let dy1 = yb - ya;
    let dx2 = xd - xc;
    let dy2 = yd - yc;

    // Cross product should be zero for parallel lines.
    let cross = dx1 * dy2 - dy1 * dx2;
    assert!(
        cross.abs() < TOL,
        "lines are not parallel: cross product = {}",
        cross
    );

    // Both should be 3 units apart vertically.
    approx(yc, 3.0, "p_c.y");
    approx(yd, 3.0, "p_d.y");
}

// ── 7. Undo/redo round-trip in integration context ──────────────────────

#[test]
fn undo_redo_preserves_solve_consistency() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(p1));
    let ch = sketch.add_constraint(Constraint::Distance(p1, p2, 3.0));

    sketch.solve().expect("initial solve");
    let (x_before, _) = sketch.point_position(p2).unwrap();
    approx(x_before, 3.0, "before set_value");

    // Change distance to 5.
    sketch.set_constraint_value(ch, 5.0).expect("set_value");
    sketch.solve().expect("solve after set");
    let (x_after, _) = sketch.point_position(p2).unwrap();
    approx(x_after, 5.0, "after set_value");

    // Undo the set_constraint_value.
    assert!(sketch.undo());
    sketch.solve().expect("solve after undo");
    let (x_undo, _) = sketch.point_position(p2).unwrap();
    approx(x_undo, 3.0, "after undo");

    // Redo the set_constraint_value.
    assert!(sketch.redo());
    sketch.solve().expect("solve after redo");
    let (x_redo, _) = sketch.point_position(p2).unwrap();
    approx(x_redo, 5.0, "after redo");
}

// ── 8. Serialize round-trip in integration context ──────────────────────

#[test]
fn serialize_roundtrip_integration() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(0.0, 0.0);
    let p1 = sketch.add_point(5.0, 0.0);
    let p2 = sketch.add_point(5.0, 3.0);
    let p3 = sketch.add_point(0.0, 3.0);
    let b = sketch.add_line(p0, p1);
    let r = sketch.add_line(p1, p2);
    let t = sketch.add_line(p2, p3);
    let l = sketch.add_line(p3, p0);
    sketch.add_constraint(Constraint::Fixed(p0));
    sketch.add_constraint(Constraint::Horizontal(b));
    sketch.add_constraint(Constraint::Vertical(r));
    sketch.add_constraint(Constraint::Horizontal(t));
    sketch.add_constraint(Constraint::Vertical(l));
    sketch.add_constraint(Constraint::Distance(p0, p1, 5.0));
    sketch.add_constraint(Constraint::Distance(p0, p3, 3.0));

    sketch.solve().expect("solve original");

    let json = serde_json::to_string(&sketch).expect("serialize");
    let mut loaded: Sketch = serde_json::from_str(&json).expect("deserialize");

    let result = loaded.solve().expect("solve loaded");
    assert!(matches!(result.status, SolveStatus::FullyConstrained));

    let (x2, y2) = loaded.point_position(p2).unwrap();
    approx(x2, 5.0, "loaded p2.x");
    approx(y2, 3.0, "loaded p2.y");
}
