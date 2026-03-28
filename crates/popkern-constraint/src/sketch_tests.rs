//! Tests for the Sketch API: element creation, constraints, driven values,
//! drag solving, and element deletion.

use crate::constraint::Constraint;
use crate::error::ConstraintError;
use crate::Sketch;

// ── Element API ──────────────────────────────────────────────────────────

#[test]
fn add_point_and_read_position() {
    let mut sketch = Sketch::new();
    let h = sketch.add_point(3.0, 4.0);
    assert_eq!(sketch.point_position(h), Some((3.0, 4.0)));
}

#[test]
fn add_line_shares_point_params() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let _l = sketch.add_line(a, b);

    // Mutate the underlying param vector directly to verify sharing.
    let idx = sketch.points[a].param_indices()[0];
    sketch.params[idx] = 5.0;

    assert_eq!(sketch.point_position(a), Some((5.0, 0.0)));
}

#[test]
fn add_and_remove_constraint() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let ch = sketch.add_constraint(Constraint::Coincident(a, b));
    assert!(sketch.constraint_status(ch).is_some());

    sketch.remove_constraint(ch);
    assert!(sketch.constraint_status(ch).is_none());
}

// ── Driven values & set_constraint_value ─────────────────────────────────

#[test]
fn driven_distance_reports_current_distance() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 4.0);
    let ch = sketch.add_driven(Constraint::Distance(a, b, 0.0));

    let val = sketch.driven_value(ch);
    assert!(val.is_some());
    assert!(
        (val.unwrap() - 5.0).abs() < 1e-8,
        "expected distance 5.0, got {:?}",
        val
    );
}

#[test]
fn driving_distance_set_value_and_solve() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(a));
    let ch = sketch.add_constraint(Constraint::Distance(a, b, 3.0));

    sketch.solve().unwrap();
    let d1 = sketch.driven_value(ch).unwrap();
    assert!((d1 - 3.0).abs() < 1e-8, "initial distance: {}", d1);

    sketch.set_constraint_value(ch, 5.0).unwrap();
    sketch.solve().unwrap();
    let d2 = sketch.driven_value(ch).unwrap();
    assert!((d2 - 5.0).abs() < 1e-8, "updated distance: {}", d2);
}

#[test]
fn set_value_on_horizontal_returns_error() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let l = sketch.add_line(a, b);
    let ch = sketch.add_constraint(Constraint::Horizontal(l));

    let result = sketch.set_constraint_value(ch, 1.0);
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ConstraintError::NotDimensional
    ));
}

#[test]
fn driven_angle_reports_current_angle() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let c = sketch.add_point(0.0, 0.0);
    let d = sketch.add_point(0.0, 1.0);
    let l1 = sketch.add_line(a, b);
    let l2 = sketch.add_line(c, d);
    let ch = sketch.add_driven(Constraint::Angle(l1, l2, 0.0));

    let val = sketch.driven_value(ch).unwrap();
    let expected = std::f64::consts::FRAC_PI_2;
    assert!(
        (val - expected).abs() < 1e-8,
        "expected angle pi/2, got {}",
        val
    );
}

// ── Drag solving ─────────────────────────────────────────────────────────

#[test]
fn drag_horizontal_line_stays_horizontal() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.0);
    let l = sketch.add_line(p1, p2);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(l));

    // Drag p2 to (5.0, 2.0) — the temporary Fixed constraint at (5,2)
    // conflicts with Horizontal, producing an over-constrained LS
    // compromise. After the temp constraint is removed, the geometry
    // reflects the solver's best-fit: p2.x near 5, p2.y pulled toward 0
    // by Horizontal but offset by the LS compromise with Fixed(5,2).
    sketch.begin_drag(p2).unwrap();
    let _result = sketch.drag_to(5.0, 2.0).unwrap();
    sketch.end_drag();

    let (x2, _y2) = sketch.point_position(p2).unwrap();
    // x should have moved towards the target.
    assert!(
        (x2 - 5.0).abs() < 1.0,
        "p2.x should be near 5.0, got {}",
        x2
    );

    // After removing the temp constraint and re-solving, we can verify the
    // horizontal constraint is truly satisfied.
    sketch.solve().unwrap();
    let (x2_final, y2_final) = sketch.point_position(p2).unwrap();
    // With only Fixed(p1) and Horizontal, the system is under-constrained;
    // the solver keeps p2 near its last position.
    assert!(
        y2_final.abs() < 1e-6,
        "p2.y should snap to 0, got {}",
        y2_final
    );
    assert!(
        (x2_final - 5.0).abs() < 1.0,
        "p2.x should stay near 5.0, got {}",
        x2_final
    );
}

#[test]
fn drag_to_without_begin_returns_error() {
    let mut sketch = Sketch::new();
    let _p = sketch.add_point(0.0, 0.0);

    let result = sketch.drag_to(1.0, 1.0);
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), ConstraintError::NoDragActive));
}

#[test]
fn begin_drag_invalid_handle_returns_error() {
    let mut sketch = Sketch::new();
    let p = sketch.add_point(0.0, 0.0);
    // Remove the point's handle by creating a fresh sketch (handle is stale).
    let mut other = Sketch::new();
    let result = other.begin_drag(p);
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        ConstraintError::InvalidHandle
    ));
}

// ── Element deletion with cascade ────────────────────────────────────────

#[test]
fn remove_point_cascades_line_and_constraints() {
    use crate::elements::ElementHandle;

    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(1.0, 0.0);
    let p3 = sketch.add_point(2.0, 0.0);
    let l12 = sketch.add_line(p1, p2);
    let _l23 = sketch.add_line(p2, p3);

    let ch_horiz = sketch.add_constraint(Constraint::Horizontal(l12));
    let ch_fix = sketch.add_constraint(Constraint::Fixed(p1));

    // Remove p1 → l12 should be removed (cascade), ch_horiz on l12 too,
    // ch_fix on p1 too. p2 and p3 survive (not endpoints of removed line
    // from their perspective, and p2/p3 are not the removed point).
    sketch.remove_element(ElementHandle::Point(p1));

    // p1 is gone.
    assert!(sketch.point_position(p1).is_none());
    // l12 is gone (p1 was an endpoint).
    assert!(sketch.lines.get(l12).is_none());
    // Constraints on l12 and p1 are gone.
    assert!(sketch.constraint_status(ch_horiz).is_none());
    assert!(sketch.constraint_status(ch_fix).is_none());
    // p2 and p3 survive.
    assert!(sketch.point_position(p2).is_some());
    assert!(sketch.point_position(p3).is_some());
}

#[test]
fn remove_line_keeps_endpoint_points() {
    use crate::elements::ElementHandle;

    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(1.0, 0.0);
    let l = sketch.add_line(p1, p2);

    let ch = sketch.add_constraint(Constraint::Horizontal(l));

    sketch.remove_element(ElementHandle::Line(l));

    // Line and its constraint are gone.
    assert!(sketch.lines.get(l).is_none());
    assert!(sketch.constraint_status(ch).is_none());
    // But the endpoint points survive.
    assert!(sketch.point_position(p1).is_some());
    assert!(sketch.point_position(p2).is_some());
}

#[test]
fn remove_arc_cascades_constraints() {
    use crate::elements::ElementHandle;

    let mut sketch = Sketch::new();
    let a = sketch.add_arc(0.0, 0.0, 1.0, 0.0, std::f64::consts::PI);
    let p = sketch.add_point(1.0, 0.0);

    let ch = sketch.add_constraint(Constraint::PointOnArc(p, a));

    sketch.remove_element(ElementHandle::Arc(a));

    assert!(sketch.arcs.get(a).is_none());
    assert!(sketch.constraint_status(ch).is_none());
    // The point survives.
    assert!(sketch.point_position(p).is_some());
}
