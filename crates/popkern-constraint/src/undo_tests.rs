//! Tests for the undo/redo system.

use crate::constraint::Constraint;
use crate::elements::ElementHandle;
use crate::Sketch;

#[test]
fn undo_add_point_removes_point_and_shrinks_params() {
    let mut sketch = Sketch::new();
    let p = sketch.add_point(3.0, 4.0);

    assert_eq!(sketch.points.len(), 1);
    assert_eq!(sketch.params.len(), 2);

    assert!(sketch.undo());
    assert!(sketch.point_position(p).is_none());
    assert_eq!(sketch.points.len(), 0);
    assert_eq!(sketch.params.len(), 0);
}

#[test]
fn redo_add_point_restores_point() {
    let mut sketch = Sketch::new();
    let _p = sketch.add_point(3.0, 4.0);

    sketch.undo();
    assert_eq!(sketch.points.len(), 0);

    assert!(sketch.redo());
    assert_eq!(sketch.points.len(), 1);
    assert_eq!(sketch.params.len(), 2);
    assert_eq!(sketch.params[0], 3.0);
    assert_eq!(sketch.params[1], 4.0);
}

#[test]
fn undo_add_constraint_removes_it() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let ch = sketch.add_constraint(Constraint::Coincident(a, b));

    assert!(sketch.constraint_status(ch).is_some());
    assert!(sketch.undo()); // undo add_constraint
    assert!(sketch.constraint_status(ch).is_none());
    assert_eq!(sketch.constraints.len(), 0);
}

#[test]
fn undo_remove_constraint_restores_it() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let ch = sketch.add_constraint(Constraint::Coincident(a, b));

    sketch.remove_constraint(ch);
    assert!(sketch.constraint_status(ch).is_none());

    assert!(sketch.undo()); // undo remove_constraint
    assert!(sketch.constraint_status(ch).is_some());
    assert_eq!(sketch.constraints.len(), 1);
}

#[test]
fn undo_set_constraint_value_restores_old_value() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    sketch.add_constraint(Constraint::Fixed(a));
    let ch = sketch.add_constraint(Constraint::Distance(a, b, 3.0));

    sketch.set_constraint_value(ch, 5.0).unwrap();
    // After undo, the distance target should be back to 3.0.
    assert!(sketch.undo());
    let entry = sketch.constraints.get(ch).unwrap();
    match &entry.constraint {
        Constraint::Distance(_, _, d) => assert_eq!(*d, 3.0),
        _ => panic!("expected Distance"),
    }
}

#[test]
fn multiple_undo_redo_cycles() {
    let mut sketch = Sketch::new();
    let _p1 = sketch.add_point(0.0, 0.0);
    let _p2 = sketch.add_point(1.0, 1.0);

    assert_eq!(sketch.points.len(), 2);
    assert_eq!(sketch.params.len(), 4);

    // Undo second add_point.
    assert!(sketch.undo());
    assert_eq!(sketch.points.len(), 1);
    assert_eq!(sketch.params.len(), 2);

    // Undo first add_point.
    assert!(sketch.undo());
    assert_eq!(sketch.points.len(), 0);
    assert_eq!(sketch.params.len(), 0);

    // Nothing more to undo.
    assert!(!sketch.undo());

    // Redo first add_point.
    assert!(sketch.redo());
    assert_eq!(sketch.points.len(), 1);

    // Redo second add_point.
    assert!(sketch.redo());
    assert_eq!(sketch.points.len(), 2);

    // Nothing more to redo.
    assert!(!sketch.redo());
}

#[test]
fn new_mutation_truncates_redo_history() {
    let mut sketch = Sketch::new();
    let _p1 = sketch.add_point(0.0, 0.0);
    let _p2 = sketch.add_point(1.0, 0.0);

    // Undo the second add.
    sketch.undo();
    assert_eq!(sketch.points.len(), 1);

    // Now add a different point -- should discard redo history.
    let _p3 = sketch.add_point(5.0, 5.0);
    assert_eq!(sketch.points.len(), 2);

    // Cannot redo (the old redo was discarded).
    assert!(!sketch.can_redo());
    assert!(!sketch.redo());

    // Can still undo back.
    assert!(sketch.undo());
    assert_eq!(sketch.points.len(), 1);
}

#[test]
fn can_undo_can_redo_initial_state() {
    let sketch = Sketch::new();
    assert!(!sketch.can_undo());
    assert!(!sketch.can_redo());
}

#[test]
fn undo_remove_element_restores_everything() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(1.0, 0.0);
    let l = sketch.add_line(p1, p2);
    let ch = sketch.add_constraint(Constraint::Horizontal(l));

    // Remove p1 -- cascades line and constraint.
    sketch.remove_element(ElementHandle::Point(p1));
    assert!(sketch.point_position(p1).is_none());
    assert!(sketch.lines.get(l).is_none());
    assert!(sketch.constraint_status(ch).is_none());

    // Undo the removal.
    assert!(sketch.undo());
    assert!(sketch.point_position(p1).is_some());
    assert!(sketch.lines.get(l).is_some());
    assert!(sketch.constraint_status(ch).is_some());
}
