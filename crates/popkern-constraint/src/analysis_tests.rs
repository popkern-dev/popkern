//! Tests for DOF analysis and constraint status detection.

use crate::Sketch;
use crate::constraint::{Constraint, ConstraintStatus};
use crate::error::SolveStatus;

#[test]
fn unconstrained_point_has_two_dof() {
    let mut sketch = Sketch::new();
    let _p = sketch.add_point(1.0, 2.0);

    let result = sketch.solve().unwrap();
    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            // One point with 2 DOF.
            assert!(!freedoms.is_empty(), "should report at least one freedom");
            let total_dirs: usize = freedoms.iter().map(|f| f.directions.len()).sum();
            assert_eq!(
                total_dirs, 2,
                "should have 2 free directions, got {}",
                total_dirs
            );
        }
        other => panic!("expected UnderConstrained, got {:?}", other),
    }
}

#[test]
fn fixed_point_zero_dof() {
    let mut sketch = Sketch::new();
    let p = sketch.add_point(1.0, 2.0);
    sketch.add_constraint(Constraint::Fixed(p));

    let result = sketch.solve().unwrap();
    assert!(
        matches!(result.status, SolveStatus::FullyConstrained),
        "expected FullyConstrained, got {:?}",
        result.status
    );
}

#[test]
fn point_on_horizontal_line_one_dof() {
    // A point with one constraint (Horizontal on its line) should have 1 DOF
    // remaining on the constrained point.
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.5);
    let l = sketch.add_line(p1, p2);

    // Fix p1, make line horizontal: p2 still free in x.
    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Horizontal(l));

    let result = sketch.solve().unwrap();
    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            assert!(!freedoms.is_empty(), "should have at least one free DOF");
            // The free DOF should be the x-direction of p2.
            let total_dirs: usize = freedoms.iter().map(|f| f.directions.len()).sum();
            assert_eq!(
                total_dirs, 1,
                "should have exactly 1 free direction, got {}",
                total_dirs
            );
        }
        other => panic!("expected UnderConstrained, got {:?}", other),
    }
}

#[test]
fn duplicate_fixed_is_over_constrained() {
    let mut sketch = Sketch::new();
    let p = sketch.add_point(1.0, 2.0);

    let _c1 = sketch.add_constraint(Constraint::Fixed(p));
    let _c2 = sketch.add_constraint(Constraint::Fixed(p));

    let result = sketch.solve().unwrap();
    match &result.status {
        SolveStatus::OverConstrained(handles) => {
            assert!(
                !handles.is_empty(),
                "should identify at least one conflicting constraint"
            );
        }
        other => panic!("expected OverConstrained, got {:?}", other),
    }
}

#[test]
fn constraint_status_satisfied_after_solve() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 0.5);
    let l = sketch.add_line(p1, p2);

    let c_fix = sketch.add_constraint(Constraint::Fixed(p1));
    let c_horiz = sketch.add_constraint(Constraint::Horizontal(l));
    let c_dist = sketch.add_constraint(Constraint::Distance(p1, p2, 3.0));

    let _result = sketch.solve().unwrap();

    assert_eq!(
        sketch.constraint_status(c_fix),
        Some(ConstraintStatus::Satisfied)
    );
    assert_eq!(
        sketch.constraint_status(c_horiz),
        Some(ConstraintStatus::Satisfied)
    );
    assert_eq!(
        sketch.constraint_status(c_dist),
        Some(ConstraintStatus::Satisfied)
    );
}

#[test]
fn two_unconstrained_points_four_dof() {
    let mut sketch = Sketch::new();
    let _p1 = sketch.add_point(0.0, 0.0);
    let _p2 = sketch.add_point(1.0, 1.0);

    let result = sketch.solve().unwrap();
    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            let total_dirs: usize = freedoms.iter().map(|f| f.directions.len()).sum();
            assert_eq!(
                total_dirs, 4,
                "two points should have 4 DOF, got {}",
                total_dirs
            );
        }
        other => panic!("expected UnderConstrained, got {:?}", other),
    }
}
