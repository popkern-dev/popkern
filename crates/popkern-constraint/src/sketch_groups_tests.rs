//! Tests for constraint group management.

use crate::Sketch;
use crate::constraint::Constraint;

/// A disabled group constraint should be excluded from the solver, leaving
/// constrained geometry unmoved.
#[test]
fn disabled_group_constraint_is_skipped_by_solver() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);

    // Fix point a so the solver has a reference.
    sketch.add_constraint(Constraint::Fixed(a));

    // Add a distance constraint in a group, then disable the group.
    let ch = sketch.add_constraint(Constraint::Distance(a, b, 10.0));
    let group = sketch.create_group("test_group");
    sketch.add_to_group(group, ch).unwrap();
    sketch.disable_group(group);

    // Solve — the distance constraint should be inactive.
    sketch.solve().unwrap();

    // b should remain at (3, 0), not pulled to (10, 0).
    let (bx, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx - 3.0).abs() < 0.1,
        "b.x should stay near 3.0 while group is disabled, got {bx}"
    );
}

/// Re-enabling a group makes the constraint active again.
#[test]
fn enabling_group_activates_constraint() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(a));

    let ch = sketch.add_constraint(Constraint::Distance(a, b, 10.0));
    let group = sketch.create_group("grp");
    sketch.add_to_group(group, ch).unwrap();
    sketch.disable_group(group);

    sketch.solve().unwrap();
    let (bx_disabled, _) = sketch.point_position(b).unwrap();
    assert!((bx_disabled - 3.0).abs() < 0.1);

    // Enable the group and solve again.
    sketch.enable_group(group);
    sketch.solve().unwrap();

    let (bx_enabled, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx_enabled - 10.0).abs() < 1e-6,
        "after enable, b.x should be 10.0, got {bx_enabled}"
    );
}

/// After deleting a group the constraint has no group and remains active.
#[test]
fn delete_group_clears_member_and_constraint_stays_active() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);

    sketch.add_constraint(Constraint::Fixed(a));

    let ch = sketch.add_constraint(Constraint::Distance(a, b, 10.0));
    let group = sketch.create_group("grp");
    sketch.add_to_group(group, ch).unwrap();
    sketch.disable_group(group);

    // Delete the group.
    sketch.delete_group(group);

    // Group handle is now invalid.
    assert_eq!(sketch.is_group_enabled(group), None);

    // Constraint has no group — it is active.
    assert_eq!(sketch.constraints.get(ch).unwrap().group, None);

    sketch.solve().unwrap();
    let (bx, _) = sketch.point_position(b).unwrap();
    assert!(
        (bx - 10.0).abs() < 1e-6,
        "constraint active after group deleted, b.x should be 10.0, got {bx}"
    );
}

/// add_to_group replaces any previous group membership.
#[test]
fn add_to_group_replaces_previous_membership() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);

    let ch = sketch.add_constraint(Constraint::Coincident(a, b));

    let g1 = sketch.create_group("g1");
    let g2 = sketch.create_group("g2");

    sketch.add_to_group(g1, ch).unwrap();
    assert_eq!(sketch.constraints.get(ch).unwrap().group, Some(g1));

    sketch.add_to_group(g2, ch).unwrap();
    assert_eq!(sketch.constraints.get(ch).unwrap().group, Some(g2));
}

/// is_group_enabled returns the correct value.
#[test]
fn is_group_enabled_returns_correct_value() {
    let mut sketch = Sketch::new();
    let g = sketch.create_group("g");

    assert_eq!(sketch.is_group_enabled(g), Some(true));
    sketch.disable_group(g);
    assert_eq!(sketch.is_group_enabled(g), Some(false));
    sketch.enable_group(g);
    assert_eq!(sketch.is_group_enabled(g), Some(true));
}

/// add_to_group with an invalid group handle returns an error.
#[test]
fn add_to_group_invalid_group_returns_error() {
    use crate::error::ConstraintError;
    use crate::groups::GroupHandle;
    use slotmap::KeyData;

    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let ch = sketch.add_constraint(Constraint::Coincident(a, b));

    // Manufacture a handle that is not in the sketch.
    let bad_group = GroupHandle::from(KeyData::from_ffi(u64::MAX));
    let result = sketch.add_to_group(bad_group, ch);
    assert!(matches!(result, Err(ConstraintError::InvalidHandle)));
}
