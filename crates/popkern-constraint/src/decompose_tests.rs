//! Tests for graph decomposition (union-find and connected components).

use super::*;
use crate::constraint::Constraint;
use crate::Sketch;

#[test]
fn union_find_basic() {
    let mut uf = UnionFind::new(5);
    assert_ne!(uf.find(0), uf.find(1));

    uf.union(0, 1);
    assert_eq!(uf.find(0), uf.find(1));

    uf.union(2, 3);
    assert_eq!(uf.find(2), uf.find(3));
    assert_ne!(uf.find(0), uf.find(2));

    uf.union(1, 3);
    assert_eq!(uf.find(0), uf.find(3));
}

#[test]
fn two_unconnected_fixed_points_two_components() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(1.0, 2.0);
    let p2 = sketch.add_point(3.0, 4.0);

    sketch.add_constraint(Constraint::Fixed(p1));
    sketch.add_constraint(Constraint::Fixed(p2));

    let components = find_components(&sketch);
    assert_eq!(
        components.len(),
        2,
        "two independent fixed points should yield 2 components"
    );

    // Each component should have exactly 2 params and 1 constraint.
    for comp in &components {
        assert_eq!(comp.global_indices.len(), 2);
        assert_eq!(comp.constraints.len(), 1);
    }
}

#[test]
fn two_points_with_distance_one_component() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 4.0);

    sketch.add_constraint(Constraint::Distance(p1, p2, 5.0));

    let components = find_components(&sketch);
    assert_eq!(
        components.len(),
        1,
        "distance between two points should merge them into 1 component"
    );

    let comp = &components[0];
    assert_eq!(comp.global_indices.len(), 4); // 2 params per point
    assert_eq!(comp.constraints.len(), 1);
}

#[test]
fn driven_constraints_excluded() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(1.0, 1.0);

    // Driven constraint should not contribute to components.
    sketch.add_driven(Constraint::Distance(p1, p2, 1.0));

    let components = find_components(&sketch);
    assert!(
        components.is_empty(),
        "driven constraint should not create any components"
    );
}

#[test]
fn empty_sketch_no_components() {
    let sketch = Sketch::new();
    let components = find_components(&sketch);
    assert!(components.is_empty());
}

#[test]
fn line_constraint_connects_endpoints() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let p2 = sketch.add_point(3.0, 2.0);
    let l = sketch.add_line(p1, p2);

    sketch.add_constraint(Constraint::Horizontal(l));

    let components = find_components(&sketch);
    assert_eq!(components.len(), 1);
    assert_eq!(components[0].global_indices.len(), 4);
}

#[test]
fn two_independent_rectangles_two_components() {
    let mut sketch = Sketch::new();

    // Rectangle 1.
    let a1 = sketch.add_point(0.0, 0.0);
    let b1 = sketch.add_point(3.0, 0.0);
    let c1 = sketch.add_point(3.0, 2.0);
    let d1 = sketch.add_point(0.0, 2.0);
    let l1_bottom = sketch.add_line(a1, b1);
    let l1_right = sketch.add_line(b1, c1);
    let _l1_top = sketch.add_line(c1, d1);
    let _l1_left = sketch.add_line(d1, a1);
    sketch.add_constraint(Constraint::Fixed(a1));
    sketch.add_constraint(Constraint::Horizontal(l1_bottom));
    sketch.add_constraint(Constraint::Vertical(l1_right));

    // Rectangle 2 (independent).
    let a2 = sketch.add_point(10.0, 10.0);
    let b2 = sketch.add_point(14.0, 10.0);
    let c2 = sketch.add_point(14.0, 13.0);
    let d2 = sketch.add_point(10.0, 13.0);
    let l2_bottom = sketch.add_line(a2, b2);
    let l2_right = sketch.add_line(b2, c2);
    let _l2_top = sketch.add_line(c2, d2);
    let _l2_left = sketch.add_line(d2, a2);
    sketch.add_constraint(Constraint::Fixed(a2));
    sketch.add_constraint(Constraint::Horizontal(l2_bottom));
    sketch.add_constraint(Constraint::Vertical(l2_right));

    let components = find_components(&sketch);
    assert_eq!(
        components.len(),
        2,
        "two independent rectangles should yield 2 components"
    );
}

#[test]
fn local_indices_are_correct() {
    let mut sketch = Sketch::new();
    let p1 = sketch.add_point(0.0, 0.0);
    let _p2 = sketch.add_point(5.0, 5.0); // unconnected
    let p3 = sketch.add_point(1.0, 1.0);

    sketch.add_constraint(Constraint::Distance(p1, p3, 2.0));
    sketch.add_constraint(Constraint::Fixed(_p2));

    let components = find_components(&sketch);
    assert_eq!(components.len(), 2);

    for comp in &components {
        // Verify local_indices are consistent with global_indices.
        for (local_idx, &global_idx) in comp.global_indices.iter().enumerate() {
            assert_eq!(
                comp.local_indices.get(&global_idx),
                Some(&local_idx),
                "local_indices mismatch for global index {}",
                global_idx
            );
        }
    }
}
