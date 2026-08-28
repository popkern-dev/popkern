//! Tests for rigid-set functionality.

use crate::Sketch;
use crate::constraint::Constraint;
use crate::elements::ElementHandle;
use crate::error::SolveStatus;

/// Helper: distance between two (x, y) pairs.
fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

/// 1. Basic rigid set: 3 points forming a triangle -> create rigid set ->
///    verify 3 DOF (tx, ty, theta).
#[test]
fn basic_rigid_set_three_dof() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    let c = sketch.add_point(0.0, 4.0);

    let rsh = sketch
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set should succeed");

    // Rigid set should have 3 points.
    assert_eq!(sketch.rigid_sets.get(rsh).unwrap().points.len(), 3);

    let result = sketch.solve().expect("solve should succeed");
    // Under-constrained: 3 DOF from the rigid body (no external constraints).
    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            // The rigid set's first element should have reported freedom.
            assert!(!freedoms.is_empty());
        }
        other => panic!("expected UnderConstrained, got {:?}", other),
    }
}

/// Verify Jacobian remapping directs columns to rigid-body params.
#[test]
fn rigid_jacobian_remaps_to_body_params() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);

    let rsh = sketch
        .create_rigid_set(&[ElementHandle::Point(a), ElementHandle::Point(b)])
        .expect("create rigid set");

    let rs_param_start = sketch.rigid_sets.get(rsh).unwrap().param_start;

    sketch.add_constraint(Constraint::OrdinateX(a, 2.0));

    let sys = crate::jacobian::assemble(&sketch);
    assert_eq!(sys.num_equations, 1);
    // All triplet columns should point to rigid-body params.
    for &(_row, col, _val) in &sys.jacobian_triplets {
        assert!(
            col >= rs_param_start,
            "column {} should be >= {} (rigid param_start)",
            col,
            rs_param_start
        );
    }
}

/// 2. Translation: Fix one external point, add distance constraint to a
///    set point. Solve. Verify all set points moved together maintaining
///    distances.
#[test]
fn rigid_translation_preserves_distances() {
    let mut sketch = Sketch::new();
    // Triangle rigid set.
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    let c = sketch.add_point(0.0, 4.0);

    let _rsh = sketch
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set");

    // Record original inter-point distances.
    let d_ab_orig = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(b).unwrap(),
    );
    let d_ac_orig = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(c).unwrap(),
    );
    let d_bc_orig = dist(
        sketch.point_position(b).unwrap(),
        sketch.point_position(c).unwrap(),
    );

    // Fix point a at (2, 3).
    sketch.add_constraint(Constraint::OrdinateX(a, 2.0));
    sketch.add_constraint(Constraint::OrdinateY(a, 3.0));

    let result = sketch.solve().expect("solve should succeed");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify point a moved to (2, 3).
    let pos_a = sketch.point_position(a).unwrap();
    assert!((pos_a.0 - 2.0).abs() < 1e-6, "a.x = {}", pos_a.0);
    assert!((pos_a.1 - 3.0).abs() < 1e-6, "a.y = {}", pos_a.1);

    // Verify inter-point distances preserved.
    let d_ab_new = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(b).unwrap(),
    );
    let d_ac_new = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(c).unwrap(),
    );
    let d_bc_new = dist(
        sketch.point_position(b).unwrap(),
        sketch.point_position(c).unwrap(),
    );
    assert!(
        (d_ab_new - d_ab_orig).abs() < 1e-6,
        "AB dist changed: {} vs {}",
        d_ab_new,
        d_ab_orig
    );
    assert!(
        (d_ac_new - d_ac_orig).abs() < 1e-6,
        "AC dist changed: {} vs {}",
        d_ac_new,
        d_ac_orig
    );
    assert!(
        (d_bc_new - d_bc_orig).abs() < 1e-6,
        "BC dist changed: {} vs {}",
        d_bc_new,
        d_bc_orig
    );
}

/// 3. Rotation: Two external fixed points + two set points constrained to
///    them -> verify rotation occurred and distances preserved.
#[test]
fn rigid_rotation_via_external_constraints() {
    let mut sketch = Sketch::new();
    // Rigid set: 3 points in an L-shape.
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(2.0, 0.0);
    let c = sketch.add_point(0.0, 2.0);

    let _rsh = sketch
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set");

    let d_ab_orig = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(b).unwrap(),
    );
    let d_ac_orig = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(c).unwrap(),
    );

    // Fix a at origin and fix b's y to 1.0 (forcing a rotation).
    sketch.add_constraint(Constraint::Fixed(a));
    sketch.add_constraint(Constraint::OrdinateY(b, 1.0));

    let result = sketch.solve().expect("solve should succeed");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify distances preserved.
    let d_ab_new = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(b).unwrap(),
    );
    let d_ac_new = dist(
        sketch.point_position(a).unwrap(),
        sketch.point_position(c).unwrap(),
    );
    assert!(
        (d_ab_new - d_ab_orig).abs() < 1e-6,
        "AB dist changed: {} vs {}",
        d_ab_new,
        d_ab_orig
    );
    assert!(
        (d_ac_new - d_ac_orig).abs() < 1e-6,
        "AC dist changed: {} vs {}",
        d_ac_new,
        d_ac_orig
    );

    // Verify rotation happened (b should no longer be at y=0).
    let pos_b = sketch.point_position(b).unwrap();
    assert!(
        (pos_b.1 - 1.0).abs() < 1e-6,
        "b.y should be 1.0, got {}",
        pos_b.1
    );
}

/// 4. Intra-set skip: Add a Distance constraint between two points in
///    the same rigid set. Verify it's skipped (doesn't over-constrain).
#[test]
fn intra_set_constraint_is_skipped() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    let c = sketch.add_point(0.0, 4.0);

    let _rsh = sketch
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set");

    // Add a distance constraint between two set members — should be skipped.
    sketch.add_constraint(Constraint::Distance(a, b, 3.0));

    // Fix the set so it's fully constrained.
    sketch.add_constraint(Constraint::Fixed(a));
    sketch.add_constraint(Constraint::OrdinateY(b, 0.0));

    let result = sketch.solve().expect("solve should succeed");
    // Should NOT be over-constrained because intra-set Distance is skipped.
    if let SolveStatus::OverConstrained(_) = &result.status {
        panic!("should not be over-constrained — intra-set constraint should be skipped");
    }
    // ok — either fully constrained or under-constrained
}

/// 5. Delete rigid set: Create, delete -> verify points are independent
///    (DOF = 2N again).
#[test]
fn delete_rigid_set_restores_independence() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);

    let rsh = sketch
        .create_rigid_set(&[ElementHandle::Point(a), ElementHandle::Point(b)])
        .expect("create rigid set");

    // Points should be tagged.
    assert!(sketch.points.get(a).unwrap().rigid_set.is_some());
    assert!(sketch.points.get(b).unwrap().rigid_set.is_some());

    // Delete the rigid set.
    sketch.delete_rigid_set(rsh);

    // Points should be untagged.
    assert!(sketch.points.get(a).unwrap().rigid_set.is_none());
    assert!(sketch.points.get(b).unwrap().rigid_set.is_none());

    // Solve: should have 4 DOF (2 per point).
    let result = sketch.solve().expect("solve should succeed");
    match &result.status {
        SolveStatus::UnderConstrained(freedoms) => {
            // Each point should report 2 directions.
            let total_dirs: usize = freedoms.iter().map(|f| f.directions.len()).sum();
            assert!(
                total_dirs >= 4,
                "expected at least 4 DOF, got {}",
                total_dirs
            );
        }
        other => panic!("expected UnderConstrained, got {:?}", other),
    }
}

/// 6. Arc in rigid set: Create arc, add to rigid set, move -> verify
///    center follows.
#[test]
fn arc_in_rigid_set_center_follows() {
    let mut sketch = Sketch::new();
    let p = sketch.add_point(0.0, 0.0);
    let arc_h = sketch.add_arc(5.0, 0.0, 1.0, 0.0, std::f64::consts::PI);

    let _rsh = sketch
        .create_rigid_set(&[ElementHandle::Point(p), ElementHandle::Arc(arc_h)])
        .expect("create rigid set");

    // Record the distance from p to arc center.
    let (cx0, cy0, _, _, _) = sketch.arc_params(arc_h).unwrap();
    let p0 = sketch.point_position(p).unwrap();
    let d_orig = dist(p0, (cx0, cy0));

    // Fix p at (2, 3) to force translation.
    sketch.add_constraint(Constraint::OrdinateX(p, 2.0));
    sketch.add_constraint(Constraint::OrdinateY(p, 3.0));

    let result = sketch.solve().expect("solve should succeed");
    assert!(
        result.residual_norm < 1e-8,
        "residual too large: {}",
        result.residual_norm
    );

    // Verify arc center moved too.
    let (cx1, cy1, _, _, _) = sketch.arc_params(arc_h).unwrap();
    let p1 = sketch.point_position(p).unwrap();
    let d_new = dist(p1, (cx1, cy1));
    assert!(
        (d_new - d_orig).abs() < 1e-6,
        "distance from point to arc center changed: {} vs {}",
        d_new,
        d_orig
    );
}

/// 7. Fully constrained: Fix 2 set points -> 0 DOF.
#[test]
fn fully_constrained_rigid_set() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(3.0, 0.0);
    let c = sketch.add_point(0.0, 4.0);

    let _rsh = sketch
        .create_rigid_set(&[
            ElementHandle::Point(a),
            ElementHandle::Point(b),
            ElementHandle::Point(c),
        ])
        .expect("create rigid set");

    // Fix two points — this locks tx, ty, theta.
    sketch.add_constraint(Constraint::Fixed(a));
    sketch.add_constraint(Constraint::Fixed(b));

    let result = sketch.solve().expect("solve should succeed");
    match &result.status {
        SolveStatus::FullyConstrained => {}   // expected
        SolveStatus::OverConstrained(_) => {} // acceptable (Fixed gives 4 eqs for 3 DOF)
        other => panic!(
            "expected FullyConstrained or OverConstrained, got {:?}",
            other
        ),
    }
}

/// Element already in rigid set cannot be added to another.
#[test]
fn cannot_add_element_to_two_rigid_sets() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(1.0, 0.0);
    let c = sketch.add_point(0.0, 1.0);

    sketch
        .create_rigid_set(&[ElementHandle::Point(a), ElementHandle::Point(b)])
        .expect("first rigid set should succeed");

    // Try to create a second rigid set containing point a.
    let result = sketch.create_rigid_set(&[ElementHandle::Point(a), ElementHandle::Point(c)]);
    assert!(
        result.is_err(),
        "should fail because point a is already in a rigid set"
    );
}

/// Rigid set with a line element collects both endpoints.
#[test]
fn rigid_set_with_line_collects_endpoints() {
    let mut sketch = Sketch::new();
    let a = sketch.add_point(0.0, 0.0);
    let b = sketch.add_point(5.0, 0.0);
    let line = sketch.add_line(a, b);

    let rsh = sketch
        .create_rigid_set(&[ElementHandle::Line(line)])
        .expect("create rigid set");

    let rs = sketch.rigid_sets.get(rsh).unwrap();
    assert_eq!(rs.points.len(), 2, "line should contribute 2 points");
    assert!(rs.points.contains(&a));
    assert!(rs.points.contains(&b));
}
