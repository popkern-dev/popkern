//! Tests for constrained repair (heal_with_constraints).

use super::*;
use crate::geometry::Point3;
use crate::healing::constraints::{ConstraintSet, GeometricConstraint, check_constraints};
use crate::healing::options::HealOptions;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn heal_with_constraints_empty_set_same_as_regular() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let options = HealOptions::new(Tolerance::default());
    let constraints = ConstraintSet::new();

    let result = heal_with_constraints(&mut brep, &options, &constraints);
    assert!(result.is_ok());
    let r = result.unwrap();
    assert_eq!(r.constraints_checked, 0);
    assert!(!r.was_modified);
}

#[test]
fn heal_with_constraints_fixed_point_preserved() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, v) = brep.vertices.iter().next().unwrap();
    let original_pos = v.position;

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: original_pos,
        max_displacement: 1e-6,
    });

    let options = HealOptions::new(Tolerance::default());
    let result = heal_with_constraints(&mut brep, &options, &cs);
    assert!(result.is_ok());

    // After healing with constraint, vertex should still be at original position.
    let actual = brep.vertices[vh].position;
    assert!(
        actual.distance_to(&original_pos) < 1e-5,
        "vertex should stay at fixed position"
    );
}

#[test]
fn heal_with_constraints_on_empty_brep() {
    let mut brep = BRep::new(Tolerance::default());
    let cs = ConstraintSet::new();
    let options = HealOptions::new(Tolerance::default());

    let result = heal_with_constraints(&mut brep, &options, &cs);
    assert!(result.is_ok());
    let r = result.unwrap();
    assert_eq!(r.constraints_checked, 0);
}

#[test]
fn heal_with_constraints_satisfied_before_remain_satisfied() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    // Coplanarity of a single face — trivially satisfied.
    let fh = brep.faces.iter().next().unwrap().0;
    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::Coplanarity {
        faces: vec![fh],
        max_deviation: 1e-3,
    });

    let options = HealOptions::new(Tolerance::default());
    let result = heal_with_constraints(&mut brep, &options, &cs).unwrap();

    assert_eq!(result.constraints_checked, 1);
    assert_eq!(result.satisfied_before, 1);
    assert_eq!(result.satisfied_after, 1);
    assert!(result.violations_prevented.is_empty());
}

#[test]
fn constrained_repair_result_fields() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, v) = brep.vertices.iter().next().unwrap();

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: v.position,
        max_displacement: 1e-6,
    });

    let options = HealOptions::new(Tolerance::default());
    let result = heal_with_constraints(&mut brep, &options, &cs).unwrap();

    // On a valid cube, healing is essentially a no-op.
    assert_eq!(result.heal_result.defects_before, 0);
    assert_eq!(result.heal_result.defects_after, 0);
    assert_eq!(result.constraints_checked, 1);
}

#[test]
fn heal_with_multiple_constraints() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let mut verts = brep.vertices.iter();
    let (v1, vd1) = verts.next().unwrap();
    let (v2, vd2) = verts.next().unwrap();

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: v1,
        position: vd1.position,
        max_displacement: 1e-6,
    });
    cs.add(GeometricConstraint::FixedPoint {
        vertex: v2,
        position: vd2.position,
        max_displacement: 1e-6,
    });

    let options = HealOptions::new(Tolerance::default());
    let result = heal_with_constraints(&mut brep, &options, &cs).unwrap();
    assert_eq!(result.constraints_checked, 2);
    assert_eq!(result.satisfied_before, 2);
    assert_eq!(result.satisfied_after, 2);
}

#[test]
fn heal_with_pre_violated_constraint_not_blamed() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, _) = brep.vertices.iter().next().unwrap();

    // Create a constraint that is already violated (wrong position).
    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: Point3::new(999.0, 999.0, 999.0),
        max_displacement: 0.001,
    });

    // Verify it's violated before repair.
    let pre_violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert_eq!(pre_violations.len(), 1);

    let options = HealOptions::new(Tolerance::default());
    let result = heal_with_constraints(&mut brep, &options, &cs).unwrap();

    // The pre-existing violation should not be in violations_prevented.
    assert!(
        result.violations_prevented.is_empty(),
        "pre-existing violations should not be blamed on repair"
    );
}
