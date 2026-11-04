//! Tests for multi-body boolean operations.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::multi_boolean::multi_body_boolean;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

#[test]
fn multi_body_two_disjoint_bodies_union() {
    let bodies = vec![cube_at(0.0, 0.0, 0.0), cube_at(5.0, 0.0, 0.0)];
    let body = multi_body_boolean(&bodies, |m| !m.is_empty(), &tol()).expect("union");
    assert_eq!(body.faces.len(), 12);
    assert!(body.validate_all().is_empty());
}

#[test]
fn multi_body_requires_two_bodies() {
    let bodies = vec![cube_at(0.0, 0.0, 0.0)];
    assert!(matches!(
        multi_body_boolean(&bodies, |_| true, &tol()),
        Err(crate::BoolError::DegenerateInput(_))
    ));
}

#[test]
fn multi_body_overlapping_cubes_union() {
    let bodies = vec![cube_at(0.0, 0.0, 0.0), cube_at(0.5, 0.3, 0.2)];
    let body = multi_body_boolean(&bodies, |m| !m.is_empty(), &tol()).expect("union");
    assert!(body.faces.len() > 6, "{} faces", body.faces.len());
    assert!(body.validate_all().is_empty());
}

/// "At least two of three" is one predicate, not a fold: the former
/// sequential pairwise evaluation inferred a union from it.
#[test]
fn multi_body_at_least_two_of_three() {
    let bodies = vec![
        cube_at(0.0, 0.0, 0.0),
        cube_at(0.5, 0.3, 0.2),
        cube_at(0.25, 0.6, 0.45),
    ];
    let at_least_two = multi_body_boolean(&bodies, |m| m.len() >= 2, &tol()).expect("select");
    let union = multi_body_boolean(&bodies, |m| !m.is_empty(), &tol()).expect("union");
    assert!(at_least_two.validate_all().is_empty());
    assert_ne!(at_least_two.faces.len(), union.faces.len());
}

/// Two identical cubes share every face: one region both hold, beside
/// a disjoint third.
#[test]
fn coincident_bodies_union_on_the_partition() {
    let bodies = vec![
        cube_at(0.0, 0.0, 0.0),
        cube_at(0.0, 0.0, 0.0),
        cube_at(5.0, 0.0, 0.0),
    ];
    let body = multi_body_boolean(&bodies, |m| !m.is_empty(), &tol()).expect("union");
    assert_eq!(body.faces.len(), 12);
    assert!(body.validate_all().is_empty());
}
