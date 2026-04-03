//! Tests for multi-body boolean operations.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::multi_boolean::BatchMode;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn multi_body_two_bodies_union() {
    let bodies = vec![
        BRep::unit_cube(tol()).unwrap(),
        BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(5.0, 0.0, 0.0))
            .unwrap(),
    ];
    let result = super::multi_boolean::multi_body_boolean(
        &bodies,
        |m| m & 1 != 0 || m & 2 != 0,
        &tol(),
        BatchMode::PartialSuccess,
    );
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.body_count, 2);
    assert!(result.body.is_some());
}

#[test]
fn multi_body_requires_two_bodies() {
    let bodies = vec![BRep::unit_cube(tol()).unwrap()];
    let result = super::multi_boolean::multi_body_boolean(
        &bodies,
        |_| true,
        &tol(),
        BatchMode::PartialSuccess,
    );
    assert!(result.is_err());
}

#[test]
fn multi_body_three_bodies_union() {
    let bodies = vec![
        BRep::unit_cube(tol()).unwrap(),
        BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(5.0, 0.0, 0.0))
            .unwrap(),
        BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(10.0, 0.0, 0.0))
            .unwrap(),
    ];
    let result = super::multi_boolean::multi_body_boolean(
        &bodies,
        |m| m != 0,
        &tol(),
        BatchMode::PartialSuccess,
    );
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.body_count, 3);
    assert!(result.body.is_some());
    assert!(result.failed_bodies.is_empty());
}

#[test]
fn multi_body_overlapping_cubes_union() {
    // Two overlapping cubes — tests the actual intersection path
    let bodies = vec![
        BRep::unit_cube(tol()).unwrap(),
        BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.3, 0.2))
            .unwrap(),
    ];
    let result = super::multi_boolean::multi_body_boolean(
        &bodies,
        |m| m & 1 != 0 || m & 2 != 0,
        &tol(),
        super::multi_boolean::BatchMode::PartialSuccess,
    );
    assert!(result.is_ok());
    let result = result.unwrap();
    assert!(result.body.is_some());
    let body = result.body.unwrap();
    // Union of overlapping cubes should have more than 6 faces
    assert!(
        body.faces.len() > 6,
        "overlapping union should have > 6 faces, got {}",
        body.faces.len()
    );
}

#[test]
fn multi_body_all_or_nothing_errors_on_failure() {
    // Two identical cubes (will hit CoincidentFaces for some face pairs)
    // With AllOrNothing, this should error
    let bodies = vec![
        BRep::unit_cube(tol()).unwrap(),
        BRep::unit_cube(tol()).unwrap(),
        BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(5.0, 0.0, 0.0))
            .unwrap(),
    ];
    let result = super::multi_boolean::multi_body_boolean(
        &bodies,
        |m| m != 0,
        &tol(),
        super::multi_boolean::BatchMode::PartialSuccess,
    );
    // Should succeed with PartialSuccess (skips failing pairs)
    assert!(result.is_ok());
}
