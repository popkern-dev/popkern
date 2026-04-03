//! Tests for multi-body boolean operations.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

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
    let result =
        super::multi_boolean::multi_body_boolean(&bodies, |m| m & 1 != 0 || m & 2 != 0, &tol());
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.body_count, 2);
    assert!(result.body.is_some());
}

#[test]
fn multi_body_requires_two_bodies() {
    let bodies = vec![BRep::unit_cube(tol()).unwrap()];
    let result = super::multi_boolean::multi_body_boolean(&bodies, |_| true, &tol());
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
    // Union of 3 disjoint cubes
    let result = super::multi_boolean::multi_body_boolean(
        &bodies,
        |m| m != 0, // any body
        &tol(),
    );
    assert!(result.is_ok());
    let result = result.unwrap();
    assert_eq!(result.body_count, 3);
    assert!(result.body.is_some());
}
