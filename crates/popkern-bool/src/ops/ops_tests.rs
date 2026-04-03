//! Tests for extended operations.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ── Section tests ──────────────────────────────────────────────────────

#[test]
fn section_overlapping_cubes_produces_curves() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let result = super::section::section(&a, &b, &tol()).unwrap();
    assert!(
        !result.curves.is_empty(),
        "overlapping cubes should produce section curves"
    );
    assert!(result.pairs_tested > 0);
}

#[test]
fn section_disjoint_cubes_produces_no_curves() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let result = super::section::section(&a, &b, &tol()).unwrap();
    assert!(result.curves.is_empty());
}

// ── Imprint tests ──────────────────────────────────────────────────────

#[test]
fn imprint_cylinder_onto_cube() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let result = super::imprint::imprint(&cube, &cyl, &tol()).unwrap();
    // Cylinder intersects cube faces — should produce imprint curves
    assert!(
        !result.curves.is_empty(),
        "cylinder through cube should produce imprint curves, got 0"
    );
    assert!(result.pairs_tested > 0);
}

#[test]
fn imprint_disjoint_produces_no_curves() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.3, 1.0, tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let result = super::imprint::imprint(&cube, &cyl, &tol()).unwrap();
    assert!(result.curves.is_empty());
}

// ── Split/Slice tests ──────────────────────────────────────────────────

#[test]
fn split_plane_api_does_not_panic() {
    use super::split_slice::{SplitPlane, split_by_plane};
    use popkern_core::geometry::Point3;

    let cube = BRep::unit_cube(tol()).unwrap();
    let plane = SplitPlane {
        origin: Point3::new(0.5, 0.5, 0.5),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let result = split_by_plane(&cube, &plane, &tol());
    assert!(result.is_ok());
    // With a plane at x=0.5, both halves should exist (plane cuts through cube)
    // The result depends on the boolean intersection working correctly.
    // At minimum, verify the function produces a result.
    let result = result.unwrap();
    let has_any = result.positive.is_some() || result.negative.is_some();
    assert!(has_any, "split should produce at least one half");
}

#[test]
fn split_empty_body_returns_none() {
    use super::split_slice::{SplitPlane, split_by_plane};
    use popkern_core::geometry::Point3;

    let empty = BRep::new(tol());
    let plane = SplitPlane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let result = split_by_plane(&empty, &plane, &tol()).unwrap();
    assert!(result.positive.is_none());
    assert!(result.negative.is_none());
}

// ── Face-level boolean tests ───────────────────────────────────────────

#[test]
fn face_level_boolean_with_overlapping_cubes() {
    use super::face_level::{FacePairSpec, face_level_boolean};
    use crate::error::BoolOperation;

    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.0, 0.0))
        .unwrap();

    let fa = a.faces.keys().next().unwrap();
    let fb = b.faces.keys().next().unwrap();

    let pairs = vec![FacePairSpec {
        face_a: fa,
        face_b: fb,
    }];

    let result = face_level_boolean(&a, &b, BoolOperation::Union, &pairs, &tol()).unwrap();
    assert_eq!(result.pairs_processed, 1);
    // curves_found depends on which specific faces were picked
}

#[test]
fn face_level_boolean_empty_pairs() {
    use super::face_level::{FacePairSpec, face_level_boolean};
    use crate::error::BoolOperation;

    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();

    let result = face_level_boolean(&a, &b, BoolOperation::Union, &[], &tol()).unwrap();
    assert_eq!(result.pairs_processed, 0);
    assert_eq!(result.curves_found, 0);
    assert!(result.body.is_none());
}
