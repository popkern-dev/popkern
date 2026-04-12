use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_has_no_self_intersections() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let pairs = detect_self_intersecting_faces(&brep);
    assert!(
        pairs.is_empty(),
        "valid cube should have no self-intersecting faces"
    );
}

#[test]
fn repair_on_valid_cube_face_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let result = repair_self_intersection(&mut brep, face);
    // Cube faces are planar (analytical), so repair returns immediately.
    assert!(
        !result.repaired,
        "analytical surface should not be repaired"
    );
    assert!(
        result.description.contains("Analytical"),
        "expected analytical surface message, got: {}",
        result.description
    );
}

#[test]
fn repair_on_missing_face_returns_not_found() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = FaceHandle::default();
    let result = repair_self_intersection(&mut brep, bogus);
    assert!(!result.repaired);
    assert!(result.description.contains("not found"));
}

#[test]
fn detect_on_empty_brep_returns_empty() {
    let brep = BRep::new(Tolerance::default());
    let pairs = detect_self_intersecting_faces(&brep);
    assert!(pairs.is_empty());
}

#[test]
fn repair_by_splitting_on_valid_cube_face_no_self_intersection() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = brep.faces.keys().next().unwrap();
    let result = repair_by_splitting(&mut brep, face, &tol);
    // Cube faces are planar (analytical), not NURBS.
    assert!(!result.repaired, "analytical surface should not be split");
    assert!(
        result.description.contains("not NURBS")
            || result.description.contains("Analytical")
            || result.description.contains("Only NURBS"),
        "expected non-NURBS message, got: {}",
        result.description
    );
}
