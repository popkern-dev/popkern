//! Tests for reference geometry import from BRep face edges.

use popkern_constraint::ElementHandle;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::plane::SketchPlane;
use crate::SketchOnPlane;

#[test]
fn import_cube_face_produces_4_lines_and_4_points() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let plane = SketchPlane::from_face(&cube, face).unwrap();
    let mut sop = SketchOnPlane::new(plane);

    let imported = sop.import_face_edges(&cube, face).unwrap();

    // A cube face has 4 edges and 4 vertices.
    let line_count = imported
        .iter()
        .filter(|e| matches!(e, ElementHandle::Line(_)))
        .count();
    let point_count = imported
        .iter()
        .filter(|e| matches!(e, ElementHandle::Point(_)))
        .count();

    assert_eq!(line_count, 4, "expected 4 construction lines");
    assert_eq!(point_count, 4, "expected 4 construction points");
}

#[test]
fn imported_elements_are_construction() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let plane = SketchPlane::from_face(&cube, face).unwrap();
    let mut sop = SketchOnPlane::new(plane);

    let imported = sop.import_face_edges(&cube, face).unwrap();

    for element in &imported {
        assert_eq!(
            sop.sketch.is_construction(element),
            Some(true),
            "imported element {element:?} should be construction"
        );
    }
}

#[test]
fn imported_points_at_correct_2d_positions() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let plane = SketchPlane::from_face(&cube, face).unwrap();
    let mut sop = SketchOnPlane::new(plane);

    let imported = sop.import_face_edges(&cube, face).unwrap();

    // All imported points should lie on the unit cube face, meaning their
    // 2D coordinates should be within the [0, 1] range (for a unit cube
    // face projected onto its own plane).
    for element in &imported {
        if let ElementHandle::Point(ph) = element {
            if let Some((x, y)) = sop.sketch.point_position(*ph) {
                // Unit cube vertex coords in the face's local frame should be
                // within [0, 1] with some tolerance.
                assert!(
                    x > -0.1 && x < 1.1,
                    "x={x} should be in [0,1] range for unit cube face"
                );
                assert!(
                    y > -0.1 && y < 1.1,
                    "y={y} should be in [0,1] range for unit cube face"
                );
            }
        }
    }
}

#[test]
fn no_duplicate_points_for_shared_vertices() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let face = cube.faces.keys().next().unwrap();

    let plane = SketchPlane::from_face(&cube, face).unwrap();
    let mut sop = SketchOnPlane::new(plane);

    let imported = sop.import_face_edges(&cube, face).unwrap();

    // 4 unique points (not 8 = 4 edges x 2 endpoints).
    let point_count = imported
        .iter()
        .filter(|e| matches!(e, ElementHandle::Point(_)))
        .count();
    assert_eq!(point_count, 4, "shared vertices should not be duplicated");
}
