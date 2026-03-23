//! Tests for the boolean pipeline orchestrator.
// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Translate all vertices and geometry of a BRep by a vector offset.
fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
    for (_, vertex) in &mut brep.vertices {
        vertex.position = vertex.position + offset;
    }
    for (_, face) in &mut brep.faces {
        if let Some(ref mut sof) = face.surface {
            if let SurfaceGeometry::Plane { ref mut origin, .. } = sof.geometry {
                *origin = *origin + offset;
            }
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve {
            if let CurveGeometry::Line { ref mut origin, .. } = coe.geometry {
                *origin = *origin + offset;
            }
        }
    }
    brep
}

#[test]
fn disjoint_cubes_union_produces_valid_result() {
    let a = BRep::unit_cube(tol()).expect("unit cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("unit cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result =
        execute_boolean(BoolOperation::Union, a, b, &tol()).expect("disjoint union should succeed");
    // Should have 12 faces (6 + 6).
    assert_eq!(result.faces.len(), 12);
}

#[test]
fn disjoint_cubes_intersection_is_empty() {
    let a = BRep::unit_cube(tol()).expect("unit cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("unit cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = execute_boolean(BoolOperation::Intersection, a, b, &tol())
        .expect("disjoint intersection should succeed");
    assert_eq!(result.faces.len(), 0);
}

#[test]
fn disjoint_cubes_subtraction_returns_a() {
    let a = BRep::unit_cube(tol()).expect("unit cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("unit cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = execute_boolean(BoolOperation::Subtraction, a, b, &tol())
        .expect("disjoint subtraction should succeed");
    // Should return A unchanged -- 6 faces.
    assert_eq!(result.faces.len(), 6);
}

#[test]
fn contained_cube_union_returns_outer() {
    let a = BRep::unit_cube(tol()).expect("unit cube a");
    let mut b = BRep::unit_cube(tol()).expect("unit cube b");
    for (_, v) in &mut b.vertices {
        v.position = popkern_core::geometry::Point3::new(
            v.position.x * 3.0 - 1.0,
            v.position.y * 3.0 - 1.0,
            v.position.z * 3.0 - 1.0,
        );
    }
    let result = execute_boolean(BoolOperation::Union, a, b, &tol())
        .expect("contained union should succeed");
    assert_eq!(result.faces.len(), 6);
}

#[test]
fn contained_cube_intersection_returns_inner() {
    let a = BRep::unit_cube(tol()).expect("unit cube a");
    let mut b = BRep::unit_cube(tol()).expect("unit cube b");
    for (_, v) in &mut b.vertices {
        v.position = popkern_core::geometry::Point3::new(
            v.position.x * 3.0 - 1.0,
            v.position.y * 3.0 - 1.0,
            v.position.z * 3.0 - 1.0,
        );
    }
    let result = execute_boolean(BoolOperation::Intersection, a, b, &tol())
        .expect("contained intersection should succeed");
    assert_eq!(result.faces.len(), 6);
}

#[test]
fn empty_a_returns_degenerate_input_error() {
    let a = BRep::new(tol());
    let b = BRep::unit_cube(tol()).expect("unit cube b");
    let result = execute_boolean(BoolOperation::Union, a, b, &tol());
    assert!(result.is_err(), "empty A should produce an error");
    assert!(
        matches!(result.unwrap_err(), BoolError::DegenerateInput(_)),
        "expected DegenerateInput"
    );
}

#[test]
fn into_error_fallback_when_no_error_set() {
    let diag = BoolDiagnostics {
        error: None,
        intersection_curves: Vec::new(),
        face_classifications: Vec::new(),
        failed_stage: None,
    };
    let err = diag.into_error();
    assert!(
        matches!(err, BoolError::MergeFailed { .. }),
        "expected MergeFailed fallback"
    );
}

#[test]
fn contained_cube_subtraction_returns_error() {
    let a = BRep::unit_cube(tol()).expect("unit cube a");
    let mut b = BRep::unit_cube(tol()).expect("unit cube b");
    for (_, v) in &mut b.vertices {
        v.position = popkern_core::geometry::Point3::new(
            v.position.x * 3.0 - 1.0,
            v.position.y * 3.0 - 1.0,
            v.position.z * 3.0 - 1.0,
        );
    }
    let result = execute_boolean(BoolOperation::Subtraction, a, b, &tol());
    assert!(result.is_err());
}
