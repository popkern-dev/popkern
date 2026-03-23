//! End-to-end integration tests for the boolean operations pipeline.
//!
//! These tests exercise the full public API (union, intersection, subtraction)
//! through the complete pipeline: intersect → split → classify → select → merge.

// Structural lint marker: this file is test code, included via #[cfg(test)] in lib.rs.
#[cfg(test)]
use popkern_core::{
    geometry::{Point3, Vec3},
    tolerance::Tolerance,
    topology::BRep,
};

use crate::csg::{CsgNode, csg_eval};
use crate::error::BoolError;

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
            match sof.geometry {
                SurfaceGeometry::Plane { ref mut origin, .. } => {
                    *origin = *origin + offset;
                }
                SurfaceGeometry::Cylinder { ref mut origin, .. } => {
                    *origin = *origin + offset;
                }
                SurfaceGeometry::Cone { ref mut apex, .. } => {
                    *apex = *apex + offset;
                }
                SurfaceGeometry::Sphere { ref mut center, .. } => {
                    *center = *center + offset;
                }
                SurfaceGeometry::Torus { ref mut center, .. } => {
                    *center = *center + offset;
                }
                _ => {} // NURBS and other variants: no origin to translate
            }
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve {
            match coe.geometry {
                CurveGeometry::Line { ref mut origin, .. } => {
                    *origin = *origin + offset;
                }
                CurveGeometry::Circle { ref mut center, .. }
                | CurveGeometry::Arc { ref mut center, .. } => {
                    *center = *center + offset;
                }
                _ => {} // NURBS and other variants: no origin to translate
            }
        }
    }
    brep
}

/// Scale a unit cube [0,1]^3 to [-1,2]^3 so it fully contains any other unit cube at origin.
fn scale_to_enclosing(mut brep: BRep) -> BRep {
    for (_, v) in &mut brep.vertices {
        v.position = Point3::new(
            v.position.x * 3.0 - 1.0,
            v.position.y * 3.0 - 1.0,
            v.position.z * 3.0 - 1.0,
        );
    }
    brep
}

// ── Overlapping cubes ──────────────────────────────────────────────────

#[test]
fn overlapping_cubes_union_succeeds() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let result = crate::union(a, b, tol());
    assert!(result.is_ok(), "Union should succeed: {:?}", result.err());
    let brep = result.unwrap();
    let errs = brep.validate_all();
    assert!(errs.is_empty(), "Union result should be valid: {:?}", errs);
    assert!(brep.faces.len() > 6);
}

#[test]
fn overlapping_cubes_intersection_succeeds() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let result = crate::intersection(a, b, tol());
    assert!(
        result.is_ok(),
        "Intersection should succeed: {:?}",
        result.err()
    );
    let brep = result.unwrap();
    let errs = brep.validate_all();
    assert!(
        errs.is_empty(),
        "Intersection result should be valid: {:?}",
        errs
    );
    assert!(!brep.faces.is_empty() && brep.faces.len() < 12);
}

#[test]
fn overlapping_cubes_subtraction_succeeds() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let result = crate::subtraction(a, b, tol());
    assert!(
        result.is_ok(),
        "Subtraction should succeed: {:?}",
        result.err()
    );
    let brep = result.unwrap();
    let errs = brep.validate_all();
    assert!(
        errs.is_empty(),
        "Subtraction result should be valid: {:?}",
        errs
    );
    assert!(!brep.faces.is_empty());
}

// ── Disjoint cubes ─────────────────────────────────────────────────────

#[test]
fn disjoint_cubes_union_has_12_faces() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::union(a, b, tol()).unwrap();
    assert_eq!(result.faces.len(), 12);
}

#[test]
fn disjoint_cubes_intersection_is_empty() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::intersection(a, b, tol()).unwrap();
    assert_eq!(result.faces.len(), 0);
}

#[test]
fn disjoint_cubes_subtraction_has_6_faces() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::subtraction(a, b, tol()).unwrap();
    assert_eq!(result.faces.len(), 6);
}

// ── Containment (A inside B) ──────────────────────────────────────────

#[test]
fn containment_union_returns_outer() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = scale_to_enclosing(BRep::unit_cube(tol()).unwrap());
    let result = crate::union(a, b, tol()).unwrap();
    assert_eq!(result.faces.len(), 6, "Union should return B (6 faces)");
}

#[test]
fn containment_intersection_returns_inner() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = scale_to_enclosing(BRep::unit_cube(tol()).unwrap());
    let result = crate::intersection(a, b, tol()).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "Intersection should return A (6 faces)"
    );
}

#[test]
fn containment_subtraction_returns_error() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = scale_to_enclosing(BRep::unit_cube(tol()).unwrap());
    let result = crate::subtraction(a, b, tol());
    assert!(
        result.is_err(),
        "Subtraction with A inside B should fail (requires KFMRH)"
    );
}

// ── _ref variants ──────────────────────────────────────────────────────

#[test]
fn ref_variants_preserve_inputs() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::union_ref(&a, &b, tol());
    assert!(result.is_ok());
    // Verify inputs are still usable after the borrow-based call.
    assert_eq!(a.faces.len(), 6);
    assert_eq!(b.faces.len(), 6);
}

#[test]
fn ref_intersection_matches_owned() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let ref_result = crate::intersection_ref(&a, &b, tol()).unwrap();
    let owned_result = crate::intersection(a.clone(), b.clone(), tol()).unwrap();
    assert_eq!(ref_result.faces.len(), owned_result.faces.len());
}

#[test]
fn ref_subtraction_matches_owned() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let ref_result = crate::subtraction_ref(&a, &b, tol()).unwrap();
    let owned_result = crate::subtraction(a.clone(), b.clone(), tol()).unwrap();
    assert_eq!(ref_result.faces.len(), owned_result.faces.len());
}

// ── CSG tree integration ───────────────────────────────────────────────

#[test]
fn csg_tree_with_three_disjoint_cubes() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let c = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(10.0, 0.0, 0.0));

    let tree = CsgNode::Union(
        Box::new(CsgNode::Union(
            Box::new(CsgNode::Leaf(Box::new(a))),
            Box::new(CsgNode::Leaf(Box::new(b))),
        )),
        Box::new(CsgNode::Leaf(Box::new(c))),
    );
    let result = csg_eval(tree, &tol()).unwrap();
    assert_eq!(result.faces.len(), 18); // 3 * 6
}

#[test]
fn csg_tree_subtraction_then_union() {
    // (A - B) ∪ C where all are disjoint: A - B = A, then A ∪ C = 12 faces.
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let c = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(10.0, 0.0, 0.0));

    let tree = CsgNode::Union(
        Box::new(CsgNode::Subtraction(
            Box::new(CsgNode::Leaf(Box::new(a))),
            Box::new(CsgNode::Leaf(Box::new(b))),
        )),
        Box::new(CsgNode::Leaf(Box::new(c))),
    );
    let result = csg_eval(tree, &tol()).unwrap();
    assert_eq!(result.faces.len(), 12); // A(6) + C(6)
}

// ── Error diagnostics ──────────────────────────────────────────────────

#[test]
fn coincident_cubes_returns_error() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap(); // Same position — faces coincide.
    let result = crate::union(a, b, tol());
    assert!(result.is_err(), "Coincident cubes should produce an error");
    // The error should be CoincidentFaces.
    let err = result.unwrap_err();
    assert!(
        matches!(err, BoolError::CoincidentFaces { .. }),
        "Expected CoincidentFaces error, got: {err:?}"
    );
}

#[test]
fn coincident_cubes_intersection_returns_error() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let result = crate::intersection(a, b, tol());
    assert!(result.is_err());
}

#[test]
fn coincident_cubes_subtraction_returns_error() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let result = crate::subtraction(a, b, tol());
    assert!(result.is_err());
}

// ── boolean() generic entry point ──────────────────────────────────────

#[test]
fn boolean_generic_union() {
    use crate::BoolOperation;
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::boolean(BoolOperation::Union, a, b, tol()).unwrap();
    assert_eq!(result.faces.len(), 12);
}

#[test]
fn boolean_ref_generic_intersection() {
    use crate::BoolOperation;
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::boolean_ref(BoolOperation::Intersection, &a, &b, tol()).unwrap();
    assert_eq!(result.faces.len(), 0);
    // Inputs still usable.
    assert_eq!(a.faces.len(), 6);
}

// ── Curved surface tests ──────────────────────────────────────────────

#[test]
fn cube_minus_sphere_diagnostic() {
    let t = Tolerance::default().with_linear(1e-3);
    // Unit cube [0,1]^3
    let cube = BRep::unit_cube(t.clone()).unwrap();
    // Sphere centered at (0.5, 0.5, 0.5) with radius 0.3 — fully inside cube
    let mut sphere = BRep::sphere(0.3, t.clone()).unwrap();
    for (_, v) in &mut sphere.vertices {
        v.position = v.position + Vec3::new(0.5, 0.5, 0.5);
    }
    // Also translate sphere surface geometry origins
    for (_, face) in &mut sphere.faces {
        if let Some(ref mut surf) = face.surface {
            match surf.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                    ref mut center, ..
                } => {
                    *center = *center + Vec3::new(0.5, 0.5, 0.5);
                }
                _ => {}
            }
        }
    }

    eprintln!(
        "Cube: {} faces, Sphere: {} faces",
        cube.faces.len(),
        sphere.faces.len()
    );
    let result = crate::subtraction(cube, sphere, t);
    match result {
        Ok(brep) => {
            eprintln!(
                "CUBE-SPHERE OK: {} faces, {} edges, {} verts",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len()
            );
            let errs = brep.validate_all();
            if errs.is_empty() {
                eprintln!("  VALIDATION: ALL PASS");
            } else {
                for e in &errs {
                    eprintln!("  VALIDATION: {e}");
                }
            }
        }
        Err(e) => eprintln!("CUBE-SPHERE ERR: {e}"),
    }
}

#[test]
fn cube_sphere_intersection_diagnostic() {
    let t = Tolerance::default().with_linear(1e-3);
    // Unit cube [0,1]^3
    let cube = BRep::unit_cube(t.clone()).unwrap();
    // Sphere centered at origin with radius 0.8 — intersects cube faces
    let sphere = BRep::sphere(0.8, t.clone()).unwrap();

    eprintln!(
        "Cube: {} faces, Sphere: {} faces",
        cube.faces.len(),
        sphere.faces.len()
    );
    let result = crate::subtraction(cube, sphere, t);
    match result {
        Ok(brep) => {
            eprintln!(
                "CUBE-SPHERE-ISECT OK: {} faces, {} edges, {} verts",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len()
            );
            let errs = brep.validate_all();
            if errs.is_empty() {
                eprintln!("  VALIDATION: ALL PASS");
            } else {
                for e in &errs {
                    eprintln!("  VALIDATION: {e}");
                }
            }
        }
        Err(e) => eprintln!("CUBE-SPHERE-ISECT ERR: {e}"),
    }
}

#[test]
fn cube_sphere_overlapping_diagnostic() {
    let t = Tolerance::default().with_linear(1e-3);
    // Unit cube [0,1]^3
    let cube = BRep::unit_cube(t.clone()).unwrap();
    // Large sphere centered at (0.5, 0.5, 0.5) radius 0.8 — pokes through all 6 faces
    let mut sphere = BRep::sphere(0.8, t.clone()).unwrap();
    for (_, v) in &mut sphere.vertices {
        v.position = v.position + Vec3::new(0.5, 0.5, 0.5);
    }
    for (_, face) in &mut sphere.faces {
        if let Some(ref mut surf) = face.surface {
            match surf.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                    ref mut center, ..
                } => {
                    *center = *center + Vec3::new(0.5, 0.5, 0.5);
                }
                _ => {}
            }
        }
    }

    // Check intersection stage directly
    let isect = crate::intersect::find_intersections(&cube, &sphere, &t);
    match &isect {
        Ok(r) => {
            eprintln!("INTERSECT: {} records", r.records.len());
            for (i, rec) in r.records.iter().enumerate() {
                eprintln!("  {i}: face_a={:?} face_b={:?}", rec.face_a, rec.face_b);
            }
        }
        Err(e) => eprintln!("INTERSECT ERR: {:?}", e),
    }

    let result = crate::subtraction(cube, sphere, t);
    match result {
        Ok(brep) => {
            eprintln!("RESULT OK: {} faces", brep.faces.len());
            let errs = brep.validate_all();
            for e in &errs {
                eprintln!("  VALIDATION: {e}");
            }
        }
        Err(e) => eprintln!("RESULT ERR: {e}"),
    }
}

#[test]
fn sphere_aabb_diagnostic() {
    let t = Tolerance::default().with_linear(1e-3);
    let mut sphere = BRep::sphere(0.8, t.clone()).unwrap();
    for (_, v) in &mut sphere.vertices {
        v.position = v.position + Vec3::new(0.5, 0.5, 0.5);
    }
    for (_, face) in &mut sphere.faces {
        if let Some(ref mut surf) = face.surface {
            match surf.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                    ref mut center, ..
                } => {
                    *center = *center + Vec3::new(0.5, 0.5, 0.5);
                }
                _ => {}
            }
        }
    }

    eprintln!(
        "Sphere: {} faces, {} verts, {} edges",
        sphere.faces.len(),
        sphere.vertices.len(),
        sphere.edges.len()
    );
    for (h, v) in &sphere.vertices {
        eprintln!(
            "  vertex {:?}: ({:.2}, {:.2}, {:.2})",
            h, v.position.x, v.position.y, v.position.z
        );
    }
    for (h, face) in &sphere.faces {
        if let Some(ref s) = face.surface {
            eprintln!("  face {:?}: {:?}", h, std::mem::discriminant(&s.geometry));
        }
    }
    let aabbs = crate::intersect::compute_face_aabbs(&sphere);
    for (h, aabb) in &aabbs {
        eprintln!(
            "  aabb {:?}: ({:.2},{:.2},{:.2}) - ({:.2},{:.2},{:.2})",
            h, aabb.min.x, aabb.min.y, aabb.min.z, aabb.max.x, aabb.max.y, aabb.max.z
        );
    }
}
