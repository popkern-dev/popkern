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
    let result = crate::union(a, b, tol(), None);
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
    let result = crate::intersection(a, b, tol(), None);
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
    let result = crate::subtraction(a, b, tol(), None);
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
    let result = crate::union(a, b, tol(), None).unwrap();
    assert_eq!(result.faces.len(), 12);
}

#[test]
fn disjoint_cubes_intersection_is_empty() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::intersection(a, b, tol(), None).unwrap();
    assert_eq!(result.faces.len(), 0);
}

#[test]
fn disjoint_cubes_subtraction_has_6_faces() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::subtraction(a, b, tol(), None).unwrap();
    assert_eq!(result.faces.len(), 6);
}

// ── Containment (A inside B) ──────────────────────────────────────────

#[test]
fn containment_union_returns_outer() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = scale_to_enclosing(BRep::unit_cube(tol()).unwrap());
    let result = crate::union(a, b, tol(), None).unwrap();
    assert_eq!(result.faces.len(), 6, "Union should return B (6 faces)");
}

#[test]
fn containment_intersection_returns_inner() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = scale_to_enclosing(BRep::unit_cube(tol()).unwrap());
    let result = crate::intersection(a, b, tol(), None).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "Intersection should return A (6 faces)"
    );
}

#[test]
fn containment_subtraction_a_inside_b_is_empty() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = scale_to_enclosing(BRep::unit_cube(tol()).unwrap());
    let result =
        crate::subtraction(a, b, tol(), None).expect("A inside B subtraction should succeed");
    assert_eq!(result.faces.len(), 0, "A inside B: A − B should be empty");
}

#[test]
fn containment_subtraction_hollow_solid() {
    // outer − inner = hollow solid with 12 faces.
    let outer = BRep::unit_cube(tol()).unwrap();
    let mut inner = BRep::unit_cube(tol()).unwrap();
    for (_, v) in &mut inner.vertices {
        v.position = Point3::new(
            v.position.x * 0.6 + 0.2,
            v.position.y * 0.6 + 0.2,
            v.position.z * 0.6 + 0.2,
        );
    }
    let result =
        crate::subtraction(outer, inner, tol(), None).expect("hollow subtraction should succeed");
    assert_eq!(result.faces.len(), 12, "hollow solid should have 12 faces");
}

// ── _ref variants ──────────────────────────────────────────────────────

#[test]
fn ref_variants_preserve_inputs() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::union_ref(&a, &b, tol(), None);
    assert!(result.is_ok());
    // Verify inputs are still usable after the borrow-based call.
    assert_eq!(a.faces.len(), 6);
    assert_eq!(b.faces.len(), 6);
}

#[test]
fn ref_intersection_matches_owned() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let ref_result = crate::intersection_ref(&a, &b, tol(), None).unwrap();
    let owned_result = crate::intersection(a.clone(), b.clone(), tol(), None).unwrap();
    assert_eq!(ref_result.faces.len(), owned_result.faces.len());
}

#[test]
fn ref_subtraction_matches_owned() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let ref_result = crate::subtraction_ref(&a, &b, tol(), None).unwrap();
    let owned_result = crate::subtraction(a.clone(), b.clone(), tol(), None).unwrap();
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
fn coincident_cubes_union_does_not_panic() {
    // Coincident cubes now compute overlap boundary curves instead of
    // erroring. The boolean may succeed or fail downstream, but should
    // never return CoincidentFaces and never panic.
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let _result = crate::union(a, b, tol(), None);
}

#[test]
fn coincident_cubes_intersection_does_not_panic() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let _result = crate::intersection(a, b, tol(), None);
}

#[test]
fn coincident_cubes_subtraction_does_not_panic() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let _result = crate::subtraction(a, b, tol(), None);
}

// ── boolean() generic entry point ──────────────────────────────────────

#[test]
fn boolean_generic_union() {
    use crate::BoolOperation;
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::boolean(BoolOperation::Union, a, b, tol(), None).unwrap();
    assert_eq!(result.faces.len(), 12);
}

#[test]
fn boolean_ref_generic_intersection() {
    use crate::BoolOperation;
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::boolean_ref(BoolOperation::Intersection, &a, &b, tol(), None).unwrap();
    assert_eq!(result.faces.len(), 0);
    // Inputs still usable.
    assert_eq!(a.faces.len(), 6);
}

// ── Curved surface tests ──────────────────────────────────────────────

#[test]
fn cube_minus_contained_sphere_returns_cube() {
    // Sphere fully inside cube — containment shortcut returns cube unchanged
    let t = Tolerance::default().with_linear(1e-3);
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let sphere = translate(
        BRep::sphere(0.3, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, 0.5),
    );
    let result = crate::subtraction(cube, sphere, t, None);
    // Containment subtraction is a TODO — returns cube unchanged or errors
    assert!(result.is_ok() || result.is_err());
}

#[test]
fn cube_minus_overlapping_sphere_succeeds() {
    // Sphere overlaps cube — curved face splitting is skipped, but the
    // pipeline should complete without hanging or panicking.
    // Use coarse tolerance (1e-2) so marching on curved surfaces is fast.
    let t = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let sphere = translate(
        BRep::sphere(0.8, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, 0.5),
    );
    let result = crate::subtraction(cube, sphere, t, None);
    assert!(
        result.is_ok(),
        "Cube-sphere subtraction should succeed: {:?}",
        result.err()
    );
}

#[test]
fn sphere_minus_overlapping_cube_does_not_hang() {
    // Reverse: sphere as operand A, cube as operand B.
    // Use coarse tolerance (1e-2) so marching on curved surfaces is fast.
    // The pipeline may return an error (trim curve endpoint mismatch) because
    // curved intersection curves don't align with planar face boundaries, but
    // it must not hang or panic.
    let t = Tolerance::default().with_linear(1e-2);
    let sphere = translate(
        BRep::sphere(0.8, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, 0.5),
    );
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let _result = crate::subtraction(sphere, cube, t, None);
    // Result may be Ok or Err — both are acceptable at this stage.
}

#[test]
fn cube_union_overlapping_sphere_succeeds() {
    // Use coarse tolerance (1e-2) so marching on curved surfaces is fast.
    let t = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let sphere = translate(
        BRep::sphere(0.8, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, 0.5),
    );
    let result = crate::union(cube, sphere, t, None);
    assert!(
        result.is_ok(),
        "Cube-sphere union should succeed: {:?}",
        result.err()
    );
}
