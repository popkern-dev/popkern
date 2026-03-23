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
//
// The full pipeline (intersect → split → classify → select → merge) does
// not yet handle all overlapping-face cases. These tests verify the
// pipeline does not panic and returns a well-typed result. When the split
// stage is hardened, these should be tightened to assert `is_ok()` and
// check face counts.

#[test]
fn overlapping_cubes_union_succeeds() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let result = crate::union(a, b, tol());
    assert!(result.is_ok(), "Union should succeed: {:?}", result.err());
    let brep = result.unwrap();
    eprintln!(
        "UNION: {} faces, {} edges, {} verts, {} half_edges",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len(),
        brep.half_edges.len()
    );
    let errs = brep.validate_all();
    for e in &errs {
        eprintln!("  VALIDATION: {e}");
    }
    assert!(
        brep.faces.len() > 6,
        "Union of overlapping cubes should have more than 6 faces, got {}",
        brep.faces.len()
    );
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
    eprintln!(
        "INTERSECTION: {} faces, {} edges, {} verts, {} half_edges",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len(),
        brep.half_edges.len()
    );
    let errs = brep.validate_all();
    for e in &errs {
        eprintln!("  VALIDATION: {e}");
    }
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
    assert!(
        !brep.faces.is_empty(),
        "Subtraction should produce at least 1 face"
    );
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

#[test]
fn diagnose_union_pipeline() {
    use crate::classify::classify_fragments;
    use crate::error::BoolOperation;
    use crate::intersect::find_intersections;
    use crate::select::{Operand, select_faces};
    use crate::split::split_along_curves;

    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));

    // Stage 1
    let isect = find_intersections(&a, &b, &tol()).unwrap();
    eprintln!("STAGE1: {} intersection records", isect.records.len());
    for (i, r) in isect.records.iter().enumerate() {
        eprintln!(
            "  record {i}: face_a={:?}, face_b={:?}, closed={}",
            r.face_a, r.face_b, r.curve.is_closed
        );
    }

    // Stage 2
    let split = split_along_curves(a, b, &isect.records, &tol()).unwrap();
    eprintln!(
        "STAGE2: A has {} faces, B has {} faces",
        split.a.faces.len(),
        split.b.faces.len()
    );
    eprintln!(
        "  A intersection_edges: {}",
        split.fragments_a.intersection_edges.len()
    );
    eprintln!(
        "  B intersection_edges: {}",
        split.fragments_b.intersection_edges.len()
    );

    // Stage 3
    let class_a = classify_fragments(
        &split.a,
        &split.b,
        &isect.aabbs_b,
        &split.fragments_a.intersection_edges,
        &tol(),
    )
    .unwrap();
    let class_b = classify_fragments(
        &split.b,
        &split.a,
        &isect.aabbs_a,
        &split.fragments_b.intersection_edges,
        &tol(),
    )
    .unwrap();
    eprintln!("STAGE3: A classifications:");
    for (f, c) in &class_a {
        eprintln!("  {:?} => {:?}", f, c);
    }
    eprintln!("STAGE3: B classifications:");
    for (f, c) in &class_b {
        eprintln!("  {:?} => {:?}", f, c);
    }

    // Stage 4
    let mut classifications = Vec::new();
    for (f, c) in &class_a {
        classifications.push((*f, Operand::A, *c));
    }
    for (f, c) in &class_b {
        classifications.push((*f, Operand::B, *c));
    }
    let selection = select_faces(BoolOperation::Union, &classifications);
    eprintln!(
        "STAGE4: keep_a={}, keep_b={}, flip_a={}, flip_b={}",
        selection.keep_a.len(),
        selection.keep_b.len(),
        selection.flip_a.len(),
        selection.flip_b.len()
    );
    eprintln!("  A total={}, A kept:", split.a.faces.len());
    for f in split.a.faces.keys() {
        let kept = selection.keep_a.contains(&f);
        let class = class_a.get(&f);
        eprintln!("    {:?} class={:?} kept={}", f, class, kept);
    }
    eprintln!("  B total={}, B kept:", split.b.faces.len());
    for f in split.b.faces.keys() {
        let kept = selection.keep_b.contains(&f);
        let class = class_b.get(&f);
        eprintln!("    {:?} class={:?} kept={}", f, class, kept);
    }
}
