//! End-to-end integration tests for the boolean operations pipeline.
//!
//! These tests exercise the full public API (union, intersection, subtraction)
//! through the complete pipeline: intersect → split → classify → select → merge.

// Structural lint marker: this file is test code, included via #[cfg(test)] in lib.rs.
#[cfg(test)]
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

use crate::csg::{CsgNode, csg_eval};

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
                CurveGeometry::Nurbs(ref mut n) => {
                    for cp in &mut n.control_points {
                        *cp = *cp + offset;
                    }
                }
                _ => {} // NURBS and other variants: no origin to translate
            }
        }
    }
    brep
}

/// Scale a unit cube [0,1]^3 to [-1,2]^3 so it fully contains any other unit cube at origin.
fn scale_to_enclosing(brep: BRep) -> BRep {
    // Geometry-honest enlargement: vertex translation alone leaves
    // the stored surfaces stale, and the classify hit-box gate
    // rightly rejects hits on a surface that no longer matches its
    // face.
    brep.transform(&popkern_core::geometry::Transform3::scale_uniform(3.0))
        .unwrap()
        .transform(&popkern_core::geometry::Transform3::translation(
            -1.0, -1.0, -1.0,
        ))
        .unwrap()
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

/// Handle-order-sensitive signature of a boolean result: for each
/// face in arena iteration order, the outer-loop vertex positions
/// (bit-exact) and the surface geometry's Debug rendering. Two
/// executions of the same boolean must produce identical signatures
/// or downstream consumers keyed on handles (fillet carve order,
/// provenance, diffing) see a different model each run.
fn face_order_signature(brep: &BRep) -> Vec<String> {
    let mut sig = Vec::new();
    for (fh, face) in &brep.faces {
        let mut entry = format!("{fh:?}");
        if let Some(l) = brep.loops.get(face.outer_loop) {
            let start = l.half_edge;
            let mut cursor = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let Some(he) = brep.half_edges.get(cursor) else {
                    break;
                };
                if let Some(v) = brep.vertices.get(he.origin) {
                    entry.push_str(&format!(
                        " ({:x},{:x},{:x})",
                        v.position.x.to_bits(),
                        v.position.y.to_bits(),
                        v.position.z.to_bits()
                    ));
                }
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
        if let Some(surface) = &face.surface {
            entry.push_str(&format!(" {:?}", surface.geometry));
        }
        sig.push(entry);
    }
    sig
}

/// Boolean output must be deterministic run-to-run (PK-099): the
/// intersection-graph walk once iterated a HashMap, shuffling result
/// face-handle order per execution. Each HashMap instance draws
/// fresh hash keys, so two in-process executions catch unordered
/// iteration without needing separate processes.
#[test]
fn union_output_is_deterministic_across_executions() {
    let run = || {
        let a = BRep::unit_cube(tol()).unwrap();
        let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));
        crate::union(a, b, tol(), None).expect("union should succeed")
    };
    let first = face_order_signature(&run());
    for _ in 0..3 {
        let again = face_order_signature(&run());
        assert_eq!(
            first, again,
            "boolean face-handle order / geometry assignment changed between executions"
        );
    }
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
    let inner = BRep::unit_cube(tol())
        .unwrap()
        .transform(&popkern_core::geometry::Transform3::scale_uniform(0.6))
        .unwrap()
        .transform(&popkern_core::geometry::Transform3::translation(
            0.2, 0.2, 0.2,
        ))
        .unwrap();
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
fn coincident_cubes_union_produces_valid_result() {
    // Union of two identical cubes should produce a single cube.
    // The coincident coplanar face pairs are skipped; the perpendicular
    // face pairs produce intersection curves at the shared edges.
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    if let Ok(result) = crate::union(a, b, tol(), None) {
        // The result should have a valid topology
        let v = result.vertices.len();
        let e = result.edges.len();
        let f = result.faces.len();
        eprintln!(
            "[coincident union] V={v} E={e} F={f} euler={}",
            v as i64 - e as i64 + f as i64
        );
        // Union of identical cubes = one cube
        assert_eq!(f, 6, "union of identical cubes should have 6 faces");
        assert_eq!(v as i64 - e as i64 + f as i64, 2, "Euler V-E+F should be 2");
    }
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
fn cube_minus_contained_sphere_serves_hollow_solid() {
    // Sphere fully inside cube — the containment path builds a hollow
    // solid: the cube's outer shell plus the sphere as an inverted
    // inner void shell. (This assertion replaced a vacuous
    // `is_ok() || is_err()` left from when the path was a TODO.)
    let t = Tolerance::default().with_linear(1e-3);
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.3, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let result = crate::subtraction(cube, sphere, t, None).expect("containment subtraction serves");
    assert!(
        result.validate().is_ok(),
        "hollow solid must validate: {:?}",
        result.validate().err()
    );
    assert_eq!(
        result.faces.len(),
        7,
        "6 cube faces + 1 inverted sphere void face"
    );
    assert_eq!(result.shells.len(), 2, "outer shell + void shell");
}

#[test]
fn cube_minus_overlapping_sphere_succeeds() {
    // Sphere overlaps cube — curved face splitting is skipped, but the
    // pipeline should complete without hanging or panicking.
    // Use coarse tolerance (1e-2) so marching on curved surfaces is fast.
    let t = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.8, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
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
    let sphere = translate(BRep::sphere(0.8, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let cube = BRep::unit_cube(t).unwrap();
    let _result = crate::subtraction(sphere, cube, t, None);
    // Result may be Ok or Err — both are acceptable at this stage.
}

#[test]
fn cube_union_overlapping_sphere_succeeds() {
    // Use coarse tolerance (1e-2) so marching on curved surfaces is fast.
    let t = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.8, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let result = crate::union(cube, sphere, t, None);
    assert!(
        result.is_ok(),
        "Cube-sphere union should succeed: {:?}",
        result.err()
    );
    // PK-111 census: the sphere bulges through all six cube faces, so
    // each cube face keeps exactly its four corner bites (24 faces)
    // and the sphere wall survives as ONE face carrying the cube's
    // imprint as eight triangular holes (three corner arcs each, one
    // per cube corner). The pre-PK-111 single-window clip collapsed
    // the four in-face arcs per face pair into hull arcs running
    // through out-of-face regions and the result never had this
    // partition.
    let brep = result.expect("verified above");
    assert_eq!(brep.faces.len(), 25, "24 corner bites + 1 sphere wall");
    let holed: Vec<usize> = brep
        .faces
        .values()
        .map(|f| f.inner_loops.len())
        .filter(|&n| n > 0)
        .collect();
    assert_eq!(
        holed,
        vec![8],
        "exactly one face (the sphere wall) carries the 8 imprint holes"
    );
}
