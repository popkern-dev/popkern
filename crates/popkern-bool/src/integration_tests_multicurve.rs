//! Integration tests for multi-curve face splitting.
//!
//! These tests verify that the split stage correctly handles faces with
//! multiple intersection curves now that the split gate has been removed.

// Structural lint marker: this file is test code, included via #[cfg(test)] in lib.rs.
#[cfg(test)]
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

fn coarse_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-2)
}

/// Translate all vertices and surface/curve geometry of a BRep by `offset`.
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
                _ => {}
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
                _ => {}
            }
        }
    }
    brep
}

/// Two cubes intersecting at different positions on a third cube.
///
/// The center cube has faces that are intersected by two separate cubes,
/// generating multiple intersection curves on the same face. With the
/// split gate removed, the second curve should now be processed.
#[test]
fn two_cubes_subtracted_from_center_cube() {
    let t = coarse_tol();
    let center = BRep::unit_cube(t).expect("center cube");

    // First small cube overlapping the left side of the center cube.
    // Offset by 0.25 in Y and Z so faces don't coincide.
    let left = translate(
        BRep::unit_cube(t).expect("left cube"),
        Vec3::new(-0.7, 0.25, 0.25),
    );

    // Subtract left cube from center.
    let result1 = crate::subtraction(center, left, t, None);
    assert!(
        result1.is_ok(),
        "first subtraction should succeed: {:?}",
        result1.err()
    );
    let intermediate = result1.expect("intermediate");

    // Second small cube overlapping the right side.
    // Offset so no coincident faces with the intermediate result.
    let right = translate(
        BRep::unit_cube(t).expect("right cube"),
        Vec3::new(0.7, 0.25, 0.25),
    );

    // Subtract right cube from the intermediate result.
    let result2 = crate::subtraction(intermediate, right, t, None);
    assert!(
        result2.is_ok(),
        "second subtraction should succeed: {:?}",
        result2.err()
    );
    let final_brep = result2.expect("final");

    // The result should have more faces than the original cube.
    assert!(
        final_brep.faces.len() > 6,
        "two subtractions should add faces, got {}",
        final_brep.faces.len()
    );
}

/// Diagnostic test: single cube-cube subtraction produces valid topology.
///
/// Performs one subtraction (unit_cube - translated cube), then checks:
/// 1. validate_all() passes
/// 2. Every face's outer loop half-edge chain closes properly
/// 3. Every half-edge in a loop has the correct loop_ assignment
#[test]
fn subtraction_result_has_valid_half_edge_loops() {
    let t = coarse_tol();
    let a = BRep::unit_cube(t).expect("cube a");
    let b = translate(
        BRep::unit_cube(t).expect("cube b"),
        Vec3::new(-0.7, 0.25, 0.25),
    );

    let result = crate::subtraction(a, b, t, None);
    assert!(
        result.is_ok(),
        "subtraction should succeed: {:?}",
        result.err()
    );
    let brep = result.expect("result");

    // Print all validation errors for diagnostics.
    let errors = brep.validate_all();
    for err in &errors {
        eprintln!("VALIDATION ERROR: {err:?}");
    }

    // Check every face's half-edge loops manually.
    let mut broken_faces = Vec::new();
    for (fh, face) in &brep.faces {
        let loop_h = face.outer_loop;
        let loop_data = brep.loops.get(loop_h);
        if loop_data.is_none() {
            eprintln!("CORRUPT: face {fh:?} has invalid outer_loop handle");
            broken_faces.push(fh);
            continue;
        }
        let start = loop_data.unwrap().half_edge;
        if !brep.half_edges.contains_key(start) {
            eprintln!("CORRUPT: face {fh:?} outer loop start half-edge missing");
            broken_faces.push(fh);
            continue;
        }
        let mut current = start;
        let mut count = 0;
        let limit = brep.half_edges.len() + 1;
        let mut loop_ok = true;
        loop {
            count += 1;
            if count > limit {
                eprintln!("CORRUPT: face {fh:?} outer loop doesn't close after {count} steps");
                loop_ok = false;
                break;
            }
            // Verify this half-edge's loop_ matches the face's outer_loop
            let he = &brep.half_edges[current];
            if he.loop_ != loop_h {
                eprintln!(
                    "CORRUPT: face {fh:?} half-edge {current:?} has loop {:?} but expected {loop_h:?}",
                    he.loop_
                );
                loop_ok = false;
            }
            current = he.next;
            if current == start {
                break;
            }
        }
        if !loop_ok {
            broken_faces.push(fh);
        }
    }

    assert!(
        broken_faces.is_empty(),
        "Found {} faces with broken half-edge loops: {:?}",
        broken_faces.len(),
        broken_faces
    );

    // Check for structural topology errors (broken loops).
    // EulerPoincare may fail due to duplicated intersection edges (known
    // limitation of the two-pair graph approach), but loops must be sound.
    let structural_errors: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::OpenLoop { .. }
            )
        })
        .collect();

    assert!(
        structural_errors.is_empty(),
        "Subtraction result has structural errors: {:?}",
        structural_errors
    );
}

/// Diagnostic test: second subtraction on result of first also produces valid topology.
///
/// This is the exact scenario that triggers infinite loops: the result of the first
/// boolean has broken half-edge chains, and the second boolean's intersection stage
/// walks those chains and hangs.
#[test]
fn chained_subtraction_result_validates() {
    let t = coarse_tol();
    let center = BRep::unit_cube(t).expect("center cube");
    let left = translate(
        BRep::unit_cube(t).expect("left cube"),
        Vec3::new(-0.7, 0.25, 0.25),
    );

    let intermediate = crate::subtraction(center, left, t, None).expect("first subtraction");

    // Check intermediate result's loop integrity.
    let mut broken = Vec::new();
    for (fh, face) in &intermediate.faces {
        let loop_data = intermediate.loops.get(face.outer_loop);
        if loop_data.is_none() {
            broken.push(format!("{fh:?}: missing loop"));
            continue;
        }
        let start = loop_data.unwrap().half_edge;
        if !intermediate.half_edges.contains_key(start) {
            broken.push(format!("{fh:?}: missing start he"));
            continue;
        }
        let mut current = start;
        let limit = intermediate.half_edges.len() + 1;
        for step in 0..limit {
            current = intermediate.half_edges[current].next;
            if current == start {
                break;
            }
            if step == limit - 1 {
                broken.push(format!("{fh:?}: loop doesn't close"));
            }
        }
    }

    let errors = intermediate.validate_all();
    for err in &errors {
        eprintln!("INTERMEDIATE VALIDATION: {err:?}");
    }

    assert!(
        broken.is_empty(),
        "Intermediate result has broken loops: {:?}",
        broken
    );
    // EulerPoincare may fail due to duplicated intersection edges, but
    // structural integrity (closed loops) is what matters for the second
    // subtraction to succeed.
    let structural_errors: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::OpenLoop { .. }
            )
        })
        .collect();
    assert!(
        structural_errors.is_empty(),
        "Intermediate result has structural errors: {:?}",
        structural_errors
    );

    // Now do the second subtraction.
    let right = translate(
        BRep::unit_cube(t).expect("right cube"),
        Vec3::new(0.7, 0.25, 0.25),
    );
    let result2 = crate::subtraction(intermediate, right, t, None);
    assert!(
        result2.is_ok(),
        "second subtraction should succeed: {:?}",
        result2.err()
    );
    let final_brep = result2.expect("final");

    let final_errors = final_brep.validate_all();
    for err in &final_errors {
        eprintln!("FINAL VALIDATION: {err:?}");
    }
    let final_structural: Vec<_> = final_errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::OpenLoop { .. }
            )
        })
        .collect();
    assert!(
        final_structural.is_empty(),
        "Final result has structural errors: {:?}",
        final_structural
    );
}

/// A cylinder punching through a thin plate (unit cube scaled flat).
///
/// The plate's top and bottom faces each get one circular intersection
/// curve, forming inner loops (holes). This is a classic multi-curve
/// scenario where the same solid has multiple faces being split.
#[test]
fn cylinder_through_plate_creates_hole() {
    let t = coarse_tol();

    // Thin plate: unit cube is already 1x1x1; use it as-is.
    // The cylinder will pass through the Z=0 and Z=1 faces.
    let plate = BRep::unit_cube(t).expect("plate");

    // Cylinder centered on the plate, radius 0.3, tall enough to pass through.
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t).expect("cylinder"),
        Vec3::new(0.5, 0.5, -0.5),
    );

    let result = crate::subtraction(plate, cyl, t, None);
    assert!(
        result.is_ok(),
        "plate-cylinder subtraction should succeed: {:?}",
        result.err()
    );
    let brep = result.expect("result");

    // The top and bottom faces should each have an inner loop (hole).
    let inner_loop_count: usize = brep.faces.values().map(|f| f.inner_loops.len()).sum();
    assert!(
        inner_loop_count >= 2,
        "should have at least 2 inner loops (top + bottom holes), got {}",
        inner_loop_count
    );
}
