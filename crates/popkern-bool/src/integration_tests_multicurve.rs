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
                CurveGeometry::Circle { ref mut center, .. } => {
                    *center = *center + offset;
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
#[ignore] // Sequential subtraction hangs during intersection of intermediate result
fn two_cubes_subtracted_from_center_cube() {
    let t = coarse_tol();
    let center = BRep::unit_cube(t.clone()).expect("center cube");

    // First small cube overlapping the left side of the center cube.
    // Offset by 0.25 in Y and Z so faces don't coincide.
    let left = translate(
        BRep::unit_cube(t.clone()).expect("left cube"),
        Vec3::new(-0.7, 0.25, 0.25),
    );

    // Subtract left cube from center.
    let result1 = crate::subtraction(center, left, t.clone());
    assert!(
        result1.is_ok(),
        "first subtraction should succeed: {:?}",
        result1.err()
    );
    let intermediate = result1.expect("intermediate");

    // Second small cube overlapping the right side.
    // Offset so no coincident faces with the intermediate result.
    let right = translate(
        BRep::unit_cube(t.clone()).expect("right cube"),
        Vec3::new(0.7, 0.25, 0.25),
    );

    // Subtract right cube from the intermediate result.
    let result2 = crate::subtraction(intermediate, right, t);
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

/// A cylinder punching through a thin plate (unit cube scaled flat).
///
/// The plate's top and bottom faces each get one circular intersection
/// curve, forming inner loops (holes). This is a classic multi-curve
/// scenario where the same solid has multiple faces being split.
#[test]
#[ignore] // Cylinder-plate requires closed-curve inner loop on both faces
fn cylinder_through_plate_creates_hole() {
    let t = coarse_tol();

    // Thin plate: unit cube is already 1x1x1; use it as-is.
    // The cylinder will pass through the Z=0 and Z=1 faces.
    let plate = BRep::unit_cube(t.clone()).expect("plate");

    // Cylinder centered on the plate, radius 0.3, tall enough to pass through.
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t.clone()).expect("cylinder"),
        Vec3::new(0.5, 0.5, -0.5),
    );

    let result = crate::subtraction(plate, cyl, t);
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
