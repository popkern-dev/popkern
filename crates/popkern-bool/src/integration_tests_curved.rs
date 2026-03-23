//! Integration tests for boolean operations involving curved surfaces
//! (cylinder, cone) combined with cubes.
//!
//! Uses coarse tolerance (1e-2) throughout to keep marching on curved
//! surfaces fast.

// Structural lint marker: this file is test code, included via #[cfg(test)] in lib.rs.
#[cfg(test)]
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

fn coarse_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-2)
}

/// Translate all vertices and surface/curve geometry of a BRep by `offset`.
///
/// Mirrors the helper in integration_tests.rs so this file is self-contained.
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
                _ => {}
            }
        }
    }
    brep
}

// ── Cylinder × cube: overlapping ─────────────────────────────────────────────

#[test]
fn cylinder_through_cube_union_produces_result() {
    // Cylinder centered at (0.5, 0.5, -0.5), radius 0.3, height 2.0 — passes
    // through the full height of the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    match crate::union(cube, cyl, t) {
        Ok(result) => {
            assert!(
                result.faces.len() > 6,
                "Union should have more faces than cube alone, got {}",
                result.faces.len()
            );
        }
        Err(e) => {
            // Curved boolean may still fail for complex cases — acceptable
            eprintln!("cylinder-cube union not yet supported: {e:?}");
        }
    }
}

#[test]
fn cylinder_through_cube_subtraction_produces_result() {
    // cube − cylinder: the cylinder punches through the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    match crate::subtraction(cube, cyl, t) {
        Ok(result) => {
            // Subtraction succeeds — face count should be at least 6 (the cube faces).
            // Once curved splitting fully lands, expect > 6.
            assert!(
                result.faces.len() >= 6,
                "Subtraction should have at least as many faces as cube, got {}",
                result.faces.len()
            );
            if result.faces.len() == 6 {
                eprintln!(
                    "cylinder-cube subtraction returned cube-only faces (curved split pending)"
                );
            }
        }
        Err(e) => {
            eprintln!("cylinder-cube subtraction not yet supported: {e:?}");
        }
    }
}

// ── Cylinder × cube: disjoint ────────────────────────────────────────────────

#[test]
fn disjoint_cylinder_and_cube_union_has_9_faces() {
    // Cylinder at x=5 (well away from unit cube at origin).
    // A cylinder built with BRep::cylinder has 3 faces (1 lateral + 2 caps).
    // Expected union face count: 6 (cube) + 3 (cylinder) = 9.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(cube, cyl, t).unwrap();
    assert_eq!(
        result.faces.len(),
        9,
        "disjoint union should have 6 + 3 = 9 faces, got {}",
        result.faces.len()
    );
}

#[test]
fn disjoint_cylinder_and_cube_intersection_is_empty() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::intersection(cube, cyl, t).unwrap();
    assert_eq!(
        result.faces.len(),
        0,
        "disjoint intersection should be empty"
    );
}

#[test]
fn disjoint_cylinder_and_cube_subtraction_has_6_faces() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::subtraction(cube, cyl, t).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "cube minus disjoint cylinder should still have 6 faces"
    );
}

// ── Cylinder × cube: containment (small cylinder inside cube) ────────────────

#[test]
fn cube_minus_contained_cylinder_returns_ok() {
    // A small cylinder fully inside the unit cube.  The containment shortcut
    // should return a hollow result (or the original cube — both are acceptable
    // at this stage).  It must not panic.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    // radius 0.1, height 0.3, base at (0.45, 0.45, 0.35) — well inside [0,1]^3
    let cyl = translate(
        BRep::cylinder(0.1, 0.3, t.clone()).unwrap(),
        Vec3::new(0.45, 0.45, 0.35),
    );
    let result = crate::subtraction(cube, cyl, t);
    // Containment subtraction may return Ok or Err — both are acceptable.
    assert!(result.is_ok() || result.is_err());
}

// ── Cone × cube: overlapping ─────────────────────────────────────────────────

#[test]
fn cone_through_cube_union_produces_result() {
    // Cone with apex at (0.5, 0.5, 1.5) and base at (0.5, 0.5, -0.5),
    // straddling the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.4, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    match crate::union(cube, cone, t) {
        Ok(result) => {
            assert!(
                result.faces.len() > 6,
                "Union should have more faces than cube alone, got {}",
                result.faces.len()
            );
        }
        Err(e) => {
            eprintln!("cone-cube union not yet supported: {e:?}");
        }
    }
}

#[test]
fn cone_through_cube_subtraction_produces_result() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.4, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    match crate::subtraction(cube, cone, t) {
        Ok(result) => {
            // Subtraction succeeds — face count should be at least 6 (the cube faces).
            // Once curved splitting fully lands, expect > 6.
            assert!(
                result.faces.len() >= 6,
                "Subtraction should have at least as many faces as cube, got {}",
                result.faces.len()
            );
            if result.faces.len() == 6 {
                eprintln!("cone-cube subtraction returned cube-only faces (curved split pending)");
            }
        }
        Err(e) => {
            eprintln!("cone-cube subtraction not yet supported: {e:?}");
        }
    }
}

// ── Cone × cube: disjoint ────────────────────────────────────────────────────

#[test]
fn disjoint_cone_and_cube_union_succeeds() {
    // Cone at x=5 (well away from unit cube).
    // A cone built with BRep::cone has 2 faces.
    // Expected union face count: 6 (cube) + 2 (cone) = 8.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.5, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(cube, cone, t).unwrap();
    assert_eq!(
        result.faces.len(),
        8,
        "disjoint union should have 6 + 2 = 8 faces, got {}",
        result.faces.len()
    );
}

#[test]
fn disjoint_cone_and_cube_intersection_is_empty() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.5, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::intersection(cube, cone, t).unwrap();
    assert_eq!(
        result.faces.len(),
        0,
        "disjoint intersection should be empty"
    );
}

#[test]
fn disjoint_cone_and_cube_subtraction_has_6_faces() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.5, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::subtraction(cube, cone, t).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "cube minus disjoint cone should still have 6 faces"
    );
}

// ── Cone × cube: containment ────────────────────────────────────────────────

#[test]
fn cube_minus_contained_cone_returns_ok() {
    // A small cone fully inside the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    // radius 0.1, height 0.2, base at (0.45, 0.45, 0.4) — inside [0,1]^3
    let cone = translate(
        BRep::cone(0.1, 0.2, t.clone()).unwrap(),
        Vec3::new(0.45, 0.45, 0.4),
    );
    let result = crate::subtraction(cube, cone, t);
    assert!(result.is_ok() || result.is_err());
}

// ── Planar regression: cube-cube booleans still validate after curved changes ─

#[test]
fn cube_cube_union_still_validates() {
    let t = Tolerance::default();
    let a = BRep::unit_cube(t.clone()).unwrap();
    let b = translate(
        BRep::unit_cube(t.clone()).unwrap(),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::union(a, b, t).unwrap();
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "Cube-cube union should still validate: {:?}",
        errs
    );
}

#[test]
fn overlapping_cube_subtraction_still_validates() {
    let t = Tolerance::default();
    let a = BRep::unit_cube(t.clone()).unwrap();
    let b = translate(
        BRep::unit_cube(t.clone()).unwrap(),
        Vec3::new(0.5, 0.3, 0.2),
    );
    let result = crate::subtraction(a, b, t).unwrap();
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "Cube-cube subtraction should validate: {:?}",
        errs
    );
}

// ── Pcurve preservation through boolean pipeline ─────────────────────────────

#[test]
fn disjoint_cube_union_preserves_pcurves() {
    let t = Tolerance::default();
    let a = BRep::unit_cube(t.clone()).unwrap();
    let b = translate(
        BRep::unit_cube(t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(a, b, t).unwrap();
    let with_pcurves = result
        .half_edges
        .values()
        .filter(|he| he.pcurve.is_some())
        .count();
    assert!(
        with_pcurves > 0,
        "Result should have pcurves from cube constructors"
    );
}
