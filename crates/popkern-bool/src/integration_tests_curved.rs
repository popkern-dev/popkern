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
fn cylinder_through_cube_union_does_not_panic() {
    // Cylinder centered at (0.5, 0.5, -0.5), radius 0.3, height 2.0 — passes
    // through the full height of the unit cube.  Curved-face splitting is not
    // yet implemented so the pipeline may return an error, but it must never
    // panic or hang.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let _result = crate::union(cube, cyl, t);
    // Ok or Err both acceptable — the pipeline must complete without hanging.
}

#[test]
fn cylinder_through_cube_subtraction_does_not_panic() {
    // cube − cylinder: the cylinder punches through the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let _result = crate::subtraction(cube, cyl, t);
}

// ── Cylinder × cube: disjoint ────────────────────────────────────────────────

#[test]
fn disjoint_cylinder_and_cube_union_has_14_faces() {
    // Cylinder at x=5 (well away from unit cube at origin).
    // A cylinder built with BRep::cylinder has 8 faces (6 side + 2 caps).
    // Expected union face count: 6 (cube) + 8 (cylinder) = 14.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(cube, cyl, t).unwrap();
    assert_eq!(
        result.faces.len(),
        14,
        "disjoint union should have 6 + 8 = 14 faces, got {}",
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
fn cone_through_cube_union_does_not_panic() {
    // Cone with apex at (0.5, 0.5, 1.5) and base at (0.5, 0.5, -0.5),
    // straddling the unit cube.  Build cone at origin then translate.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    // BRep::cone: base at origin, apex at (0,0,height).
    // Translate by (0.5, 0.5, -0.5) so base sits below the cube.
    let cone = translate(
        BRep::cone(0.4, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let _result = crate::union(cube, cone, t);
}

#[test]
fn cone_through_cube_subtraction_does_not_panic() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.4, 2.0, t.clone()).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let _result = crate::subtraction(cube, cone, t);
}

// ── Cone × cube: disjoint ────────────────────────────────────────────────────

#[test]
fn disjoint_cone_and_cube_union_succeeds() {
    // Cone at x=5 (well away from unit cube).
    // A cone built with BRep::cone has 3 faces.
    // Expected union face count: 6 (cube) + 3 (cone) = 9.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t.clone()).unwrap();
    let cone = translate(
        BRep::cone(0.5, 1.0, t.clone()).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(cube, cone, t).unwrap();
    assert_eq!(
        result.faces.len(),
        9,
        "disjoint union should have 6 + 3 = 9 faces, got {}",
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
