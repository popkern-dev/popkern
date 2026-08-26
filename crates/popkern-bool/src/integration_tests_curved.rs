//! Integration tests for boolean operations involving curved surfaces
//! (cylinder, cone, sphere) combined with cubes.
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

// ── Cylinder × cube: overlapping ─────────────────────────────────────────────

#[test]
fn cylinder_through_cube_union_produces_result() {
    // Cylinder centered at (0.5, 0.5, -0.5), radius 0.3, height 2.0 — passes
    // through the full height of the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let result = crate::union(cube, cyl, t, None).expect("cylinder-cube union should succeed");
    assert!(
        result.faces.len() > 6,
        "Union should have more faces than cube alone, got {}",
        result.faces.len()
    );
}

#[test]
fn cylinder_through_cube_subtraction_produces_result() {
    // cube - cylinder: the cylinder punches through the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, t).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let result =
        crate::subtraction(cube, cyl, t, None).expect("cylinder-cube subtraction should succeed");
    // Subtraction: cube faces with inner loops (cylinder holes).
    // The cylinder lateral face may or may not be included depending on
    // sample point selection.
    assert!(
        result.faces.len() >= 6,
        "Subtraction should have at least cube faces, got {}",
        result.faces.len()
    );
    // At least some cube faces should have inner loops (holes).
    let inner_loop_count: usize = result.faces.values().map(|f| f.inner_loops.len()).sum();
    assert!(
        inner_loop_count >= 2,
        "Should have at least 2 inner loops (cylinder holes), got {}",
        inner_loop_count
    );
}

// ── Cylinder × cube: disjoint ────────────────────────────────────────────────

#[test]
fn disjoint_cylinder_and_cube_union_has_9_faces() {
    // Cylinder at x=5 (well away from unit cube at origin).
    // A cylinder built with BRep::cylinder has 3 faces (1 lateral + 2 caps).
    // Expected union face count: 6 (cube) + 3 (cylinder) = 9.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::union(cube, cyl, t, None).unwrap();
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
    let cube = BRep::unit_cube(t).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::intersection(cube, cyl, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        0,
        "disjoint intersection should be empty"
    );
}

#[test]
fn disjoint_cylinder_and_cube_subtraction_has_6_faces() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, t).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = crate::subtraction(cube, cyl, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "cube minus disjoint cylinder should still have 6 faces"
    );
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "disjoint cylinder-cube subtraction should validate: {:?}",
        errs
    );
}

// ── Cylinder × cube: containment (small cylinder inside cube) ────────────────

#[test]
fn cube_minus_contained_cylinder_returns_ok() {
    // A small cylinder fully inside the unit cube. The containment shortcut
    // returns the original cube (hollow subtraction not yet implemented).
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    // radius 0.1, height 0.3, base at (0.45, 0.45, 0.35) — well inside [0,1]^3
    let cyl = translate(
        BRep::cylinder(0.1, 0.3, t).unwrap(),
        Vec3::new(0.45, 0.45, 0.35),
    );
    let result = crate::subtraction(cube, cyl, t, None)
        .expect("contained cylinder subtraction should succeed");
    assert!(
        result.faces.len() >= 6,
        "result should have at least cube faces, got {}",
        result.faces.len()
    );
}

// ── Cone × cube: overlapping ─────────────────────────────────────────────────

#[test]
fn cone_through_cube_union_produces_result() {
    // Cone with apex at (0.5, 0.5, 1.5) and base at (0.5, 0.5, -0.5),
    // straddling the unit cube.
    //
    // PK-343: this union had ALWAYS produced a torn cone-top loop
    // (the walk refused the apex turnaround — the seam's twins in one
    // loop — and the election filed the small top region as a "hole"
    // of the disjoint bottom frustum). Served silently until the
    // PK-342 closure gate; now serves whole: 6 cube faces (top and
    // bottom carrying the cone holes) + frustum lateral + base cap +
    // apex piece.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cone = translate(BRep::cone(0.4, 2.0, t).unwrap(), Vec3::new(0.5, 0.5, -0.5));
    let result = crate::union(cube, cone, t, None).expect("cone-cube union serves");
    assert_eq!(
        result.faces.len(),
        9,
        "cube 6 + frustum + base cap + apex piece"
    );
    assert!(result.validate_all().is_empty());
}

#[test]
fn cone_through_cube_subtraction_produces_result() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cone = translate(BRep::cone(0.4, 2.0, t).unwrap(), Vec3::new(0.5, 0.5, -0.5));
    let result =
        crate::subtraction(cube, cone, t, None).expect("cone-cube subtraction should succeed");
    assert!(
        result.faces.len() >= 6,
        "Subtraction should have at least cube faces, got {}",
        result.faces.len()
    );
    let inner_loop_count: usize = result.faces.values().map(|f| f.inner_loops.len()).sum();
    assert!(
        inner_loop_count >= 2,
        "Should have at least 2 inner loops (cone holes), got {}",
        inner_loop_count
    );
}

// ── Cone × cube: disjoint ────────────────────────────────────────────────────

#[test]
fn disjoint_cone_and_cube_union_succeeds() {
    // Cone at x=5 (well away from unit cube).
    // A cone built with BRep::cone has 2 faces.
    // Expected union face count: 6 (cube) + 2 (cone) = 8.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cone = translate(BRep::cone(0.5, 1.0, t).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::union(cube, cone, t, None).unwrap();
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
    let cube = BRep::unit_cube(t).unwrap();
    let cone = translate(BRep::cone(0.5, 1.0, t).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::intersection(cube, cone, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        0,
        "disjoint intersection should be empty"
    );
}

#[test]
fn disjoint_cone_and_cube_subtraction_has_6_faces() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let cone = translate(BRep::cone(0.5, 1.0, t).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::subtraction(cube, cone, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "cube minus disjoint cone should still have 6 faces"
    );
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "disjoint cone-cube subtraction should validate: {:?}",
        errs
    );
}

// ── Cone × cube: containment ────────────────────────────────────────────────

#[test]
fn cube_minus_contained_cone_returns_ok() {
    // A small cone fully inside the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    // radius 0.1, height 0.2, base at (0.45, 0.45, 0.4) — inside [0,1]^3
    let cone = translate(BRep::cone(0.1, 0.2, t).unwrap(), Vec3::new(0.45, 0.45, 0.4));
    let result =
        crate::subtraction(cube, cone, t, None).expect("contained cone subtraction should succeed");
    assert!(
        result.faces.len() >= 6,
        "result should have at least cube faces, got {}",
        result.faces.len()
    );
}

// ── Sphere × cube: disjoint ─────────────────────────────────────────────────

#[test]
fn disjoint_sphere_and_cube_union_has_7_faces() {
    // Sphere at x=5 (well away from unit cube at origin).
    // A sphere built with BRep::sphere has 1 face.
    // Expected union face count: 6 (cube) + 1 (sphere) = 7.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.4, t).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::union(cube, sphere, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        7,
        "disjoint union should have 6 + 1 = 7 faces, got {}",
        result.faces.len()
    );
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "disjoint sphere-cube union should validate: {:?}",
        errs
    );
}

#[test]
fn disjoint_sphere_and_cube_intersection_is_empty() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.4, t).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::intersection(cube, sphere, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        0,
        "disjoint intersection should be empty"
    );
}

#[test]
fn disjoint_sphere_and_cube_subtraction_has_6_faces() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.4, t).unwrap(), Vec3::new(5.0, 0.0, 0.0));
    let result = crate::subtraction(cube, sphere, t, None).unwrap();
    assert_eq!(
        result.faces.len(),
        6,
        "cube minus disjoint sphere should still have 6 faces"
    );
    let errs = result.validate_all();
    assert!(
        errs.is_empty(),
        "disjoint sphere-cube subtraction should validate: {:?}",
        errs
    );
}

// ── Sphere × cube: overlapping ──────────────────────────────────────────────

#[test]
fn sphere_through_cube_union_produces_result() {
    // Sphere centered at (0.5, 0.5, 0.5) with radius 0.8 — protrudes from
    // all faces of the unit cube (sphere encloses cube center).
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.8, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let result = crate::union(cube, sphere, t, None).expect("sphere-cube union should succeed");
    // For r=0.8 sphere enclosing cube, the union is the sphere.
    // The empty-selection fallback detects all-A-inside and returns B (sphere).
    assert!(
        !result.faces.is_empty(),
        "Union should produce at least 1 face (the sphere), got {}",
        result.faces.len()
    );
}

#[test]
fn sphere_through_cube_subtraction_produces_result() {
    // r=0.8 sphere encloses cube — subtraction leaves sphere shell only.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.8, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let result =
        crate::subtraction(cube, sphere, t, None).expect("sphere-cube subtraction should succeed");
    assert!(
        !result.faces.is_empty(),
        "Subtraction should produce faces, got {}",
        result.faces.len()
    );
}

#[test]
fn sphere_medium_cube_subtraction_has_inner_loops() {
    // r=0.6 sphere centered in cube: intersection circles fit inside each face.
    // Expected: 6 cube faces with inner loops + 1 sphere face = 7 faces.
    //
    // DEFAULT tolerance, not coarse: at 1e-2 this subtraction is
    // broken in every era (measured: volume 1.098 vs the closed-form
    // 0.202 before PK-341, 0.895 after — both with validator errors
    // and hundreds of non-manifold mesh edges) and the historical
    // loop-count pin was satisfied by ACCIDENT of the broken output's
    // debris. At the default tolerance the result is the true shape
    // (volume 0.212 vs 0.202 closed form, validate_all clean) and is
    // identical before/after the PK-341 walker fixes.
    let t = Tolerance::default();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.6, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let result = crate::subtraction(cube, sphere, t, None).expect("medium sphere-cube subtraction");
    assert!(
        result.faces.len() >= 7,
        "Should have 6 cube faces + sphere face, got {}",
        result.faces.len()
    );
    let inner_loops: usize = result.faces.values().map(|f| f.inner_loops.len()).sum();
    assert!(
        inner_loops >= 6,
        "Should have 6 inner loops (sphere holes on each cube face), got {}",
        inner_loops
    );
    assert!(
        result.validate_all().is_empty(),
        "default-tolerance sphere-cube must serve validator-clean"
    );
}

#[test]
fn sphere_through_cube_intersection_produces_result() {
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.8, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let _result = crate::intersection(cube, sphere, t, None)
        .expect("sphere-cube intersection should succeed");
}

// ── Sphere × cube: containment ──────────────────────────────────────────────

#[test]
fn cube_minus_contained_sphere_returns_ok() {
    // Small sphere fully inside the unit cube.
    let t = coarse_tol();
    let cube = BRep::unit_cube(t).unwrap();
    let sphere = translate(BRep::sphere(0.2, t).unwrap(), Vec3::new(0.5, 0.5, 0.5));
    let result = crate::subtraction(cube, sphere, t, None)
        .expect("contained sphere subtraction should succeed");
    assert!(
        result.faces.len() >= 6,
        "result should have at least cube faces, got {}",
        result.faces.len()
    );
}

// Planar regression and pcurve preservation tests are in the
// `integration_tests_curved_regression` sibling module.
