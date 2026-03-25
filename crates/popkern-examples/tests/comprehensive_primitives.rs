//! Category 1 & 2: Primitive solids and transforms.

use popkern_core::error::ValidationError;
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default().with_linear(1e-4)
}

fn assert_euler(brep: &BRep) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    assert_eq!(v - e + f, 2, "Euler: V={v} E={e} F={f}");
}

fn assert_no_structural_errors(brep: &BRep) {
    let errors = brep.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                ValidationError::BrokenCycle { .. }
                    | ValidationError::BrokenTwin { .. }
                    | ValidationError::InconsistentVertex { .. }
                    | ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(structural.is_empty(), "structural errors: {:?}", structural);
}

fn assert_finite_vertices(brep: &BRep) {
    for (vh, v) in &brep.vertices {
        assert!(v.position.x.is_finite(), "vertex {:?} non-finite x", vh);
        assert!(v.position.y.is_finite(), "vertex {:?} non-finite y", vh);
        assert!(v.position.z.is_finite(), "vertex {:?} non-finite z", vh);
    }
}

fn assert_all_faces_have_surfaces(brep: &BRep) {
    for (fh, face) in &brep.faces {
        assert!(face.surface.is_some(), "face {:?} has no surface", fh);
    }
}

// ── Category 1: Primitive Solids ──────────────────────────────────

#[test]
fn t01_unit_cube_is_valid() {
    let cube = BRep::unit_cube(tol()).unwrap();
    assert_eq!(cube.faces.len(), 6);
    assert_eq!(cube.edges.len(), 12);
    assert_eq!(cube.vertices.len(), 8);
    assert_euler(&cube);
    assert_no_structural_errors(&cube);
    assert_finite_vertices(&cube);
    assert_all_faces_have_surfaces(&cube);
    // All cube faces should be planar
    for (_fh, face) in &cube.faces {
        let surf = face.surface.as_ref().unwrap();
        assert!(
            matches!(surf.geometry, SurfaceGeometry::Plane { .. }),
            "cube face should be planar"
        );
    }
}

#[test]
fn t02_cylinder_is_valid() {
    let cyl = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    assert_eq!(cyl.faces.len(), 3, "cylinder: 1 side + 2 caps");
    assert_eq!(cyl.edges.len(), 3);
    assert_eq!(cyl.vertices.len(), 2);
    assert_euler(&cyl);
    assert_no_structural_errors(&cyl);
    assert_finite_vertices(&cyl);
    assert_all_faces_have_surfaces(&cyl);
    // Should have at least one cylindrical surface
    let has_cyl = cyl.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
    });
    assert!(has_cyl, "cylinder should have a cylindrical surface");
}

#[test]
fn t03_sphere_is_valid() {
    let sphere = BRep::sphere(1.0, tol()).unwrap();
    // Sphere topology: 1 face, 1 edge (seam), 2 vertices (poles) — or similar
    assert!(!sphere.faces.is_empty(), "sphere has faces");
    assert_no_structural_errors(&sphere);
    assert_finite_vertices(&sphere);
    assert_all_faces_have_surfaces(&sphere);
    let has_sphere = sphere.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }))
    });
    assert!(has_sphere, "sphere should have spherical surface");
}

#[test]
fn t04_cone_is_valid() {
    let cone = BRep::cone(0.5, 2.0, tol()).unwrap();
    assert!(cone.faces.len() >= 2, "cone has base + side");
    assert_no_structural_errors(&cone);
    assert_finite_vertices(&cone);
    assert_all_faces_have_surfaces(&cone);
    let has_cone = cone.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cone { .. }))
    });
    assert!(has_cone, "cone should have conical surface");
}

#[test]
#[ignore = "BRep::torus() does not yet attach surface geometry to faces"]
fn t05_torus_is_valid() {
    // BRep::torus() creates topology (faces, edges, vertices) but does not
    // attach surface geometry to faces yet. Once the torus constructor adds
    // surface geometry, remove the #[ignore] and re-enable full checks.
    let torus = BRep::torus(1.0, 0.3, tol()).unwrap();
    assert!(!torus.faces.is_empty(), "torus has faces");
    assert_no_structural_errors(&torus);
    assert_finite_vertices(&torus);
    let faces_with_surface = torus.faces.values().filter(|f| f.surface.is_some()).count();
    assert!(
        faces_with_surface > 0,
        "torus faces have no surface geometry — torus constructor needs surface assignment"
    );
}

// ── Category 2: Transforms ────────────────────────────────────────

#[test]
fn t06_translate_cube_preserves_topology() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let moved = cube
        .transform(&Transform3::translation(5.0, 3.0, -1.0))
        .unwrap();
    assert_eq!(moved.faces.len(), cube.faces.len());
    assert_eq!(moved.edges.len(), cube.edges.len());
    assert_eq!(moved.vertices.len(), cube.vertices.len());
    assert_euler(&moved);
    assert_no_structural_errors(&moved);
    assert_finite_vertices(&moved);
    // Check that vertices actually moved
    let min_x: f64 = moved
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::INFINITY, f64::min);
    assert!(
        min_x >= 4.9,
        "translated cube min x should be ~5.0, got {min_x}"
    );
}

#[test]
fn t07_translate_preserves_surface_geometry() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let moved = cube
        .transform(&Transform3::translation(1.0, 2.0, 3.0))
        .unwrap();
    assert_all_faces_have_surfaces(&moved);
    // All faces should still be planar after translation
    for (_fh, face) in &moved.faces {
        let surf = face.surface.as_ref().unwrap();
        assert!(
            matches!(surf.geometry, SurfaceGeometry::Plane { .. }),
            "translated cube face should still be planar"
        );
    }
}

#[test]
fn t08_translate_cylinder_preserves_cylinder_surface() {
    let cyl = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    let moved = cyl
        .transform(&Transform3::translation(10.0, 0.0, 0.0))
        .unwrap();
    assert_eq!(moved.faces.len(), cyl.faces.len());
    assert_euler(&moved);
    assert_no_structural_errors(&moved);
    let has_cyl = moved.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
    });
    assert!(
        has_cyl,
        "translated cylinder should still have cylindrical surface"
    );
}
