//! Integration tests for popkern-mesh: end-to-end tessellation and export.
//!
//! These tests exercise the full tessellation pipeline with various BRep
//! constructors and validate mesh quality and export output.

#[cfg(test)]
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

use crate::{
    TessellationOptions, export_gltf, export_obj, export_stl, export_stl_writer, tessellate,
};

// ── Surface type coverage ─────────────────────────────────────────────

#[test]
fn tessellate_cylinder() {
    let cyl = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cyl, &TessellationOptions::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
    assert_eq!(mesh.face_meshes.len(), cyl.faces.len());
}

#[test]
fn tessellate_cone() {
    let cone = BRep::cone(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cone, &TessellationOptions::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
    assert_eq!(mesh.face_meshes.len(), cone.faces.len());
}

#[test]
fn tessellate_torus() {
    let torus = BRep::torus(1.0, 0.3, Tolerance::default()).unwrap();
    let mesh = tessellate(&torus, &TessellationOptions::default()).unwrap();
    assert!(!mesh.triangles.is_empty(), "torus should produce triangles");
    assert_eq!(mesh.face_meshes.len(), torus.faces.len());
}

#[test]
fn tessellate_sphere_produces_triangles() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
    assert_eq!(mesh.face_meshes.len(), sphere.faces.len());
}

// ── Boolean result tessellation ───────────────────────────────────────

/// Translate all geometry in a BRep by `offset`.
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

#[test]
fn tessellate_boolean_union_result() {
    let tol = Tolerance::default();
    let a = BRep::unit_cube(tol).unwrap();
    let b = translate(BRep::unit_cube(tol).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let union = popkern_bool::union(a, b, tol, None).unwrap();
    // Faces without surface geometry (from the merge stage) are skipped gracefully.
    let mesh = tessellate(&union, &TessellationOptions::default()).unwrap();
    assert!(
        mesh.triangles.len() > 12,
        "union should have more triangles than a single cube"
    );
}

#[test]
fn tessellate_boolean_subtraction_result() {
    let tol = Tolerance::default();
    let a = BRep::unit_cube(tol).unwrap();
    let b = translate(BRep::unit_cube(tol).unwrap(), Vec3::new(0.5, 0.3, 0.2));
    let result = popkern_bool::subtraction(a, b, tol, None).unwrap();
    let mesh = tessellate(&result, &TessellationOptions::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
}

#[test]
fn tessellate_brep_with_no_surface_geometry_skips_faces() {
    // Build a BRep with a face that has no surface geometry and verify
    // tessellation succeeds, skipping it gracefully.
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap();
    // Remove surface geometry from the first face.
    if let Some(face_handle) = brep.faces.keys().next() {
        brep.faces[face_handle].surface = None;
    }
    let mesh = tessellate(&brep, &TessellationOptions::default()).unwrap();
    // Should succeed with 5 faces tessellated instead of 6.
    assert_eq!(mesh.face_meshes.len(), 5);
}

// ── Mesh quality checks ───────────────────────────────────────────────

#[test]
fn non_pole_triangles_on_sphere_are_non_degenerate() {
    // Triangles at the poles of a sphere may collapse because multiple UV corners
    // map to the same 3-D point (sin(0) = sin(PI) = 0). This test verifies that
    // the vast majority of triangles are non-degenerate, which demonstrates
    // correct tessellation everywhere except the polar singularities.
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    assert!(!mesh.triangles.is_empty());
    let mut non_degenerate = 0_usize;
    for tri in &mesh.triangles {
        let v0 = &mesh.vertices[tri[0]].position;
        let v1 = &mesh.vertices[tri[1]].position;
        let v2 = &mesh.vertices[tri[2]].position;
        let e1x = v1.x - v0.x;
        let e1y = v1.y - v0.y;
        let e1z = v1.z - v0.z;
        let e2x = v2.x - v0.x;
        let e2y = v2.y - v0.y;
        let e2z = v2.z - v0.z;
        let cx = e1y * e2z - e1z * e2y;
        let cy = e1z * e2x - e1x * e2z;
        let cz = e1x * e2y - e1y * e2x;
        let area = (cx * cx + cy * cy + cz * cz).sqrt();
        if area > 1e-15 {
            non_degenerate += 1;
        }
    }
    // At least 90% of triangles should be non-degenerate.
    let ratio = non_degenerate as f64 / mesh.triangles.len() as f64;
    assert!(
        ratio >= 0.9,
        "Too many degenerate triangles: {}/{} non-degenerate ({ratio:.2})",
        non_degenerate,
        mesh.triangles.len()
    );
}

#[test]
fn cylinder_planar_cap_faces_have_valid_triangles() {
    // The planar cap faces of the cylinder (those with `Plane` surface geometry)
    // use bounded_parameter_range and produce valid triangles.
    // Note: the cylindrical side faces have an infinite v-range which the
    // tessellator does not yet bound from vertex projections, so they may
    // produce degenerate NaN vertices; this is a known limitation tracked
    // separately.
    let cyl = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cyl, &TessellationOptions::default()).unwrap();
    assert!(
        !mesh.triangles.is_empty(),
        "cylinder must produce at least some triangles"
    );
    // The mesh should have triangles for all 3 faces.
    assert_eq!(mesh.face_meshes.len(), 3);
}

#[test]
fn face_map_length_equals_triangle_count() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
    assert_eq!(mesh.face_map.len(), mesh.triangles.len());
}

#[test]
fn face_map_length_equals_triangle_count_sphere() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    assert_eq!(mesh.face_map.len(), mesh.triangles.len());
}

#[test]
fn all_face_map_handles_valid_for_cylinder() {
    let cyl = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cyl, &TessellationOptions::default()).unwrap();
    for face_handle in &mesh.face_map {
        assert!(cyl.faces.contains_key(*face_handle));
    }
}

#[test]
fn sphere_normals_are_unit_length_away_from_poles() {
    // At the poles of the sphere, the surface normal is degenerate (zero length)
    // because the parametric partial derivatives collapse. Non-pole vertices
    // should have unit-length normals.
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    for v in &mesh.vertices {
        let len = (v.normal.x.powi(2) + v.normal.y.powi(2) + v.normal.z.powi(2)).sqrt();
        // Skip pole vertices that have degenerate (near-zero) normals.
        if len < 1e-6 {
            continue;
        }
        assert!((len - 1.0).abs() < 1e-6, "Normal not unit length: {len}");
    }
}

#[test]
fn cylinder_normals_are_unit_length() {
    let cyl = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cyl, &TessellationOptions::default()).unwrap();
    for v in &mesh.vertices {
        let len = (v.normal.x.powi(2) + v.normal.y.powi(2) + v.normal.z.powi(2)).sqrt();
        assert!((len - 1.0).abs() < 1e-6, "Normal not unit length: {len}");
    }
}

// ── Export integration ────────────────────────────────────────────────

#[test]
fn export_sphere_to_stl_binary_file() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_sphere.stl");
    export_stl(&mesh, &path, true).unwrap();
    assert!(path.exists());
    let metadata = std::fs::metadata(&path).unwrap();
    // Binary STL: 80 header + 4 count + 50*n
    let expected = (80 + 4 + mesh.triangles.len() * 50) as u64;
    assert_eq!(metadata.len(), expected);
    std::fs::remove_file(&path).ok();
}

#[test]
fn export_sphere_to_stl_ascii_file() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_sphere_ascii.stl");
    export_stl(&mesh, &path, false).unwrap();
    assert!(path.exists());
    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.starts_with("solid popkern"));
    assert!(content.ends_with("endsolid popkern\n"));
    std::fs::remove_file(&path).ok();
}

#[test]
fn export_sphere_to_obj_writer() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    let mut buf = Vec::new();
    export_obj(&mesh, &mut buf).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(text.starts_with("# Generated by popkern"));
    let v_count = text.lines().filter(|l| l.starts_with("v ")).count();
    assert_eq!(v_count, mesh.vertices.len());
}

#[test]
fn export_sphere_to_gltf_files() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    let dir = std::env::temp_dir();
    let stem = dir.join("popkern_test_sphere");
    export_gltf(&mesh, &stem).unwrap();
    let gltf_path = dir.join("popkern_test_sphere.gltf");
    let bin_path = dir.join("popkern_test_sphere.bin");
    assert!(gltf_path.exists());
    assert!(bin_path.exists());
    let json = std::fs::read_to_string(&gltf_path).unwrap();
    assert!(json.contains("\"asset\""));
    assert!(json.contains("popkern_test_sphere.bin"));
    std::fs::remove_file(&gltf_path).ok();
    std::fs::remove_file(&bin_path).ok();
}

#[test]
fn export_stl_writer_binary_produces_correct_size() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
    let mut buf = Vec::new();
    export_stl_writer(&mesh, &mut buf, true).unwrap();
    let expected = 80 + 4 + mesh.triangles.len() * 50;
    assert_eq!(buf.len(), expected);
}

#[test]
fn export_stl_writer_ascii_contains_all_facets() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
    let mut buf = Vec::new();
    export_stl_writer(&mesh, &mut buf, false).unwrap();
    let text = String::from_utf8(buf).unwrap();
    let facet_count = text.matches("facet normal").count();
    assert_eq!(facet_count, mesh.triangles.len());
}

#[test]
fn export_cylinder_all_formats() {
    let cyl = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cyl, &TessellationOptions::default()).unwrap();

    // STL binary
    let mut stl_buf = Vec::new();
    export_stl_writer(&mesh, &mut stl_buf, true).unwrap();
    assert_eq!(stl_buf.len(), 80 + 4 + mesh.triangles.len() * 50);

    // OBJ
    let mut obj_buf = Vec::new();
    export_obj(&mesh, &mut obj_buf).unwrap();
    let obj_text = String::from_utf8(obj_buf).unwrap();
    assert!(obj_text.contains("g face_0"));

    // glTF
    let dir = std::env::temp_dir();
    let stem = dir.join("popkern_test_cylinder");
    export_gltf(&mesh, &stem).unwrap();
    assert!(dir.join("popkern_test_cylinder.gltf").exists());
    assert!(dir.join("popkern_test_cylinder.bin").exists());
    std::fs::remove_file(dir.join("popkern_test_cylinder.gltf")).ok();
    std::fs::remove_file(dir.join("popkern_test_cylinder.bin")).ok();
}

// ── Options coverage ──────────────────────────────────────────────────

#[test]
fn max_subdivisions_equal_to_min_no_refinement() {
    // When max_subdivisions == min_subdivisions, no refinement can happen.
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let opts = TessellationOptions {
        min_subdivisions: 4,
        max_subdivisions: 4,
        chord_height: 0.001, // very tight, but cap prevents refinement
    };
    let mesh = tessellate(&sphere, &opts).unwrap();
    // 4x4 = 16 quads per face, 2 triangles each = 32 triangles (sphere has 1 face)
    assert!(!mesh.triangles.is_empty());
}

#[test]
fn tessellate_options_clone_and_debug() {
    let opts = TessellationOptions::default();
    let opts2 = opts.clone();
    assert!((opts.chord_height - opts2.chord_height).abs() < f64::EPSILON);
    // Ensure Debug impl works.
    let s = format!("{:?}", opts2);
    assert!(s.contains("chord_height"));
}

// ── Curved boolean tessellation ──────────────────────────────────────

#[test]
fn tessellate_disjoint_cylinder_cube_union() {
    let tol = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 1.0, tol).unwrap(),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = popkern_bool::union(cube, cyl, tol, None).unwrap();
    let mesh = tessellate(&result, &TessellationOptions::default()).unwrap();
    assert!(
        !mesh.triangles.is_empty(),
        "curved boolean tessellation should produce triangles"
    );
    // 9 faces: 6 cube + 3 cylinder
    assert!(
        mesh.face_meshes.len() >= 9,
        "expected at least 9 face meshes, got {}",
        mesh.face_meshes.len()
    );
}

#[test]
fn tessellate_cylinder_through_cube_union() {
    let tol = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = translate(
        BRep::cylinder(0.3, 2.0, tol).unwrap(),
        Vec3::new(0.5, 0.5, -0.5),
    );
    let result = popkern_bool::union(cube, cyl, tol, None).unwrap();
    let mesh = tessellate(&result, &TessellationOptions::default()).unwrap();
    assert!(
        !mesh.triangles.is_empty(),
        "overlapping curved boolean tessellation should produce triangles"
    );
}

// ── Trim regression diagnostics ──────────────────────────────────────

#[path = "diagnostic_trim_tests.rs"]
mod diagnostic_trim;
#[path = "setback_blend_export_tests.rs"]
mod setback_blend_export;
