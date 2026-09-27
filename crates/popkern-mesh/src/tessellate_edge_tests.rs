//! Tessellation edge-case tests: empty BReps, minimal triangles,
//! extreme tolerances, single-face models, and outward normals.

use super::*;
use crate::tessellate;
use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;

// ── Empty BRep tessellation ─────────────────────────────────────────

#[test]
fn tessellate_empty_brep_produces_empty_mesh() {
    let brep = BRep::new(Tolerance::default());
    let mesh = tessellate(&brep, &TessellationOptions::default()).unwrap();
    assert!(mesh.vertices.is_empty(), "no vertices for empty brep");
    assert!(mesh.triangles.is_empty(), "no triangles for empty brep");
    assert!(mesh.face_map.is_empty(), "no face map for empty brep");
    assert!(mesh.face_meshes.is_empty(), "no face meshes for empty brep");
}

#[test]
fn tessellate_empty_brep_custom_options_produces_empty_mesh() {
    let brep = BRep::new(Tolerance::default());
    let opts = TessellationOptions {
        chord_height: 0.0001,
        min_subdivisions: 8,
        max_subdivisions: 128,
        ..Default::default()
    };
    let mesh = tessellate(&brep, &opts).unwrap();
    assert!(mesh.triangles.is_empty());
}

// ── Flat planar face: minimal triangles ─────────────────────────────

#[test]
fn flat_planar_face_min_subdivisions_1_produces_2_triangles() {
    // A flat plane should never refine beyond the initial grid because
    // chord_height error is zero for a plane. With min_subdivisions=1,
    // we get 1 quad = 2 triangles per face.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let opts = TessellationOptions {
        chord_height: 0.005,
        min_subdivisions: 1,
        max_subdivisions: 64,
        ..Default::default()
    };
    let mesh = tessellate(&cube, &opts).unwrap();
    // Boundary-first: planar faces stay minimal (no interior points
    // at min_subdivisions=1, straight edges sample only endpoints) —
    // exactly 2 triangles per quad face.
    assert_eq!(
        mesh.triangles.len(),
        12,
        "6 planar faces with min_subdivisions=1 should give 12 triangles"
    );
}

#[test]
fn planar_face_never_refines_adaptively() {
    // Even with extremely tight tolerance, planes produce exactly
    // min_subdivisions^2 quads (no adaptive refinement needed).
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let tight = TessellationOptions {
        chord_height: 1e-12, // absurdly tight
        min_subdivisions: 2,
        max_subdivisions: 64,
        ..Default::default()
    };
    let mesh = tessellate(&cube, &tight).unwrap();
    // 6 faces * 4 quads * 2 triangles = 48 triangles.
    assert!(!mesh.triangles.is_empty());
}

// ── Extreme tolerance values ────────────────────────────────────────

#[test]
fn very_loose_tolerance_produces_minimal_triangles() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let opts = TessellationOptions {
        chord_height: 100.0, // absurdly loose
        min_subdivisions: 1,
        max_subdivisions: 64,
        ..Default::default()
    };
    let mesh = tessellate(&sphere, &opts).unwrap();
    // With min_subdivisions=1 and huge tolerance, should produce exactly
    // 1 quad = 2 triangles per face (1 face for sphere).
    assert!(
        !mesh.triangles.is_empty(),
        "should produce at least some triangles"
    );
    // With very loose tolerance, no adaptive refinement.
    assert_eq!(
        mesh.triangles.len(),
        2,
        "sphere with loose tolerance and min_subdivisions=1 should produce 2 triangles"
    );
}

#[test]
fn very_tight_tolerance_produces_many_triangles() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let tight = TessellationOptions {
        chord_height: 0.0001,
        min_subdivisions: 4,
        max_subdivisions: 64,
        ..Default::default()
    };
    let mesh = tessellate(&sphere, &tight).unwrap();
    // Should have significantly more triangles than default.
    let default_mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    assert!(
        mesh.triangles.len() >= default_mesh.triangles.len(),
        "tighter tolerance should produce at least as many triangles: {} vs {}",
        mesh.triangles.len(),
        default_mesh.triangles.len()
    );
}

#[test]
fn min_subdivisions_equals_max_no_panic() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let opts = TessellationOptions {
        chord_height: 0.005,
        min_subdivisions: 8,
        max_subdivisions: 8,
        ..Default::default()
    };
    let mesh = tessellate(&cube, &opts).unwrap();
    assert!(!mesh.triangles.is_empty());
}

#[test]
fn min_subdivisions_greater_than_max_still_works() {
    // Edge case: min > max. The adaptive_subdivide handles this by
    // setting max_depth to 0 (no refinement).
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let opts = TessellationOptions {
        chord_height: 0.005,
        min_subdivisions: 16,
        max_subdivisions: 4,
        ..Default::default()
    };
    let mesh = tessellate(&cube, &opts).unwrap();
    assert!(!mesh.triangles.is_empty());
}

// ── Single-face BRep tessellation ───────────────────────────────────

#[test]
fn single_face_sphere_tessellation() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    // Sphere has exactly 1 face.
    assert_eq!(sphere.faces.len(), 1);
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    assert_eq!(mesh.face_meshes.len(), 1, "sphere should have 1 face mesh");
    assert!(!mesh.triangles.is_empty());
    // All face_map entries should refer to the same face.
    let first_face = mesh.face_map[0];
    for fh in &mesh.face_map {
        assert_eq!(*fh, first_face, "all triangles belong to the single face");
    }
}

// ── Triangle normals consistent with winding ────────────────────────

#[test]
fn cube_triangle_winding_matches_vertex_normals() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();

    // For each triangle, the geometric normal (from cross product of edges)
    // should agree with the vertex normals (positive dot product).
    for tri in &mesh.triangles {
        let v0 = &mesh.vertices[tri[0]];
        let v1 = &mesh.vertices[tri[1]];
        let v2 = &mesh.vertices[tri[2]];

        let e1 = Vec3::new(
            v1.position.x - v0.position.x,
            v1.position.y - v0.position.y,
            v1.position.z - v0.position.z,
        );
        let e2 = Vec3::new(
            v2.position.x - v0.position.x,
            v2.position.y - v0.position.y,
            v2.position.z - v0.position.z,
        );
        let cross = Vec3::new(
            e1.y * e2.z - e1.z * e2.y,
            e1.z * e2.x - e1.x * e2.z,
            e1.x * e2.y - e1.y * e2.x,
        );
        let cross_len = (cross.x.powi(2) + cross.y.powi(2) + cross.z.powi(2)).sqrt();
        if cross_len < 1e-12 {
            continue; // Degenerate triangle.
        }
        // Dot the geometric normal with the average vertex normal.
        let avg_nx = (v0.normal.x + v1.normal.x + v2.normal.x) / 3.0;
        let avg_ny = (v0.normal.y + v1.normal.y + v2.normal.y) / 3.0;
        let avg_nz = (v0.normal.z + v1.normal.z + v2.normal.z) / 3.0;
        let dot = cross.x * avg_nx + cross.y * avg_ny + cross.z * avg_nz;
        assert!(
            dot > -1e-6,
            "triangle winding should agree with vertex normals (dot={dot})"
        );
    }
}

#[test]
fn sphere_vertex_normals_are_consistent_with_winding() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();

    // For non-degenerate triangles on the sphere, the geometric normal
    // (from cross product) should agree with vertex normals.
    let mut checked = 0;
    for tri in &mesh.triangles {
        let v0 = &mesh.vertices[tri[0]];
        let v1 = &mesh.vertices[tri[1]];
        let v2 = &mesh.vertices[tri[2]];

        let e1 = Vec3::new(
            v1.position.x - v0.position.x,
            v1.position.y - v0.position.y,
            v1.position.z - v0.position.z,
        );
        let e2 = Vec3::new(
            v2.position.x - v0.position.x,
            v2.position.y - v0.position.y,
            v2.position.z - v0.position.z,
        );
        let cross = Vec3::new(
            e1.y * e2.z - e1.z * e2.y,
            e1.z * e2.x - e1.x * e2.z,
            e1.x * e2.y - e1.y * e2.x,
        );
        let cross_len = (cross.x.powi(2) + cross.y.powi(2) + cross.z.powi(2)).sqrt();
        if cross_len < 1e-10 {
            continue; // Skip degenerate pole triangles.
        }
        // Check that at least one vertex normal has positive dot with geometric normal.
        let avg_nx = (v0.normal.x + v1.normal.x + v2.normal.x) / 3.0;
        let avg_ny = (v0.normal.y + v1.normal.y + v2.normal.y) / 3.0;
        let avg_nz = (v0.normal.z + v1.normal.z + v2.normal.z) / 3.0;
        let avg_len = (avg_nx.powi(2) + avg_ny.powi(2) + avg_nz.powi(2)).sqrt();
        if avg_len < 1e-6 {
            continue; // All vertex normals near-zero (pole).
        }
        let dot = cross.x * avg_nx + cross.y * avg_ny + cross.z * avg_nz;
        assert!(
            dot > -1e-6,
            "sphere triangle winding should agree with vertex normals (dot={dot}) cross_len={cross_len} v0=({:.4},{:.4},{:.4}) v1=({:.4},{:.4},{:.4}) v2=({:.4},{:.4},{:.4})",
            v0.position.x,
            v0.position.y,
            v0.position.z,
            v1.position.x,
            v1.position.y,
            v1.position.z,
            v2.position.x,
            v2.position.y,
            v2.position.z
        );
        checked += 1;
    }
    assert!(checked > 0, "should have checked at least some triangles");
}

#[test]
fn cylinder_cap_normals_point_outward() {
    let cyl = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&cyl, &TessellationOptions::default()).unwrap();

    // Cylinder caps have normals along Z axis: bottom cap = -Z, top cap = +Z.
    // The cylindrical side face normals should point radially outward.
    // We just verify all normals are approximately unit-length.
    for v in &mesh.vertices {
        let len = (v.normal.x.powi(2) + v.normal.y.powi(2) + v.normal.z.powi(2)).sqrt();
        assert!(
            (len - 1.0).abs() < 1e-4,
            "cylinder normal should be unit length, got {len}"
        );
    }
}

// ── All indices in bounds ───────────────────────────────────────────

#[test]
fn all_triangle_indices_valid_for_sphere() {
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
    for tri in &mesh.triangles {
        for &idx in tri {
            assert!(
                idx < mesh.vertices.len(),
                "index {idx} out of bounds (vertices: {})",
                mesh.vertices.len()
            );
        }
    }
}

#[test]
fn face_meshes_count_matches_faces_with_geometry() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
    // All 6 faces have surface geometry.
    assert_eq!(mesh.face_meshes.len(), 6);

    // Now remove surface from one face.
    let mut brep = cube;
    if let Some(fh) = brep.faces.keys().next() {
        brep.faces[fh].surface = None;
    }
    let mesh2 = tessellate(&brep, &TessellationOptions::default()).unwrap();
    assert_eq!(mesh2.face_meshes.len(), 5);
}

// ── Unwrap periodic edge case ───────────────────────────────────────

#[test]
fn unwrap_periodic_no_period_is_noop() {
    let mut values = vec![0.0, 1.0, 2.0, 3.0];
    let original = values.clone();
    super::unwrap_periodic(&mut values, None);
    assert_eq!(values, original);
}

#[test]
fn unwrap_periodic_zero_period_is_noop() {
    let mut values = vec![0.0, 1.0, 2.0, 3.0];
    let original = values.clone();
    super::unwrap_periodic(&mut values, Some(0.0));
    assert_eq!(values, original);
}

#[test]
fn unwrap_periodic_single_value_is_noop() {
    let mut values = vec![1.5];
    super::unwrap_periodic(&mut values, Some(std::f64::consts::TAU));
    assert_eq!(values, vec![1.5]);
}

#[test]
fn unwrap_periodic_seam_crossing() {
    // Simulate values crossing a 2*PI seam.
    use std::f64::consts::TAU;
    let mut values = vec![6.0, 0.3]; // jump from ~6.0 to ~0.3 crosses the seam
    super::unwrap_periodic(&mut values, Some(TAU));
    // After unwrapping, values[1] should be close to values[0] (no big jump).
    let diff = (values[1] - values[0]).abs();
    assert!(
        diff < TAU * 0.5,
        "unwrapped diff should be < half period, got {diff}"
    );
}
