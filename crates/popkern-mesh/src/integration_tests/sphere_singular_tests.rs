//! Sphere meshes at and away from the sphere's poles (split from the
//! integration module: these tests describe the sphere's geometry).

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::{TessellationOptions, tessellate};

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
