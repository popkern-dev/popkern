// Satisfy structural lint: this entire file is test code.
// #[cfg(test)] — included via `#[cfg(test)] mod tests;` in mod.rs.

use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn unit_cube_topology() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    assert_eq!(cube.vertices.len(), 8);
    assert_eq!(cube.edges.len(), 12);
    assert_eq!(cube.faces.len(), 6);
}

#[test]
fn unit_cube_passes_validation() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    assert!(cube.validate().is_ok(), "Cube should pass full validation");
}

#[test]
fn unit_cube_has_geometry() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    // Every edge should have a curve
    for edge in cube.edges.values() {
        assert!(edge.curve.is_some(), "Every edge should have geometry");
    }
    // Every face should have a surface
    for face in cube.faces.values() {
        assert!(face.surface.is_some(), "Every face should have geometry");
    }
}

#[test]
fn prism_triangle_topology() {
    let profile = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
    ];
    let prism = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .unwrap();
    // Triangular prism: V=6, E=9, F=5
    assert_eq!(prism.vertices.len(), 6);
    assert_eq!(prism.edges.len(), 9);
    assert_eq!(prism.faces.len(), 5);
}

#[test]
fn prism_passes_validation() {
    let profile = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
    ];
    let prism = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .unwrap();
    assert!(prism.validate().is_ok(), "Prism should pass validation");
}

#[test]
fn prism_too_few_vertices_fails() {
    let profile = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    let result = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    );
    assert!(result.is_err());
}

#[test]
fn prism_pentagon_topology() {
    let profile = vec![
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.309, 0.951, 0.0),
        Point3::new(-0.809, 0.588, 0.0),
        Point3::new(-0.809, -0.588, 0.0),
        Point3::new(0.309, -0.951, 0.0),
    ];
    let prism = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        2.0,
        Tolerance::default(),
    )
    .unwrap();
    // Pentagonal prism: V=10, E=15, F=7
    assert_eq!(prism.vertices.len(), 10);
    assert_eq!(prism.edges.len(), 15);
    assert_eq!(prism.faces.len(), 7);
}

// --- Cylinder tests ---

#[test]
fn cylinder_topology() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Hexagonal prism: V=12, E=18, F=8
    assert_eq!(cyl.vertices.len(), 12);
    assert_eq!(cyl.edges.len(), 18);
    assert_eq!(cyl.faces.len(), 8);
}

#[test]
fn cylinder_passes_validation() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    assert!(cyl.validate().is_ok(), "Cylinder should pass validation");
}

#[test]
fn cylinder_has_geometry() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    for edge in cyl.edges.values() {
        assert!(
            edge.curve.is_some(),
            "Every edge should have curve geometry"
        );
    }
    for face in cyl.faces.values() {
        assert!(
            face.surface.is_some(),
            "Every face should have surface geometry"
        );
    }
}

#[test]
fn cylinder_has_cylindrical_side_faces() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let mut cylinder_faces = 0;
    let mut plane_faces = 0;
    for face in cyl.faces.values() {
        if let Some(ref sof) = face.surface {
            match sof.geometry {
                SurfaceGeometry::Cylinder { .. } => cylinder_faces += 1,
                SurfaceGeometry::Plane { .. } => plane_faces += 1,
                _ => {}
            }
        }
    }
    assert_eq!(cylinder_faces, 6, "Should have 6 cylindrical side faces");
    assert_eq!(plane_faces, 2, "Should have 2 planar cap faces");
}

#[test]
fn cylinder_invalid_radius() {
    let result = BRep::cylinder(0.0, 1.0, Tolerance::default());
    assert!(result.is_err());
    let result = BRep::cylinder(-1.0, 1.0, Tolerance::default());
    assert!(result.is_err());
}

#[test]
fn cylinder_invalid_height() {
    let result = BRep::cylinder(1.0, 0.0, Tolerance::default());
    assert!(result.is_err());
    let result = BRep::cylinder(1.0, -1.0, Tolerance::default());
    assert!(result.is_err());
}

// --- Sphere tests ---

#[test]
fn sphere_topology() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    assert_eq!(sph.vertices.len(), 2);
    assert_eq!(sph.edges.len(), 3);
    assert_eq!(sph.faces.len(), 3);
}

#[test]
fn sphere_passes_validation() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    assert!(sph.validate().is_ok(), "Sphere should pass validation");
}

#[test]
fn sphere_has_edge_geometry() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for edge in sph.edges.values() {
        assert!(
            edge.curve.is_some(),
            "Every edge should have curve geometry"
        );
    }
}

#[test]
fn sphere_vertex_positions() {
    let sph = BRep::sphere(2.5, Tolerance::default()).unwrap();
    let positions: Vec<_> = sph.vertices.values().map(|v| v.position).collect();
    // One vertex at (0,0,-2.5) and one at (0,0,2.5)
    let has_south = positions
        .iter()
        .any(|p| (p.z + 2.5).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
    let has_north = positions
        .iter()
        .any(|p| (p.z - 2.5).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
    assert!(has_south, "Should have south pole vertex");
    assert!(has_north, "Should have north pole vertex");
}

#[test]
fn sphere_invalid_radius() {
    let result = BRep::sphere(0.0, Tolerance::default());
    assert!(result.is_err());
    let result = BRep::sphere(-1.0, Tolerance::default());
    assert!(result.is_err());
}

#[test]
fn sphere_faces_have_surface_geometry() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for (_, face) in &sph.faces {
        let surface = face
            .surface
            .as_ref()
            .expect("every sphere face should have surface geometry");
        match &surface.geometry {
            SurfaceGeometry::Sphere { center, radius } => {
                assert!(
                    center.x.abs() < 1e-10 && center.y.abs() < 1e-10 && center.z.abs() < 1e-10,
                    "sphere center should be at origin"
                );
                assert!(
                    (*radius - 1.0).abs() < 1e-10,
                    "sphere radius should match constructor arg"
                );
            }
            other => panic!(
                "expected SurfaceGeometry::Sphere, got {:?}",
                std::mem::discriminant(other)
            ),
        }
    }
}

#[test]
fn sphere_with_custom_radius_has_correct_surface() {
    let sph = BRep::sphere(3.7, Tolerance::default()).unwrap();
    for (_, face) in &sph.faces {
        let surface = face.surface.as_ref().unwrap();
        match &surface.geometry {
            SurfaceGeometry::Sphere { radius, .. } => {
                assert!((*radius - 3.7).abs() < 1e-10);
            }
            _ => panic!("expected SurfaceGeometry::Sphere"),
        }
    }
}

// --- Pcurve tests ---

#[test]
fn sphere_all_half_edges_have_pcurves() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for (_, he) in &sph.half_edges {
        assert!(
            he.pcurve.is_some(),
            "Every sphere half-edge should have a pcurve"
        );
    }
}

#[test]
fn sphere_pcurves_are_degree_1() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for (_, he) in &sph.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        assert_eq!(pcurve.degree, 1, "Sphere pcurves should be degree 1");
        assert_eq!(pcurve.control_points.len(), 2);
        assert_eq!(pcurve.weights.len(), 2);
    }
}

#[test]
fn sphere_pcurves_have_valid_v_range() {
    // Sphere v ranges from 0 (north) to PI (south).
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for (_, he) in &sph.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        for cp in &pcurve.control_points {
            assert!(
                cp.y >= -1e-10 && cp.y <= std::f64::consts::PI + 1e-10,
                "Sphere pcurve v should be in [0, PI], got {}",
                cp.y
            );
        }
    }
}

#[test]
fn cylinder_all_half_edges_have_pcurves() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    for (_, he) in &cyl.half_edges {
        assert!(
            he.pcurve.is_some(),
            "Every cylinder half-edge should have a pcurve"
        );
    }
}

#[test]
fn cylinder_pcurves_are_degree_1() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    for (_, he) in &cyl.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        assert_eq!(pcurve.degree, 1, "Cylinder pcurves should be degree 1");
        assert_eq!(pcurve.control_points.len(), 2);
    }
}

#[test]
fn cone_all_half_edges_have_pcurves() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    for (_, he) in &cone.half_edges {
        assert!(
            he.pcurve.is_some(),
            "Every cone half-edge should have a pcurve"
        );
    }
}

#[test]
fn cone_pcurves_are_degree_1() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    for (_, he) in &cone.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        assert_eq!(pcurve.degree, 1, "Cone pcurves should be degree 1");
        assert_eq!(pcurve.control_points.len(), 2);
    }
}

#[test]
fn cone_pcurves_have_valid_v_range() {
    // Cone v ranges from 0 (apex) to height (base).
    let height = 2.0;
    let cone = BRep::cone(1.0, height, Tolerance::default()).unwrap();
    for (_, he) in &cone.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        for cp in &pcurve.control_points {
            assert!(
                cp.y >= -1e-10 && cp.y <= height + 1e-10,
                "Cone pcurve v should be in [0, height], got {}",
                cp.y
            );
        }
    }
}

// --- Torus tests ---

#[test]
fn torus_basic_succeeds() {
    let result = BRep::torus(2.0, 0.5, Tolerance::default());
    assert!(result.is_ok(), "Torus constructor should succeed");
}
