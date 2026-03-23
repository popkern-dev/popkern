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
    // Proper 3-face cylinder: V=2, E=3, F=3
    assert_eq!(cyl.vertices.len(), 2);
    assert_eq!(cyl.edges.len(), 3);
    assert_eq!(cyl.faces.len(), 3);
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
fn cylinder_has_surface_types() {
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
    assert_eq!(cylinder_faces, 1, "Should have 1 cylindrical face");
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

#[test]
fn cylinder_seam_edge_both_half_edges_on_lateral_face() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // The seam edge has both half-edges on the same (lateral) face.
    // Find it: the edge whose two half-edges share the same face.
    let mut found_seam = false;
    for edge in cyl.edges.values() {
        let he = &cyl.half_edges[edge.half_edge];
        let twin = &cyl.half_edges[he.twin];
        let face_he = cyl.loops[he.loop_].face;
        let face_twin = cyl.loops[twin.loop_].face;
        if face_he == face_twin {
            found_seam = true;
        }
    }
    assert!(
        found_seam,
        "Should have a seam edge with both halves on the lateral face"
    );
}

#[test]
fn cylinder_cap_loops_are_self_loops() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let mut single_he_loops = 0;
    for loop_ in cyl.loops.values() {
        let he = &cyl.half_edges[loop_.half_edge];
        if he.next == loop_.half_edge {
            single_he_loops += 1;
        }
    }
    assert_eq!(
        single_he_loops, 2,
        "Should have 2 single-half-edge cap loops"
    );
}

#[test]
fn cylinder_vertex_positions() {
    let cyl = BRep::cylinder(2.5, 3.0, Tolerance::default()).unwrap();
    let positions: Vec<_> = cyl.vertices.values().map(|v| v.position).collect();
    let has_bot = positions
        .iter()
        .any(|p| (p.x - 2.5).abs() < 1e-10 && p.y.abs() < 1e-10 && p.z.abs() < 1e-10);
    let has_top = positions
        .iter()
        .any(|p| (p.x - 2.5).abs() < 1e-10 && p.y.abs() < 1e-10 && (p.z - 3.0).abs() < 1e-10);
    assert!(has_bot, "Should have bottom vertex at (radius, 0, 0)");
    assert!(has_top, "Should have top vertex at (radius, 0, height)");
}

// --- Sphere tests ---

#[test]
fn sphere_topology() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    // Proper 1-face sphere: V=2, E=1, F=1
    assert_eq!(sph.vertices.len(), 2);
    assert_eq!(sph.edges.len(), 1);
    assert_eq!(sph.faces.len(), 1);
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
    assert_eq!(sph.faces.len(), 1, "Should have exactly 1 face");
    for (_, face) in &sph.faces {
        let surface = face
            .surface
            .as_ref()
            .expect("sphere face should have surface geometry");
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
fn sphere_seam_edge_both_half_edges_on_same_face() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for edge in sph.edges.values() {
        let he = &sph.half_edges[edge.half_edge];
        let twin = &sph.half_edges[he.twin];
        let face_he = sph.loops[he.loop_].face;
        let face_twin = sph.loops[twin.loop_].face;
        assert_eq!(
            face_he, face_twin,
            "Seam edge should have both halves on the same face"
        );
    }
}

#[test]
fn sphere_loop_has_two_half_edges() {
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    for loop_ in sph.loops.values() {
        let start = loop_.half_edge;
        let next = sph.half_edges[start].next;
        let after_next = sph.half_edges[next].next;
        assert_ne!(start, next, "Loop should have more than 1 half-edge");
        assert_eq!(after_next, start, "Loop should have exactly 2 half-edges");
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
fn cylinder_lateral_pcurves_are_degree_1() {
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // The lateral face has 4 half-edges with degree-1 pcurves.
    let mut deg1_count = 0;
    let mut deg2_count = 0;
    for (_, he) in &cyl.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        if pcurve.degree == 1 {
            deg1_count += 1;
        } else if pcurve.degree == 2 {
            deg2_count += 1;
        }
    }
    assert_eq!(deg1_count, 4, "Lateral face should have 4 degree-1 pcurves");
    assert_eq!(deg2_count, 2, "Cap faces should have 2 degree-2 pcurves");
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
fn cone_lateral_pcurves_are_degree_1() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    // The lateral face has 3 half-edges with degree-1 pcurves.
    // The cap face has 1 half-edge with degree-2 pcurve.
    let mut deg1_count = 0;
    let mut deg2_count = 0;
    for (_, he) in &cone.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        if pcurve.degree == 1 {
            deg1_count += 1;
        } else if pcurve.degree == 2 {
            deg2_count += 1;
        }
    }
    assert_eq!(deg1_count, 3, "Lateral face should have 3 degree-1 pcurves");
    assert_eq!(deg2_count, 1, "Cap face should have 1 degree-2 pcurve");
}

#[test]
fn cone_pcurves_have_valid_v_range() {
    // Cone v ranges from 0 (apex) to height (base).
    let height = 2.0;
    let cone = BRep::cone(1.0, height, Tolerance::default()).unwrap();
    // Check only the degree-1 (lateral) pcurves for v-range.
    for (_, he) in &cone.half_edges {
        let pcurve = he.pcurve.as_ref().unwrap();
        if pcurve.degree == 1 {
            for cp in &pcurve.control_points {
                assert!(
                    cp.y >= -1e-10 && cp.y <= height + 1e-10,
                    "Cone pcurve v should be in [0, height], got {}",
                    cp.y
                );
            }
        }
    }
}

#[test]
fn cone_has_surface_types() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    let mut cone_faces = 0;
    let mut plane_faces = 0;
    for face in cone.faces.values() {
        if let Some(ref sof) = face.surface {
            match sof.geometry {
                SurfaceGeometry::Cone { .. } => cone_faces += 1,
                SurfaceGeometry::Plane { .. } => plane_faces += 1,
                _ => {}
            }
        }
    }
    assert_eq!(cone_faces, 1, "Should have 1 conical face");
    assert_eq!(plane_faces, 1, "Should have 1 planar cap face");
}

// --- Torus tests ---

#[test]
fn torus_basic_succeeds() {
    let result = BRep::torus(2.0, 0.5, Tolerance::default());
    assert!(result.is_ok(), "Torus constructor should succeed");
}
