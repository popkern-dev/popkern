//! Tests for side wall construction (split under the 500-line cap).

use super::side_walls_support::plane_from_quad;
use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn find_boundary_edges_one_removed_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let removed_set: HashSet<FaceHandle> = [first_face].into_iter().collect();
    let edges = find_boundary_edges(&brep, first_face, &removed_set).unwrap();
    // A cube face has 4 edges, all boundary (adjacent faces are kept).
    assert_eq!(edges.len(), 4);
}

#[test]
fn find_boundary_edges_two_adjacent_removed() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    // Remove first two faces. Some edges between them won't be boundary.
    let removed_set: HashSet<FaceHandle> = [face_handles[0], face_handles[1]].into_iter().collect();

    let edges_f0 = find_boundary_edges(&brep, face_handles[0], &removed_set).unwrap();
    let edges_f1 = find_boundary_edges(&brep, face_handles[1], &removed_set).unwrap();

    // Each face has 4 edges, but edges shared between the two removed faces
    // are not boundary edges. Depending on adjacency, we expect < 4 for at
    // least one of them (if adjacent) or 4 each (if non-adjacent).
    let total = edges_f0.len() + edges_f1.len();
    assert!((6..=8).contains(&total), "total boundary edges: {total}");
}

#[test]
fn plane_from_quad_basic() {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let p2 = Point3::new(1.0, 0.0, 1.0);
    let p3 = Point3::new(0.0, 0.0, 1.0);
    let surface = plane_from_quad(p0, p1, p2, p3).unwrap();
    match surface {
        SurfaceGeometry::Plane { normal, .. } => {
            // Normal should be along y-axis.
            assert!(normal.y.abs() > 0.9);
        }
        _ => panic!("expected Plane"),
    }
}

#[test]
fn plane_from_quad_collinear_returns_none() {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let p2 = Point3::new(2.0, 0.0, 0.0);
    let p3 = Point3::new(3.0, 0.0, 0.0);
    assert!(plane_from_quad(p0, p1, p2, p3).is_none());
}
