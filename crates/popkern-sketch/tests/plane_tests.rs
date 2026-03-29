//! Tests for the SketchPlane type.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::SketchPlane;

#[test]
fn to_3d_to_2d_round_trip_xy() {
    let plane = SketchPlane::xy();
    let p = plane.to_3d(3.0, 4.0);
    assert!((p.x - 3.0).abs() < 1e-10);
    assert!((p.y - 4.0).abs() < 1e-10);
    assert!((p.z - 0.0).abs() < 1e-10);

    let (x, y) = plane.to_2d(p);
    assert!((x - 3.0).abs() < 1e-10);
    assert!((y - 4.0).abs() < 1e-10);
}

#[test]
fn to_3d_to_2d_round_trip_xz() {
    let plane = SketchPlane::xz();
    let p = plane.to_3d(2.0, 5.0);
    // u_axis = (1,0,0), v_axis = (0,0,1)
    assert!((p.x - 2.0).abs() < 1e-10);
    assert!((p.y - 0.0).abs() < 1e-10);
    assert!((p.z - 5.0).abs() < 1e-10);

    let (x, y) = plane.to_2d(p);
    assert!((x - 2.0).abs() < 1e-10);
    assert!((y - 5.0).abs() < 1e-10);
}

#[test]
fn to_3d_to_2d_round_trip_yz() {
    let plane = SketchPlane::yz();
    let p = plane.to_3d(2.0, 5.0);
    // u_axis = (0,1,0), v_axis = (0,0,1)
    assert!((p.x - 0.0).abs() < 1e-10);
    assert!((p.y - 2.0).abs() < 1e-10);
    assert!((p.z - 5.0).abs() < 1e-10);

    let (x, y) = plane.to_2d(p);
    assert!((x - 2.0).abs() < 1e-10);
    assert!((y - 5.0).abs() < 1e-10);
}

#[test]
fn tilted_plane_from_normal_round_trip() {
    let normal = Vec3::new(1.0, 1.0, 1.0);
    let origin = Point3::new(10.0, 20.0, 30.0);
    let plane = SketchPlane::from_normal(origin, normal).unwrap();

    // Verify orthogonality
    assert!(plane.u_axis.dot(&plane.v_axis).abs() < 1e-10);
    assert!(plane.u_axis.dot(&plane.normal).abs() < 1e-10);
    assert!(plane.v_axis.dot(&plane.normal).abs() < 1e-10);

    // Round-trip
    let p3d = plane.to_3d(1.5, -2.3);
    let (x, y) = plane.to_2d(p3d);
    assert!((x - 1.5).abs() < 1e-10);
    assert!((y - (-2.3)).abs() < 1e-10);
}

#[test]
fn from_face_on_unit_cube() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();

    // Pick the first face and try to extract a plane
    let face_handle = cube.faces.keys().next().unwrap();
    let plane = SketchPlane::from_face(&cube, face_handle).unwrap();

    // Verify the plane has orthonormal axes
    assert!((plane.u_axis.length() - 1.0).abs() < 1e-10);
    assert!((plane.v_axis.length() - 1.0).abs() < 1e-10);
    assert!((plane.normal.length() - 1.0).abs() < 1e-10);
    assert!(plane.u_axis.dot(&plane.v_axis).abs() < 1e-10);
}

#[test]
fn from_face_offset_shifts_origin() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();

    let face_handle = cube.faces.keys().next().unwrap();
    let base_plane = SketchPlane::from_face(&cube, face_handle).unwrap();
    let offset_plane = SketchPlane::from_face_offset(&cube, face_handle, 5.0).unwrap();

    // Origin should be shifted by 5.0 along normal
    let expected = base_plane.origin + base_plane.normal * 5.0;
    assert!(offset_plane.origin.distance_to(&expected) < 1e-10);

    // Axes and normal should be the same
    assert!((offset_plane.u_axis.dot(&base_plane.u_axis) - 1.0).abs() < 1e-10);
    assert!((offset_plane.v_axis.dot(&base_plane.v_axis) - 1.0).abs() < 1e-10);
    assert!((offset_plane.normal.dot(&base_plane.normal) - 1.0).abs() < 1e-10);
}

#[test]
fn from_points_three_noncollinear() {
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let p2 = Point3::new(1.0, 2.0, 0.0);
    let p3 = Point3::new(1.0, 0.0, 3.0);
    let plane = SketchPlane::from_points(p1, p2, p3).unwrap();

    // Origin should be p1
    assert!(plane.origin.distance_to(&p1) < 1e-10);

    // u_axis should be along p1->p2 direction
    let expected_u = (p2 - p1).normalized().unwrap();
    assert!((plane.u_axis.dot(&expected_u) - 1.0).abs() < 1e-10);

    // Normal should be perpendicular to both u and v
    assert!(plane.u_axis.dot(&plane.normal).abs() < 1e-10);
    assert!(plane.v_axis.dot(&plane.normal).abs() < 1e-10);

    // Round-trip test
    let p3d = plane.to_3d(0.5, 0.3);
    let (x, y) = plane.to_2d(p3d);
    assert!((x - 0.5).abs() < 1e-10);
    assert!((y - 0.3).abs() < 1e-10);
}

#[test]
fn from_points_collinear_returns_error() {
    let p1 = Point3::new(0.0, 0.0, 0.0);
    let p2 = Point3::new(1.0, 0.0, 0.0);
    let p3 = Point3::new(2.0, 0.0, 0.0); // collinear!
    let result = SketchPlane::from_points(p1, p2, p3);
    assert!(result.is_err());
}

#[test]
fn from_points_coincident_returns_error() {
    let p = Point3::new(1.0, 2.0, 3.0);
    let result = SketchPlane::from_points(p, p, Point3::new(4.0, 5.0, 6.0));
    assert!(result.is_err());
}

#[test]
fn flip_reverses_normal() {
    let mut plane = SketchPlane::xy();
    let original_normal = plane.normal;
    let original_u = plane.u_axis;
    let original_v = plane.v_axis;

    plane.flip();

    // Normal should be negated
    assert!((plane.normal.x - (-original_normal.x)).abs() < 1e-10);
    assert!((plane.normal.y - (-original_normal.y)).abs() < 1e-10);
    assert!((plane.normal.z - (-original_normal.z)).abs() < 1e-10);

    // u and v should be swapped
    assert!((plane.u_axis.dot(&original_v) - 1.0).abs() < 1e-10);
    assert!((plane.v_axis.dot(&original_u) - 1.0).abs() < 1e-10);
}

#[test]
fn from_normal_with_z_aligned_normal() {
    // Normal along Z should still produce valid u/v axes
    let plane =
        SketchPlane::from_normal(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap();
    assert!(plane.u_axis.dot(&plane.v_axis).abs() < 1e-10);
    assert!(plane.u_axis.dot(&plane.normal).abs() < 1e-10);
    assert!((plane.normal.z - 1.0).abs() < 1e-10);
}

#[test]
fn from_normal_with_x_aligned_normal() {
    // Normal along X should pick appropriate candidate
    let plane =
        SketchPlane::from_normal(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)).unwrap();
    assert!(plane.u_axis.dot(&plane.v_axis).abs() < 1e-10);
    assert!(plane.u_axis.dot(&plane.normal).abs() < 1e-10);
    assert!((plane.normal.x - 1.0).abs() < 1e-10);
}

#[test]
fn from_normal_zero_vector_returns_error() {
    let result = SketchPlane::from_normal(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
    assert!(result.is_err());
}

#[test]
fn new_with_degenerate_axes_returns_error() {
    let result = SketchPlane::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0), // parallel to u -> cross product is zero
    );
    assert!(result.is_err());
}

#[test]
fn to_2d_projects_correctly_off_plane_point() {
    let plane = SketchPlane::xy();
    // A point above the XY plane should project to its XY coordinates
    let p = Point3::new(3.0, 4.0, 100.0);
    let (x, y) = plane.to_2d(p);
    assert!((x - 3.0).abs() < 1e-10);
    assert!((y - 4.0).abs() < 1e-10);
}
