//! Tests for SketchPlane construction, flip, and from_face_offset error cases.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::SketchError;
use crate::plane::SketchPlane;

// ---------------------------------------------------------------------------
// from_face_offset with invalid face handle
// ---------------------------------------------------------------------------

#[test]
fn from_face_offset_invalid_face_returns_error() {
    let tol = Tolerance::default();
    let brep = BRep::new(tol);
    let bogus = FaceHandle::default();

    let result = SketchPlane::from_face_offset(&brep, bogus, 1.0);
    assert!(
        result.is_err(),
        "invalid face handle should return an error"
    );
    assert!(
        matches!(result, Err(SketchError::NonPlanarFace)),
        "expected NonPlanarFace error for invalid handle"
    );
}

// ---------------------------------------------------------------------------
// from_face_offset with zero offset
// ---------------------------------------------------------------------------

#[test]
fn from_face_offset_zero_keeps_origin() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit cube should build");

    // Find a planar face.
    let face = cube.faces.keys().next().expect("cube has faces");
    let base_plane = SketchPlane::from_face(&cube, face).expect("cube face should yield a plane");

    let offset_plane =
        SketchPlane::from_face_offset(&cube, face, 0.0).expect("zero offset should succeed");

    let d = offset_plane.origin - base_plane.origin;
    assert!(
        d.length() < 1e-12,
        "zero-offset origin should match base: delta = {d:?}"
    );
    // Normal, u, v should be identical.
    assert!((offset_plane.normal - base_plane.normal).length() < 1e-12);
    assert!((offset_plane.u_axis - base_plane.u_axis).length() < 1e-12);
    assert!((offset_plane.v_axis - base_plane.v_axis).length() < 1e-12);
}

// ---------------------------------------------------------------------------
// from_face_offset with positive offset
// ---------------------------------------------------------------------------

#[test]
fn from_face_offset_positive_shifts_along_normal() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit cube should build");
    let face = cube.faces.keys().next().expect("cube has faces");

    let base = SketchPlane::from_face(&cube, face).expect("plane ok");
    let dist = 3.5;
    let shifted = SketchPlane::from_face_offset(&cube, face, dist).expect("offset ok");

    let expected_origin = base.origin + base.normal * dist;
    let delta = shifted.origin - expected_origin;
    assert!(
        delta.length() < 1e-12,
        "positive offset origin mismatch: {delta:?}"
    );
}

// ---------------------------------------------------------------------------
// from_face_offset with negative offset
// ---------------------------------------------------------------------------

#[test]
fn from_face_offset_negative_shifts_opposite_normal() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit cube should build");
    let face = cube.faces.keys().next().expect("cube has faces");

    let base = SketchPlane::from_face(&cube, face).expect("plane ok");
    let dist = -2.0;
    let shifted = SketchPlane::from_face_offset(&cube, face, dist).expect("offset ok");

    let expected_origin = base.origin + base.normal * dist;
    let delta = shifted.origin - expected_origin;
    assert!(
        delta.length() < 1e-12,
        "negative offset origin mismatch: {delta:?}"
    );
    // Axes and normal should be unchanged.
    assert!((shifted.normal - base.normal).length() < 1e-12);
    assert!((shifted.u_axis - base.u_axis).length() < 1e-12);
    assert!((shifted.v_axis - base.v_axis).length() < 1e-12);
}

// ---------------------------------------------------------------------------
// from_face with non-planar (cylindrical) surface
// ---------------------------------------------------------------------------

#[test]
fn from_face_non_planar_returns_error() {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).expect("cylinder should build");

    // Find the cylindrical face (non-planar).
    let cyl_face = cyl
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .map(|(h, _)| h)
        .expect("cylinder BRep should have a cylindrical face");

    let result = SketchPlane::from_face(&cyl, cyl_face);
    assert!(
        matches!(result, Err(SketchError::NonPlanarFace)),
        "cylindrical face should yield NonPlanarFace error"
    );
}

// ---------------------------------------------------------------------------
// from_face_offset inherits the non-planar error
// ---------------------------------------------------------------------------

#[test]
fn from_face_offset_non_planar_returns_error() {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).expect("cylinder should build");

    let cyl_face = cyl
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .map(|(h, _)| h)
        .expect("cylinder should have a cylindrical face");

    let result = SketchPlane::from_face_offset(&cyl, cyl_face, 5.0);
    assert!(
        matches!(result, Err(SketchError::NonPlanarFace)),
        "from_face_offset on cylindrical face should return NonPlanarFace"
    );
}

// ---------------------------------------------------------------------------
// from_face with face that has no surface geometry
// ---------------------------------------------------------------------------

#[test]
fn from_face_no_surface_returns_error() {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);

    // Create minimal topology with a face that has no surface.
    let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs ok");
    // mvfs creates a face with no surface geometry.
    let face = brep.faces.keys().next().expect("should have a face");

    let result = SketchPlane::from_face(&brep, face);
    assert!(
        matches!(result, Err(SketchError::NonPlanarFace)),
        "face without surface should yield NonPlanarFace"
    );
}

// ---------------------------------------------------------------------------
// flip() reverses orientation
// ---------------------------------------------------------------------------

#[test]
fn flip_reverses_normal() {
    let mut plane = SketchPlane::xy();
    let original_normal = plane.normal;

    plane.flip();

    let sum = plane.normal + original_normal;
    assert!(
        sum.length() < 1e-12,
        "flipped normal should be negation of original: {sum:?}"
    );
}

#[test]
fn flip_swaps_axes() {
    let mut plane = SketchPlane::xy();
    let original_u = plane.u_axis;
    let original_v = plane.v_axis;

    plane.flip();

    assert!(
        (plane.u_axis - original_v).length() < 1e-12,
        "after flip, u_axis should be original v_axis"
    );
    assert!(
        (plane.v_axis - original_u).length() < 1e-12,
        "after flip, v_axis should be original u_axis"
    );
}

#[test]
fn flip_is_involution() {
    let mut plane = SketchPlane::xz();
    let original_origin = plane.origin;
    let original_u = plane.u_axis;
    let original_v = plane.v_axis;
    let original_normal = plane.normal;

    plane.flip();
    plane.flip();

    assert!((plane.origin - original_origin).length() < 1e-12);
    assert!((plane.u_axis - original_u).length() < 1e-12);
    assert!((plane.v_axis - original_v).length() < 1e-12);
    assert!((plane.normal - original_normal).length() < 1e-12);
}

#[test]
fn flip_preserves_right_hand_rule() {
    let mut plane = SketchPlane::xy();
    plane.flip();

    // After flip, u x v should equal the (new) normal.
    let cross = plane.u_axis.cross(&plane.v_axis);
    let diff = cross - plane.normal;
    assert!(
        diff.length() < 1e-12,
        "u x v should equal normal after flip: diff = {diff:?}"
    );
}

// ---------------------------------------------------------------------------
// from_face on a cube preserves orientation consistency
// ---------------------------------------------------------------------------

#[test]
fn from_face_cube_preserves_axes_orientation() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit cube should build");

    for (fh, _) in &cube.faces {
        let plane = SketchPlane::from_face(&cube, fh).expect("cube face should be planar");
        // u x v should equal normal.
        let cross = plane.u_axis.cross(&plane.v_axis);
        let diff = cross - plane.normal;
        assert!(
            diff.length() < 1e-6,
            "u x v should equal normal for face {fh:?}: diff = {diff:?}"
        );
    }
}

// ---------------------------------------------------------------------------
// new() with degenerate axes
// ---------------------------------------------------------------------------

#[test]
fn new_with_zero_u_axis_returns_error() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let zero = Vec3::new(0.0, 0.0, 0.0);
    let v = Vec3::new(0.0, 1.0, 0.0);
    let result = SketchPlane::new(origin, zero, v);
    assert!(matches!(result, Err(SketchError::DegenerateNormal)));
}

#[test]
fn new_with_parallel_axes_returns_error() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let u = Vec3::new(1.0, 0.0, 0.0);
    let v = Vec3::new(2.0, 0.0, 0.0); // parallel to u
    let result = SketchPlane::new(origin, u, v);
    assert!(matches!(result, Err(SketchError::DegenerateNormal)));
}

// ---------------------------------------------------------------------------
// to_3d / to_2d round-trip
// ---------------------------------------------------------------------------

#[test]
fn to_3d_and_to_2d_round_trip() {
    let plane = SketchPlane::xy();
    let x = 3.0;
    let y = -7.5;
    let p3d = plane.to_3d(x, y);
    let (rx, ry) = plane.to_2d(p3d);
    assert!((rx - x).abs() < 1e-12);
    assert!((ry - y).abs() < 1e-12);
}
