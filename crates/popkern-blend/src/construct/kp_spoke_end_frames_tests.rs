//! Spoke end frames: plane normals and leg projection.

use super::*;

#[test]
fn a_corner_plane_is_spanned_by_the_hex_edges_and_declines_when_smooth() {
    let n = corner_plane_normal(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.3, 2.0, 0.0))
        .unwrap_or_else(|| panic!("a real corner has a plane"));
    assert!((n.z.abs() - 1.0).abs() < 1e-12);
    // The rim runs smoothly through the corner: no plane is imposed.
    assert!(corner_plane_normal(Vec3::new(1.0, 0.0, 0.0), Vec3::new(-2.0, 1e-5, 0.0)).is_none());
}

#[test]
fn the_center_plane_of_coplanar_legs_is_their_plane() {
    let legs: Vec<Vec3> = (0..3)
        .map(|k| {
            let a = std::f64::consts::TAU * k as f64 / 3.0 + 0.2 * k as f64;
            Vec3::new(a.cos(), a.sin(), 0.0) * (1.0 + k as f64)
        })
        .collect();
    let n = center_plane_normal(&legs).unwrap_or_else(|| panic!("plane"));
    assert!((n.z - 1.0).abs() < 1e-12, "{n:?}");
    // A tilted leg is outvoted; the normal stays near the majority plane.
    let mut tilted = legs.clone();
    tilted[1] = tilted[1] + Vec3::new(0.0, 0.0, 0.3);
    let n = center_plane_normal(&tilted).unwrap_or_else(|| panic!("plane"));
    assert!(n.z > 0.9, "{n:?}");
}

#[test]
fn a_leg_is_rotated_into_the_plane_keeping_its_length() {
    let anchor = Point3::new(1.0, 2.0, 3.0);
    let handle = Point3::new(1.3, 2.4, 3.5);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let moved = leg_into_plane(anchor, handle, normal);
    let (before, after) = (handle - anchor, moved - anchor);
    assert!((before.length() - after.length()).abs() < 1e-12);
    assert!(after.dot(&normal).abs() < 1e-12);
    // The in-plane direction is kept.
    assert!(after.cross(&Vec3::new(0.3, 0.4, 0.0)).length() < 1e-12);
    // A leg already in the plane is unchanged.
    let flat = Point3::new(1.3, 2.4, 3.0);
    let same = leg_into_plane(anchor, flat, normal);
    assert!((same - flat).length() < 1e-15);
}
