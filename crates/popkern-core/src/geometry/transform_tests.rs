use super::*;
use crate::geometry::Point3;
use std::f64::consts::FRAC_PI_2;

#[test]
fn identity_transform_point_unchanged() {
    let t = Transform3::identity();
    let p = Point3::new(1.0, 2.0, 3.0);
    let q = t.apply_point(&p);
    assert!((q.x - 1.0).abs() < 1e-12);
    assert!((q.y - 2.0).abs() < 1e-12);
    assert!((q.z - 3.0).abs() < 1e-12);
}

#[test]
fn translation_moves_point() {
    let t = Transform3::translation(10.0, 20.0, 30.0);
    let p = Point3::new(1.0, 2.0, 3.0);
    let q = t.apply_point(&p);
    assert!((q.x - 11.0).abs() < 1e-12);
    assert!((q.y - 22.0).abs() < 1e-12);
    assert!((q.z - 33.0).abs() < 1e-12);
}

#[test]
fn translation_does_not_affect_direction() {
    let t = Transform3::translation(10.0, 20.0, 30.0);
    let v = Vec3::new(1.0, 0.0, 0.0);
    let u = t.apply_direction(&v);
    assert!((u.x - 1.0).abs() < 1e-12);
    assert!((u.y).abs() < 1e-12);
    assert!((u.z).abs() < 1e-12);
}

#[test]
fn rotation_z_90_degrees() {
    let t = Transform3::rotation_z(FRAC_PI_2);
    let p = Point3::new(1.0, 0.0, 0.0);
    let q = t.apply_point(&p);
    assert!((q.x).abs() < 1e-12);
    assert!((q.y - 1.0).abs() < 1e-12);
    assert!((q.z).abs() < 1e-12);
}

#[test]
fn rotation_x_90_degrees() {
    let t = Transform3::rotation_x(FRAC_PI_2);
    let p = Point3::new(0.0, 1.0, 0.0);
    let q = t.apply_point(&p);
    assert!((q.x).abs() < 1e-12);
    assert!((q.y).abs() < 1e-12);
    assert!((q.z - 1.0).abs() < 1e-12);
}

#[test]
fn rotation_y_90_degrees() {
    let t = Transform3::rotation_y(FRAC_PI_2);
    let p = Point3::new(0.0, 0.0, 1.0);
    let q = t.apply_point(&p);
    assert!((q.x - 1.0).abs() < 1e-12);
    assert!((q.y).abs() < 1e-12);
    assert!((q.z).abs() < 1e-12);
}

#[test]
fn uniform_scale() {
    let t = Transform3::scale_uniform(2.0);
    let p = Point3::new(1.0, 2.0, 3.0);
    let q = t.apply_point(&p);
    assert!((q.x - 2.0).abs() < 1e-12);
    assert!((q.y - 4.0).abs() < 1e-12);
    assert!((q.z - 6.0).abs() < 1e-12);
}

#[test]
fn non_uniform_scale() {
    let t = Transform3::scale(2.0, 3.0, 4.0);
    let p = Point3::new(1.0, 1.0, 1.0);
    let q = t.apply_point(&p);
    assert!((q.x - 2.0).abs() < 1e-12);
    assert!((q.y - 3.0).abs() < 1e-12);
    assert!((q.z - 4.0).abs() < 1e-12);
}

#[test]
fn compose_translation_then_rotation() {
    let t = Transform3::translation(1.0, 0.0, 0.0);
    let r = Transform3::rotation_z(FRAC_PI_2);
    let combined = r.compose(&t);
    let p = Point3::new(1.0, 0.0, 0.0);
    let q = combined.apply_point(&p);
    assert!((q.x).abs() < 1e-12);
    assert!((q.y - 2.0).abs() < 1e-12);
}

#[test]
fn inverse_undoes_transform() {
    let t = Transform3::translation(3.0, 4.0, 5.0);
    let inv = t.inverse().expect("translation is invertible");
    let combined = inv.compose(&t);
    let p = Point3::new(7.0, 8.0, 9.0);
    let q = combined.apply_point(&p);
    assert!((q.x - 7.0).abs() < 1e-12);
    assert!((q.y - 8.0).abs() < 1e-12);
    assert!((q.z - 9.0).abs() < 1e-12);
}

#[test]
fn inverse_of_rotation() {
    let t = Transform3::rotation_z(1.23);
    let inv = t.inverse().expect("rotation is invertible");
    let p = Point3::new(1.0, 2.0, 3.0);
    let q = inv.apply_point(&t.apply_point(&p));
    assert!((q.x - 1.0).abs() < 1e-10);
    assert!((q.y - 2.0).abs() < 1e-10);
    assert!((q.z - 3.0).abs() < 1e-10);
}

#[test]
fn singular_matrix_has_no_inverse() {
    let mut t = Transform3::identity();
    t.matrix[0] = [0.0; 4];
    assert!(t.inverse().is_err());
}

#[test]
fn from_axis_placement_identity_axes() {
    let t = Transform3::from_axis_placement(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
    )
    .expect("valid axes");
    let p = Point3::new(1.0, 2.0, 3.0);
    let q = t.apply_point(&p);
    assert!((q.x - 1.0).abs() < 1e-12);
    assert!((q.y - 2.0).abs() < 1e-12);
    assert!((q.z - 3.0).abs() < 1e-12);
}

#[test]
fn from_axis_placement_with_offset() {
    let t = Transform3::from_axis_placement(
        Point3::new(10.0, 20.0, 30.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
    )
    .expect("valid axes");
    let p = Point3::new(0.0, 0.0, 0.0);
    let q = t.apply_point(&p);
    assert!((q.x - 10.0).abs() < 1e-12);
    assert!((q.y - 20.0).abs() < 1e-12);
    assert!((q.z - 30.0).abs() < 1e-12);
}

#[test]
fn approx_eq_identical() {
    let t = Transform3::rotation_z(1.0);
    assert!(t.approx_eq(&t, 1e-12));
}

#[test]
fn approx_eq_different() {
    let a = Transform3::translation(1.0, 0.0, 0.0);
    let b = Transform3::translation(2.0, 0.0, 0.0);
    assert!(!a.approx_eq(&b, 1e-6));
}

#[test]
fn rotation_around_z_axis_at_origin() {
    // 90 degrees around Z at origin: (1,0,0) -> (0,1,0)
    let t = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::FRAC_PI_2,
    )
    .unwrap();
    let p = t.apply_point(&Point3::new(1.0, 0.0, 0.0));
    assert!((p.x).abs() < 1e-10);
    assert!((p.y - 1.0).abs() < 1e-10);
    assert!((p.z).abs() < 1e-10);
}

#[test]
fn rotation_around_offset_axis() {
    // 180 degrees around Z axis at (1,0,0): (2,0,0) -> (0,0,0)
    let t = Transform3::rotation_around_axis(
        Point3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::PI,
    )
    .unwrap();
    let p = t.apply_point(&Point3::new(2.0, 0.0, 0.0));
    assert!((p.x).abs() < 1e-10);
    assert!((p.y).abs() < 1e-10);
}

#[test]
fn rotation_origin_is_fixed() {
    let origin = Point3::new(3.0, 4.0, 5.0);
    let t = Transform3::rotation_around_axis(origin, Vec3::new(1.0, 1.0, 1.0), 1.23).unwrap();
    let p = t.apply_point(&origin);
    assert!((p.x - origin.x).abs() < 1e-10);
    assert!((p.y - origin.y).abs() < 1e-10);
    assert!((p.z - origin.z).abs() < 1e-10);
}

#[test]
fn reflection_across_yz_plane() {
    // Mirror across YZ at origin: (3,4,5) -> (-3,4,5)
    let t = Transform3::reflection(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let p = t.apply_point(&Point3::new(3.0, 4.0, 5.0));
    assert!((p.x - (-3.0)).abs() < 1e-10);
    assert!((p.y - 4.0).abs() < 1e-10);
    assert!((p.z - 5.0).abs() < 1e-10);
}

#[test]
fn reflection_origin_is_fixed() {
    let origin = Point3::new(1.0, 2.0, 3.0);
    let t = Transform3::reflection(origin, Vec3::new(0.0, 1.0, 0.0)).unwrap();
    let p = t.apply_point(&origin);
    assert!((p.x - origin.x).abs() < 1e-10);
    assert!((p.y - origin.y).abs() < 1e-10);
    assert!((p.z - origin.z).abs() < 1e-10);
}

#[test]
fn reflection_is_involution() {
    // Reflecting twice returns to original
    let t = Transform3::reflection(Point3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let t2 = t.compose(&t);
    let p = Point3::new(7.0, 8.0, 9.0);
    let p2 = t2.apply_point(&p);
    assert!((p2.x - p.x).abs() < 1e-10);
    assert!((p2.y - p.y).abs() < 1e-10);
    assert!((p2.z - p.z).abs() < 1e-10);
}

#[test]
fn rotation_zero_axis_errors() {
    let r =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0), 1.0);
    assert!(r.is_err());
}

#[test]
fn reflection_zero_normal_errors() {
    let r = Transform3::reflection(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
    assert!(r.is_err());
}
