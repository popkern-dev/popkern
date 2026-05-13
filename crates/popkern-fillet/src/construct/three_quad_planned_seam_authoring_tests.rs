use popkern_core::geometry::Vec3;

use super::*;

#[test]
fn sphere_plane_intersection_triple_validates_for_non_tangent_plane() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let seam = build_sphere_plane_intersection_triple(
        Point3::new(0.0, 0.0, 0.0),
        1.0,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        &sphere,
        &plane,
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    )
    .expect("equatorial sphere-plane circle should author");

    validate_sampled_agreement(
        &seam.curve,
        &seam.cap_pcurve,
        &seam.other_pcurve,
        &sphere,
        &plane,
    )
    .expect("authored seam should pass sampled agreement");
}

#[test]
fn sphere_plane_intersection_triple_declines_tangent_plane() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let result = build_sphere_plane_intersection_triple(
        Point3::new(0.0, 0.0, 0.0),
        1.0,
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, 1.0),
        &sphere,
        &plane,
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
    );
    let err = match result {
        Ok(_) => panic!("tangent sphere-plane contact should not author a boundary curve"),
        Err(err) => err,
    };
    assert!(err.contains("empty or tangent"));
}
