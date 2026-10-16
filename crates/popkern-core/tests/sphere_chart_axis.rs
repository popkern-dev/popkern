//! The sphere chart's pole axis (`SurfaceGeometry::Sphere::axis`).
//!
//! The chart turns with its axis exactly as torus charts do, the world
//! +Z axis reproduces the chart spheres had before the field existed
//! bit for bit, and documents written without the field read as +Z.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

fn sphere(axis: Vec3) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.3, -1.2, 2.5),
        radius: 1.7,
        axis,
    }
}

#[test]
fn a_world_pole_evaluates_bit_identically_to_the_world_chart() {
    let s = sphere(Vec3::Z);
    for (u, v) in [(0.0, 0.0), (0.4, 0.9), (2.1, 1.3), (5.9, 3.0)] {
        let world = Point3::new(
            0.3 + 1.7 * v.psin() * u.pcos(),
            -1.2 + 1.7 * v.psin() * u.psin(),
            2.5 + 1.7 * v.pcos(),
        );
        let p = s.point_at(u, v);
        assert_eq!(
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()],
            [world.x.to_bits(), world.y.to_bits(), world.z.to_bits()],
            "({u}, {v})"
        );
    }
}

#[test]
fn a_tilted_pole_carries_the_chart_with_it() {
    let axis = Vec3::new(0.6, 0.0, 0.8);
    let s = sphere(axis);
    let center = Point3::new(0.3, -1.2, 2.5);
    // The pole sits on the axis, the equator is perpendicular to it,
    // and u = 0 starts at the first basis vector, as on a torus.
    let pole = s.point_at(1.3, 0.0);
    assert!(pole.distance_to(&(center + axis * 1.7)) < 1e-12);
    let (e1, _) = orthonormal_basis(axis);
    let origin = s.point_at(0.0, std::f64::consts::FRAC_PI_2);
    assert!(origin.distance_to(&(center + e1 * 1.7)) < 1e-12);
    // closest_point inverts point_at off the poles, and the analytic
    // partials agree with central differences.
    for (u, v) in [(0.4, 0.9), (2.1, 1.5), (5.9, 2.8)] {
        let p = s.point_at(u, v);
        let (cu, cv) = s.closest_point(&p).expect("projects");
        assert!(s.point_at(cu, cv).distance_to(&p) < 1e-12, "({u}, {v})");
        let h = 1e-6;
        let du = (s.point_at(u + h, v) - s.point_at(u - h, v)) * (0.5 / h);
        let dv = (s.point_at(u, v + h) - s.point_at(u, v - h)) * (0.5 / h);
        assert!((du - s.partial_u(u, v)).length() < 1e-8);
        assert!((dv - s.partial_v(u, v)).length() < 1e-8);
    }
}

#[test]
fn a_rotation_turns_the_sphere_axis() {
    let rotated = popkern_core::topology::BRep::sphere(1.0, Default::default())
        .expect("sphere")
        .transform(&Transform3::rotation_x(std::f64::consts::FRAC_PI_2))
        .expect("transform");
    let axes: Vec<Vec3> = rotated
        .faces
        .values()
        .filter_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { axis, .. }) => Some(*axis),
            _ => None,
        })
        .collect();
    assert!(!axes.is_empty());
    for axis in axes {
        assert!(
            (axis - Vec3::new(0.0, -1.0, 0.0)).length() < 1e-12,
            "{axis:?}"
        );
    }
    assert!(rotated.validate_all().is_empty());
}

#[test]
fn a_sphere_stored_without_an_axis_reads_as_world_z() {
    let json = r#"{"Sphere":{"center":{"x":0.0,"y":0.0,"z":0.0},"radius":2.0}}"#;
    let s: SurfaceGeometry = serde_json::from_str(json).expect("reads");
    assert_eq!(
        s,
        SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 2.0,
            axis: Vec3::Z,
        }
    );
}
