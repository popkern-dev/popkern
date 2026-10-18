use super::*;

fn z_cylinder(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

fn groove(center_z: f64, offset_x: f64) -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(offset_x, 0.0, center_z),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 0.7,
        minor_radius: 0.16,
    }
}

#[test]
fn a_groove_beside_its_wall_rolls_on_an_exact_tube_centre() {
    let (cyl, tor) = (z_cylinder(0.58), groove(0.2, 0.0));
    let (origin, axis) = shared_axis(&cyl, &tor).expect("coaxial");
    let s1 = section(&cyl, origin, axis).unwrap();
    let s2 = section(&tor, origin, axis).unwrap();
    let rho = 0.012;
    let mut all = Vec::new();
    for a in offsets(s1, rho) {
        for b in offsets(s2, rho) {
            all.extend(crossings(a, b));
        }
    }
    // The wall offset inward (r 0.568) crosses the groove offset into
    // the material (tube radius 0.172) below and above the ball ring.
    let dz = (0.172_f64 * 0.172 - 0.132 * 0.132).sqrt();
    for h in [0.2 - dz, 0.2 + dz] {
        assert!(
            all.iter()
                .any(|q| (q.0 - 0.568).abs() < 1e-14 && (q.1 - h).abs() < 1e-14),
            "{all:?}"
        );
    }
}

#[test]
fn supports_off_one_axis_are_not_coaxial() {
    assert!(shared_axis(&z_cylinder(0.58), &groove(0.2, 1e-6)).is_none());
    let tilted = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1e-4, 0.0, 1.0),
        radius: 0.58,
    };
    assert!(shared_axis(&tilted, &groove(0.2, 0.0)).is_none());
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.3),
        radius: 1.0,
        axis: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(shared_axis(&sphere, &groove(0.2, 0.0)).is_some());
}

#[test]
fn a_cone_section_is_its_generator() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.5,
    };
    let (origin, axis) = shared_axis(&cone, &z_cylinder(0.3)).unwrap();
    let Some(Section::Line { p, d }) = section(&cone, origin, axis) else {
        panic!("a cone's section is a line");
    };
    // The generator 2 from the apex (meridian origin), opening along
    // the cone's own axis (+h).
    let q = (p.0 + 2.0 * d.0, p.1 + 2.0 * d.1);
    assert_eq!(p, (0.0, 0.0));
    assert!(
        q.1 > 0.0 && (q.0 - q.1 * 0.5_f64.tan()).abs() < 1e-14,
        "{q:?}"
    );
}
