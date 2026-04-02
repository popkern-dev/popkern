use super::*;
use popkern_core::geometry::Vec3;

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn xz_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn test_tolerance() -> Tolerance {
    Tolerance {
        linear: 1e-3,
        angular: 1e-3,
        parametric: 1e-6,
    }
}

/// Plane at 45 degrees for variable-angle tests. Normal = (0, sin45, cos45).
fn plane_45() -> SurfaceGeometry {
    let c = std::f64::consts::FRAC_1_SQRT_2;
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, c, c),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, -c, c),
    }
}

#[test]
fn constant_center_curve_two_perpendicular_planes() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let range = (-5.0, 5.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);

    let result = constant_radius_center_curve(
        &s1,
        range,
        range,
        &s2,
        range,
        range,
        1.0,
        &midpoint,
        &test_tolerance(),
        true,
        true,
    );

    assert!(result.is_ok(), "center curve failed: {:?}", result.err());
    let cr = result.unwrap();
    assert!(!cr.footprint_1.is_empty(), "footprint_1 should have points");
    assert!(!cr.footprint_2.is_empty(), "footprint_2 should have points");

    let (t0, t1) = cr.center_curve.parameter_range();
    let mid = cr.center_curve.point_at((t0 + t1) * 0.5);
    assert!(
        mid.z.abs().max(0.5) > 0.5,
        "center z should be near +/-1.0, got {}",
        mid.z
    );
    assert!(
        mid.y.abs().max(0.5) > 0.5,
        "center y should be near +/-1.0, got {}",
        mid.y
    );
}

#[test]
fn variable_center_curve_two_planes() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let range = (-5.0, 5.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);

    let result = variable_radius_center_curve(
        &s1,
        range,
        range,
        &s2,
        range,
        range,
        0.5,
        1.5,
        &midpoint,
        &test_tolerance(),
        true,
        true,
    );

    assert!(
        result.is_ok(),
        "variable center curve failed: {:?}",
        result.err()
    );
    let cr = result.unwrap();
    assert!(!cr.footprint_1.is_empty(), "footprint_1 should have points");
    assert!(!cr.footprint_2.is_empty(), "footprint_2 should have points");
}

#[test]
fn variable_center_curve_endpoints_at_expected_distance() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let range = (-5.0, 5.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);

    let cr = variable_radius_center_curve(
        &s1,
        range,
        range,
        &s2,
        range,
        range,
        0.5,
        1.5,
        &midpoint,
        &test_tolerance(),
        true,
        true,
    )
    .expect("variable center curve should succeed");

    assert_eq!(cr.footprint_1.len(), cr.footprint_2.len());
    assert!(
        cr.footprint_1.len() >= 2,
        "need at least 2 footprint points"
    );

    for fp in &cr.footprint_1 {
        assert!(
            fp.z.abs() < 0.5,
            "footprint_1 should lie near the xy-plane (z~0), got z={}",
            fp.z,
        );
    }
    for fp in &cr.footprint_2 {
        assert!(
            fp.y.abs() < 0.5,
            "footprint_2 should lie near the xz-plane (y~0), got y={}",
            fp.y,
        );
    }
}

#[test]
fn variable_center_curve_equal_radii_degenerates() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let range = (-5.0, 5.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);
    let tol = test_tolerance();

    let const_result = constant_radius_center_curve(
        &s1, range, range, &s2, range, range, 1.0, &midpoint, &tol, true, true,
    );
    assert!(const_result.is_ok(), "constant path should succeed");

    let var_result = variable_radius_center_curve(
        &s1, range, range, &s2, range, range, 1.0, 1.0, &midpoint, &tol, true, true,
    );
    assert!(
        var_result.is_err(),
        "equal-radius variable path should degenerate for two perpendicular planes"
    );
}

#[test]
fn variable_center_curve_with_45_degree_planes() {
    let s1 = xy_plane();
    let s2 = plane_45();
    let range = (-5.0, 5.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);

    let result = variable_radius_center_curve(
        &s1,
        range,
        range,
        &s2,
        range,
        range,
        0.3,
        0.8,
        &midpoint,
        &test_tolerance(),
        true,
        true,
    );

    assert!(
        result.is_ok(),
        "variable center curve at 45deg failed: {:?}",
        result.err(),
    );
    let cr = result.unwrap();
    assert!(cr.footprint_1.len() >= 2);
    assert_eq!(cr.footprint_1.len(), cr.footprint_2.len());
}

#[test]
fn variable_center_curve_small_radius_range() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let range = (-5.0, 5.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);

    let result = variable_radius_center_curve(
        &s1,
        range,
        range,
        &s2,
        range,
        range,
        0.49,
        0.51,
        &midpoint,
        &test_tolerance(),
        true,
        true,
    );

    assert!(
        result.is_ok(),
        "small-range variable center curve failed: {:?}",
        result.err(),
    );
}
