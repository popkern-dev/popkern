#[cfg(test)]
use super::*;
use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use std::f64::consts::TAU;

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

fn z5_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 5.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// The solver finds the intersection of the XY plane (z=0) and the XZ
/// plane (y=0), which is the X axis.  The advancing-front constraint
/// `(p - target) · X = 0` fixes x along the march, pinning the solution
/// to the unique point on the X axis with x=0.5, y=0, z=0.
#[test]
fn project_onto_two_perpendicular_planes() {
    let s1 = xy_plane(); // z = 0
    let s2 = xz_plane(); // y = 0

    // Use march_dir = X so the front constraint has non-zero Jacobian entries
    // for the XY plane (whose partials are in the XY plane).
    let march_dir = Vec3::new(1.0, 0.0, 0.0);
    let target = Point3::new(0.5, 0.0, 0.0);

    let result = march_project(
        &s1,
        0.1,
        0.1,
        &s2,
        0.1,
        0.4,
        &target,
        &march_dir,
        &Tolerance::default(),
    );

    assert!(result.is_some(), "expected convergence");
    let s = result.unwrap();
    assert!(s.point.y.abs() < 1e-8, "point should lie on XZ plane (y=0)");
    assert!(s.point.z.abs() < 1e-8, "point should lie on XY plane (z=0)");
}

/// Two parallel planes (z=0 and z=5) never intersect; the Jacobian for
/// the surface-difference rows is singular.
#[test]
fn fails_for_parallel_planes() {
    let s1 = xy_plane(); // z = 0
    let s2 = z5_plane(); // z = 5

    let result = march_project(
        &s1,
        0.0,
        0.0,
        &s2,
        0.0,
        0.0,
        &Point3::new(0.0, 0.0, 2.5),
        &Vec3::new(1.0, 0.0, 0.0),
        &Tolerance::default(),
    );

    assert!(result.is_none(), "parallel planes should not converge");
}

/// Two perpendicular planes (z=0 normal and x=0 normal) intersect along the
/// y-axis.  Marching from the origin should produce samples along y with
/// x ~ 0 and z ~ 0.
#[test]
fn march_two_perpendicular_planes() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let (samples, is_closed) = march(
        &s1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &s2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &start,
        &Tolerance::default(),
        &DefaultDegeneracyHandler,
    );
    assert!(samples.len() > 2, "Should produce multiple samples");
    assert!(!is_closed, "Planes produce open intersection");
    // All samples should have x ~ 0 and z ~ 0 (intersection is y-axis)
    for s in &samples {
        assert!(
            (s.point.x).abs() < 1e-4,
            "x should be ~0, got {}",
            s.point.x
        );
        assert!(
            (s.point.z).abs() < 1e-4,
            "z should be ~0, got {}",
            s.point.z
        );
    }
}

/// `wrap_parameter` should fold values outside the range back in for
/// periodic parameters, and leave non-periodic values unchanged.
#[test]
fn wrap_parameter_basic() {
    // Just past the top of [0, TAU]: wraps back near the bottom.
    let wrapped_high = wrap_parameter(6.5, (0.0, TAU), true);
    assert!(
        (wrapped_high - (6.5 - TAU)).abs() < 1e-10,
        "expected {}, got {}",
        6.5 - TAU,
        wrapped_high
    );

    // Just below 0 on [0, TAU]: wraps to near TAU.
    let wrapped_low = wrap_parameter(-0.5, (0.0, TAU), true);
    assert!(
        (wrapped_low - (TAU - 0.5)).abs() < 1e-10,
        "expected {}, got {}",
        TAU - 0.5,
        wrapped_low
    );

    // Non-periodic: value is returned unchanged.
    assert_eq!(
        wrap_parameter(6.5, (0.0, TAU), false),
        6.5,
        "non-periodic value should be unchanged"
    );
}

/// The intersection of a unit cylinder (u periodic, [0, TAU]) and a
/// horizontal plane at z=0.5 is a closed circle.  Marching should wrap
/// around the seam at u=0/2π and detect the closed curve.
///
/// We supply the starting point directly to avoid the expensive
/// `find_starting_points` subdivision and keep the test fast.
#[test]
fn march_crosses_cylinder_seam() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // At u1=0, v1=0.5 the cylinder gives point (1, 0, 0.5).
    // On the plane at z=0.5, (1, 0, 0.5) maps to u2=1, v2=0.
    let start = MarchSample {
        u1: 0.0,
        v1: 0.5,
        u2: 1.0,
        v2: 0.0,
        point: Point3::new(1.0, 0.0, 0.5),
    };
    let tol = Tolerance::default().with_linear(1e-3);
    let (samples, is_closed) = march(
        &cylinder,
        (0.0, TAU),
        (0.0, 1.0),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    assert!(
        samples.len() > 3,
        "should produce multiple samples, got {}",
        samples.len()
    );
    assert!(
        is_closed,
        "cylinder-plane intersection should be a closed curve"
    );
    // All samples should be at z ~ 0.5 (on the plane) and radius ~ 1.0 (on the cylinder).
    for s in &samples {
        assert!(
            (s.point.z - 0.5).abs() < 1e-2,
            "z should be ~0.5, got {}",
            s.point.z
        );
        let r = (s.point.x * s.point.x + s.point.y * s.point.y).sqrt();
        assert!((r - 1.0).abs() < 1e-2, "radius should be ~1.0, got {}", r);
    }
}
