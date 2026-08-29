//! End-to-end Newton convergence tests for `solve_rolling_ball_contact`
//! against closed-form analytic solutions on plane-plane,
//! plane-cylinder, and plane-sphere (in both convex/outside and
//! concave/inside-cap configurations).

#[cfg(test)]
use crate::contact::{
    BallSide, ContactOptions, ContactRequest, ContactTrace, solve_rolling_ball_contact,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

fn unit_tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn newton_plane_plane_perpendicular_converges() {
    // Two perpendicular planes (x=0 and z=0) meeting along the y-axis.
    // For radius r, the rolling-ball center sits at (r, y, r).
    let plane_z0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let plane_x0 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let radius = 0.5;
    let request = ContactRequest {
        surface_0: &plane_z0,
        surface_1: &plane_x0,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius,
        // Far-from-truth seed exercises convergence.
        initial_uv0: (2.0, 0.0),
        uv0_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_uv1: (0.0, 2.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(0.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let mut trace = ContactTrace::default();
    let solution = solve_rolling_ball_contact(
        &request,
        &ContactOptions::default(),
        &unit_tol(),
        Some(&mut trace),
    )
    .expect("plane-plane should converge");
    let expected_center = Point3::new(radius, 0.0, radius);
    assert!(
        (solution.center - expected_center).length() < 1e-9,
        "center {:?} expected {:?}",
        solution.center,
        expected_center
    );
    assert!(
        solution.iterations <= 5,
        "expected ≤5 iter, got {}",
        solution.iterations
    );
    assert!(solution.residual < 1e-9);
    assert!(!trace.steps.is_empty());
    assert_eq!(trace.steps[0].iteration, 0);
}

#[test]
fn newton_plane_cylinder_converges() {
    // Cylinder along z-axis radius 2, plane z=0. They intersect in
    // x²+y²=4, z=0. A rolling ball of radius r at the convex corner
    // outside both surfaces has center ((R+r)cosθ, (R+r)sinθ, r).
    // At θ = 0: (2.3, 0, 0.3) for r = 0.3.
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.3;
    let request = ContactRequest {
        surface_0: &cyl,
        surface_1: &plane,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius: r,
        initial_uv0: (0.0, 0.0),
        uv0_range: ((-PI, 2.0 * PI), (-10.0, 10.0)),
        initial_uv1: (2.0, 0.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(2.0, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let solution =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None)
            .expect("plane-cylinder should converge");
    let expected = Point3::new(2.0 + r, 0.0, r);
    assert!(
        (solution.center - expected).length() < 1e-6,
        "plane-cylinder center {:?} expected {:?}",
        solution.center,
        expected
    );
    assert!(solution.iterations <= 10);
}

#[test]
fn newton_plane_sphere_converges_outside() {
    // Sphere R=2 at (0,0,1) intersects plane z=0 in x²+y²=3.
    // Convex (outside) ball of radius r=0.5 has center
    // (sqrt(6), 0, 0.5) on the spine plane y=0.
    //
    // popkern's sphere has natural normal pointing INWARD at the
    // equator (longitude×colatitude cross-product), so the convex
    // configuration uses BallSide::Inside on the sphere.
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.5;
    let spine_x = 3.0_f64.sqrt();
    let request = ContactRequest {
        surface_0: &sphere,
        surface_1: &plane,
        side_0: BallSide::Inside,
        side_1: BallSide::Outside,
        radius: r,
        initial_uv0: (0.0, PI / 2.0),
        uv0_range: ((-PI, 2.0 * PI), (0.0, PI)),
        initial_uv1: (spine_x, 0.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(spine_x, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let solution =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None)
            .expect("plane-sphere should converge");
    let expected = Point3::new(6.0_f64.sqrt(), 0.0, r);
    assert!(
        (solution.center - expected).length() < 1e-6,
        "plane-sphere center {:?} expected {:?}",
        solution.center,
        expected
    );
    assert!(solution.iterations <= 10);
}

#[test]
fn newton_plane_sphere_converges_inside_cap() {
    // Same sphere and plane as the outside test, but the ball sits
    // in the concave cap (inside the sphere AND above the plane).
    // Center at distance R-r = 1.5 from sphere center, z = +r:
    //   x² + (0.5 - 1)² = 2.25 → x = sqrt(2).
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let r = 0.5;
    let spine_x = 3.0_f64.sqrt();
    let request = ContactRequest {
        surface_0: &sphere,
        surface_1: &plane,
        side_0: BallSide::Outside, // +inward → ball inside sphere
        side_1: BallSide::Outside, // +z → ball above plane
        radius: r,
        initial_uv0: (0.0, PI / 2.0),
        uv0_range: ((-PI, 2.0 * PI), (0.0, PI)),
        initial_uv1: (spine_x, 0.0),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        spine_point: Point3::new(spine_x, 0.0, 0.0),
        spine_plane_normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let solution =
        solve_rolling_ball_contact(&request, &ContactOptions::default(), &unit_tol(), None)
            .expect("plane-sphere concave cap should converge");
    let expected = Point3::new(2.0_f64.sqrt(), 0.0, r);
    assert!(
        (solution.center - expected).length() < 1e-6,
        "plane-sphere concave center {:?} expected {:?}",
        solution.center,
        expected
    );
}
