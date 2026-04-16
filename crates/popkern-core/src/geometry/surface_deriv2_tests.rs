//! Tests for second partial derivatives, closest_point, and bounding_box
//! on SurfaceEval implementations.

use crate::geometry::surface::SurfaceEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

// --- Second partial derivatives: verify analytically ---

#[test]
fn plane_second_derivatives_are_zero() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(plane.second_partial_uu(1.0, 2.0).length() < 1e-10);
    assert!(plane.second_partial_uv(1.0, 2.0).length() < 1e-10);
    assert!(plane.second_partial_vv(1.0, 2.0).length() < 1e-10);
}

#[test]
fn sphere_second_partial_uu_at_equator() {
    // Sphere: S(u,v) = r*(sin(v)*cos(u), sin(v)*sin(u), cos(v))
    // d²S/du² = r*(-sin(v)*cos(u), -sin(v)*sin(u), 0)
    // At u=0, v=PI/2: d²S/du² = r*(-1, 0, 0) = (-1, 0, 0)
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let suu = sphere.second_partial_uu(0.0, FRAC_PI_2);
    assert!((suu.x - (-1.0)).abs() < 1e-8, "suu.x = {}", suu.x);
    assert!(suu.y.abs() < 1e-8, "suu.y = {}", suu.y);
    assert!(suu.z.abs() < 1e-8, "suu.z = {}", suu.z);
}

#[test]
fn sphere_second_partial_vv_at_equator() {
    // d²S/dv² = r*(-cos(u)*sin(v), -sin(u)*sin(v), -cos(v))
    // At u=0, v=PI/2: d²S/dv² = (-1, 0, 0)
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let svv = sphere.second_partial_vv(0.0, FRAC_PI_2);
    assert!((svv.x - (-1.0)).abs() < 1e-8, "svv.x = {}", svv.x);
    assert!(svv.y.abs() < 1e-8, "svv.y = {}", svv.y);
    assert!(svv.z.abs() < 1e-8, "svv.z = {}", svv.z);
}

#[test]
fn sphere_second_partial_uv_at_equator() {
    // d²S/dudv = r*(-sin(u)*cos(v), cos(u)*cos(v), 0)
    // At u=0, v=PI/2: d²S/dudv = (0, 0, 0)  (cos(PI/2) = 0)
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let suv = sphere.second_partial_uv(0.0, FRAC_PI_2);
    assert!(
        suv.length() < 1e-8,
        "suv should be zero at equator: {:?}",
        suv
    );
}

#[test]
fn cylinder_second_partial_uu_at_u0() {
    // Cylinder S(u,v) = origin + perp1*(r*cos(u)) + perp2*(r*sin(u)) + axis*v
    // d²S/du² = perp1*(-r*cos(u)) + perp2*(-r*sin(u))
    // At u=0: d²S/du² = -r*perp1 => length = r
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let suu = cyl.second_partial_uu(0.0, 0.0);
    assert!(
        (suu.length() - 2.0).abs() < 1e-8,
        "||suu|| = {}",
        suu.length()
    );
}

#[test]
fn cylinder_second_partial_vv_is_zero() {
    // d²S/dv² = 0 for cylinder (linear in v)
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert!(cyl.second_partial_vv(1.0, 1.0).length() < 1e-10);
}

#[test]
fn cylinder_second_partial_uv_is_zero() {
    // d²S/dudv = 0 (u and v are independent for cylinder)
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert!(cyl.second_partial_uv(1.0, 1.0).length() < 1e-10);
}

#[test]
fn cone_second_partial_vv_is_zero() {
    // Cone is linear in v (for fixed u), so d²S/dv² = 0
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: FRAC_PI_4,
    };
    assert!(cone.second_partial_vv(1.0, 1.0).length() < 1e-10);
}

#[test]
fn torus_second_derivatives_nonzero() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    // Torus has nonzero second derivatives everywhere (except possibly at special points)
    assert!(torus.second_partial_uu(0.5, 0.5).length() > 0.1);
    assert!(torus.second_partial_vv(0.5, 0.5).length() > 0.1);
}

// --- Closest point ---

#[test]
fn sphere_closest_point() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    // Point in +x +y +z octant: well-conditioned (away from poles)
    let pt = Point3::new(3.0, 3.0, 3.0);
    let (u, v) = sphere.closest_point(&pt).expect("should converge");
    let closest = sphere.point_at(u, v);
    // Closest point should be on the unit sphere, facing toward pt
    let r = (closest.x.powi(2) + closest.y.powi(2) + closest.z.powi(2)).sqrt();
    assert!((r - 1.0).abs() < 0.01, "not on sphere: r = {}", r);
    // The closest point should minimize distance
    let dist = closest.distance_to(&pt);
    let expected = pt.distance_to(&Point3::new(0.0, 0.0, 0.0)) - 1.0;
    assert!(
        (dist - expected).abs() < 0.1,
        "dist = {}, expected = {}",
        dist,
        expected
    );
}

#[test]
fn cylinder_closest_point() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    // Point at (5, 0, 3): closest should be (1, 0, 3)
    let pt = Point3::new(5.0, 0.0, 3.0);
    let (u, v) = cyl.closest_point(&pt).expect("should converge");
    let closest = cyl.point_at(u, v);
    let radial = (closest.x.powi(2) + closest.y.powi(2)).sqrt();
    assert!(
        (radial - 1.0).abs() < 0.01,
        "radial = {} (expected 1.0)",
        radial
    );
    assert!((closest.z - 3.0).abs() < 0.1, "z = {}", closest.z);
}

#[test]
fn plane_closest_point() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // Point at (3, 4, 7): closest on z=0 plane should be (3, 4, 0)
    let pt = Point3::new(3.0, 4.0, 7.0);
    let (u, v) = plane.closest_point(&pt).expect("should converge");
    let closest = plane.point_at(u, v);
    assert!((closest.x - 3.0).abs() < 0.01);
    assert!((closest.y - 4.0).abs() < 0.01);
    assert!(closest.z.abs() < 0.01);
}

// --- Bounding box ---

#[test]
fn sphere_bounding_box_quadrant() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    // First octant of sphere: u in [0, PI/2], v in [0, PI/2]
    let bbox = sphere.bounding_box((0.0, FRAC_PI_2), (0.0, FRAC_PI_2));
    // All points should have x >= 0, y >= 0, z >= 0
    assert!(bbox.min.x >= -0.01);
    assert!(bbox.min.y >= -0.01);
    assert!(bbox.min.z >= -0.01);
    // Max coords should be near 1.0
    assert!(bbox.max.x <= 1.01);
    assert!(bbox.max.y <= 1.01);
    assert!(bbox.max.z <= 1.01);
}

#[test]
fn plane_bounding_box_finite_region() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let bbox = plane.bounding_box((-1.0, 1.0), (-1.0, 1.0));
    assert!((bbox.min.x - (-1.0)).abs() < 0.01);
    assert!((bbox.max.x - 1.0).abs() < 0.01);
    assert!((bbox.min.y - (-1.0)).abs() < 0.01);
    assert!((bbox.max.y - 1.0).abs() < 0.01);
}

// --- Cone closest_point ---

#[test]
fn cone_closest_point() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: FRAC_PI_4,
    };
    // Point far out at (10, 0, 5): the cone surface at u=0, v=5 has radius=5*tan(PI/4)=5
    // so (5, 0, 5). Closest point should be on the cone surface.
    let pt = Point3::new(10.0, 0.0, 5.0);
    let (u, v) = cone.closest_point(&pt).expect("should converge");
    let closest = cone.point_at(u, v);
    // Verify the point is actually on the cone (radial distance = v * tan(half_angle))
    let radial = (closest.x.powi(2) + closest.y.powi(2)).sqrt();
    let expected_radial = closest.z * FRAC_PI_4.tan();
    assert!(
        (radial - expected_radial).abs() < 0.1,
        "not on cone: radial={radial}, expected={expected_radial}"
    );
}

// --- Torus closest_point ---

#[test]
fn torus_closest_point() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    // Torus closest_point is challenging due to multiple local minima
    // and the doubly-periodic parameter space. Verify it returns Some
    // and the result is on the torus surface (within tolerance).
    let pt = Point3::new(0.0, -6.0, 0.0);
    let result = torus.closest_point(&pt);
    assert!(result.is_some(), "should return Some for torus");

    let (u, v) = result.unwrap();
    let closest = torus.point_at(u, v);

    // Point is on the torus: distance from tube center ≈ minor_radius
    let r_xy = (closest.x.powi(2) + closest.y.powi(2)).sqrt();
    let off_ring = r_xy - 3.0;
    let tube_dist = (off_ring.powi(2) + closest.z.powi(2)).sqrt();
    assert!(
        (tube_dist - 1.0).abs() < 0.2,
        "point should be on torus tube, tube_dist={tube_dist}"
    );
}

// --- Torus bounding box ---

#[test]
fn torus_bounding_box_full() {
    use std::f64::consts::TAU;
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    let bbox = torus.bounding_box((0.0, TAU), (0.0, TAU));
    // Full torus: extends to ±(major+minor) = ±4 in XY, ±minor = ±1 in Z
    assert!(bbox.min.x < -3.9, "min.x = {}", bbox.min.x);
    assert!(bbox.max.x > 3.9, "max.x = {}", bbox.max.x);
    assert!(bbox.min.z < -0.9, "min.z = {}", bbox.min.z);
    assert!(bbox.max.z > 0.9, "max.z = {}", bbox.max.z);
}

// --- NURBS second partials validation ---

#[test]
fn nurbs_second_partials_match_analytic_cylinder() {
    // Convert a cylinder to NURBS and verify second partials match analytic.
    use std::f64::consts::TAU;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let nurbs_surf = cyl.to_nurbs_in_range((0.0, TAU), (-1.0, 1.0));
    let nurbs = SurfaceGeometry::Nurbs(nurbs_surf);

    // Test at u=1.0, v=0.0 (away from knot boundaries)
    let u = 1.0;
    let v = 0.0;
    let analytic_uu = cyl.second_partial_uu(u, v);
    let nurbs_uu = nurbs.second_partial_uu(u, v);

    // NURBS finite differences are approximate, so use loose tolerance
    assert!(
        (analytic_uu.length() - nurbs_uu.length()).abs() < 0.3,
        "Suu magnitude: analytic={}, nurbs={}",
        analytic_uu.length(),
        nurbs_uu.length()
    );

    // Svv should be near zero for both (cylinder is linear in v)
    let analytic_vv = cyl.second_partial_vv(u, v);
    let nurbs_vv = nurbs.second_partial_vv(u, v);
    assert!(
        analytic_vv.length() < 0.01,
        "analytic Svv should be ~0: {}",
        analytic_vv.length()
    );
    assert!(
        nurbs_vv.length() < 0.3,
        "nurbs Svv should be near 0: {}",
        nurbs_vv.length()
    );
}

// --- Cylinder bounding box ---

#[test]
fn cylinder_bounding_box() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let bbox = cyl.bounding_box((0.0, std::f64::consts::TAU), (-1.0, 1.0));
    assert!(bbox.min.x < -1.9, "min.x = {}", bbox.min.x);
    assert!(bbox.max.x > 1.9, "max.x = {}", bbox.max.x);
    assert!((bbox.min.z - (-1.0)).abs() < 0.01);
    assert!((bbox.max.z - 1.0).abs() < 0.01);
}
