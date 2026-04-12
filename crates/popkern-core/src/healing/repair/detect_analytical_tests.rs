use super::*;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

#[test]
fn detect_cone_from_analytical() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_6, // 30 degrees
    };
    // Sample a patch of the cone away from the apex.
    let n = 6;
    let u_range = (0.2, 2.5);
    let v_range = (1.0, 3.0);
    let mut points = Vec::with_capacity(n * n);
    let mut normals = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u_range.0 + (u_range.1 - u_range.0) * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let v = v_range.0 + (v_range.1 - v_range.0) * (j as f64) / ((n - 1) as f64);
            points.push(cone.point_at(u, v));
            normals.push(cone.normal_at(u, v));
        }
    }

    let result = try_detect_cone(&cone, &points, &normals, 1e-6);
    assert!(
        result.is_some(),
        "should detect cone from analytical surface, points={}, normals={}",
        points.len(),
        normals.len()
    );
    let (_, _, half_angle, dev) = result.unwrap();
    assert!(dev < 1e-6, "deviation {dev} should be within tolerance");
    assert!(
        (half_angle - std::f64::consts::FRAC_PI_6).abs() < 0.05,
        "half_angle should be ~pi/6, got {half_angle}"
    );
}

#[test]
fn detect_cone_from_nurbs() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(1.0, 2.0, 3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_4, // 45 degrees
    };
    let nurbs = SurfaceGeometry::Nurbs(cone.to_nurbs_in_range((0.2, 2.5), (1.0, 3.0)));
    let param_range = nurbs.parameter_range();
    let ((u0, u1), (v0, v1)) = param_range;
    let n = 6;
    let mut points = Vec::with_capacity(n * n);
    let mut normals = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u0 + (u1 - u0) * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let v = v0 + (v1 - v0) * (j as f64) / ((n - 1) as f64);
            points.push(nurbs.point_at(u, v));
            normals.push(nurbs.normal_at(u, v));
        }
    }

    let result = try_detect_cone(&nurbs, &points, &normals, 1.0);
    assert!(result.is_some(), "should detect cone from NURBS");
    let (_, _, _, dev) = result.unwrap();
    assert!(dev < 1.0, "deviation {dev} should be within tolerance");
}

#[test]
fn detect_torus_from_analytical() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    let n = 6;
    let u_range = (0.2, 2.5);
    let v_range = (0.2, 2.5);
    let mut points = Vec::with_capacity(n * n);
    let mut normals = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u_range.0 + (u_range.1 - u_range.0) * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let v = v_range.0 + (v_range.1 - v_range.0) * (j as f64) / ((n - 1) as f64);
            points.push(torus.point_at(u, v));
            normals.push(torus.normal_at(u, v));
        }
    }

    let result = try_detect_torus(&torus, &points, &normals, 0.5);
    assert!(
        result.is_some(),
        "should detect torus from analytical surface"
    );
    let (_, _, major_r, minor_r, dev) = result.unwrap();
    assert!(dev < 0.5, "deviation {dev} should be within tolerance");
    assert!(
        (major_r - 3.0).abs() < 0.5,
        "major_radius should be ~3.0, got {major_r}"
    );
    assert!(
        (minor_r - 1.0).abs() < 0.1,
        "minor_radius should be ~1.0, got {minor_r}"
    );
}

#[test]
fn detect_torus_from_nurbs() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(1.0, 2.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 4.0,
        minor_radius: 1.5,
    };
    let nurbs = SurfaceGeometry::Nurbs(torus.to_nurbs_in_range((0.2, 2.5), (0.2, 2.5)));
    let param_range = nurbs.parameter_range();
    let ((u0, u1), (v0, v1)) = param_range;
    let n = 6;
    let mut points = Vec::with_capacity(n * n);
    let mut normals = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u0 + (u1 - u0) * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let v = v0 + (v1 - v0) * (j as f64) / ((n - 1) as f64);
            points.push(nurbs.point_at(u, v));
            normals.push(nurbs.normal_at(u, v));
        }
    }

    let result = try_detect_torus(&nurbs, &points, &normals, 1.0);
    assert!(result.is_some(), "should detect torus from NURBS");
    let (_, _, _, _, dev) = result.unwrap();
    assert!(dev < 0.25, "deviation {dev} should be within tolerance");
}
