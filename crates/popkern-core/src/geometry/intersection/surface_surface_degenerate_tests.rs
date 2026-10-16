use super::*;
use crate::geometry::{Point3, Vec3};
use std::f64::consts::{PI, TAU};

/// A coarser tolerance so marching takes fewer, larger steps in tests.
fn test_tolerance() -> Tolerance {
    Tolerance::default()
        .with_linear(1e-3)
        .with_angular(1e-3)
        .with_parametric(1e-6)
}

// ── Parallel planes (no intersection) ──────────────────────────────────

#[test]
fn parallel_planes_separated_no_intersection() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 5.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let result = intersect_surfaces(
        &s1,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &s2,
        (-10.0, 10.0),
        (-10.0, 10.0),
        &test_tolerance(),
    )
    .unwrap();
    assert!(
        result.intersections.is_empty(),
        "Separated parallel planes should have no intersection"
    );
}

#[test]
fn anti_parallel_planes_separated_no_intersection() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 3.0),
        normal: Vec3::new(0.0, 0.0, -1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, -1.0, 0.0),
    };
    let result = intersect_surfaces(
        &s1,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &s2,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &test_tolerance(),
    )
    .unwrap();
    assert!(
        result.intersections.is_empty(),
        "Anti-parallel separated planes should have no intersection"
    );
}

// ── Near-coincident planes ─────────────────────────────────────────────

#[test]
fn near_coincident_planes_within_tolerance() {
    let tol = test_tolerance();
    let tiny_offset = tol.linear() * 0.1;
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, tiny_offset),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let result = intersect_surfaces(
        &s1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &s2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &tol,
    )
    .unwrap();
    assert!(
        result
            .intersections
            .iter()
            .any(|i| matches!(i, SurfaceSurfaceIntersection::CoincidentRegion { .. })),
        "Near-coincident planes should produce a CoincidentRegion"
    );
}

// ── Perpendicular planes with known intersection line ──────────────────

#[test]
fn perpendicular_planes_intersection_along_x_axis() {
    // XY plane (z=0) and XZ plane (y=0) intersect along the X axis.
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let result = intersect_surfaces(
        &s1,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &s2,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &test_tolerance(),
    )
    .unwrap();
    assert_eq!(result.intersections.len(), 1);
    let curve = match &result.intersections[0] {
        SurfaceSurfaceIntersection::Curve(c) => c,
        other => panic!("Expected Curve, got {:?}", other),
    };
    assert!(!curve.is_closed);

    // All points should have z~0 and y~0.
    use crate::geometry::curve::CurveEval;
    let (t0, t1) = curve.curve_3d.parameter_range();
    for i in 0..=4 {
        let t = t0 + (t1 - t0) * (i as f64) / 4.0;
        let p = curve.curve_3d.point_at(t);
        assert!(p.z.abs() < 0.01, "z={} at t={t}", p.z);
        assert!(p.y.abs() < 0.01, "y={} at t={t}", p.y);
    }
}

#[test]
fn angled_planes_45_degrees() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, -s, s),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, s, s),
    };
    let result = intersect_surfaces(
        &s1,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &s2,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &test_tolerance(),
    )
    .unwrap();
    assert_eq!(result.intersections.len(), 1);
    let curve = match &result.intersections[0] {
        SurfaceSurfaceIntersection::Curve(c) => c,
        other => panic!("Expected Curve, got {:?}", other),
    };
    use crate::geometry::curve::CurveEval;
    let (t0, t1) = curve.curve_3d.parameter_range();
    let mid = curve.curve_3d.point_at(0.5 * (t0 + t1));
    assert!(mid.y.abs() < 0.01, "y={}", mid.y);
    assert!(mid.z.abs() < 0.01, "z={}", mid.z);
}

// ── Tangent surface cases ──────────────────────────────────────────────

#[test]
fn cylinder_tangent_to_plane_no_curve() {
    // Plane at y=1 tangent to cylinder (radius 1, axis Z).
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let result = intersect_surfaces(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cyl,
        (0.0, TAU),
        (-5.0, 5.0),
        &test_tolerance(),
    )
    .unwrap();
    assert!(
        result.intersections.is_empty(),
        "Tangent cylinder-plane should be empty, got {} items",
        result.intersections.len()
    );
}

#[test]
fn sphere_tangent_to_plane_no_curve() {
    // Plane at z=1 tangent to sphere (radius 1, origin).
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let result = intersect_surfaces(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &test_tolerance(),
    )
    .unwrap();
    assert!(result.intersections.is_empty());
}

// ── Plane-cylinder intersection geometry ───────────────────────────────

#[test]
fn plane_perpendicular_to_cylinder_produces_circle() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 2.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let result = intersect_surfaces(
        &plane,
        (-3.0, 3.0),
        (-3.0, 3.0),
        &cyl,
        (0.0, TAU),
        (-5.0, 5.0),
        &test_tolerance(),
    )
    .unwrap();
    assert_eq!(result.intersections.len(), 1);
    let curve = match &result.intersections[0] {
        SurfaceSurfaceIntersection::Curve(c) => c,
        other => panic!("Expected Curve, got {:?}", other),
    };
    assert!(curve.is_closed);

    use crate::geometry::curve::CurveEval;
    let (t0, t1) = curve.curve_3d.parameter_range();
    for i in 0..=8 {
        let t = t0 + (t1 - t0) * (i as f64) / 8.0;
        let p = curve.curve_3d.point_at(t);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 1.0).abs() < 0.01, "r={r}");
        assert!((p.z - 2.0).abs() < 0.01, "z={}", p.z);
    }
}

#[test]
fn plane_parallel_to_cylinder_axis_two_lines() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let result = intersect_surfaces(
        &plane,
        (-3.0, 3.0),
        (-10.0, 10.0),
        &cyl,
        (0.0, TAU),
        (-10.0, 10.0),
        &test_tolerance(),
    )
    .unwrap();
    assert_eq!(result.intersections.len(), 2);
    for int in &result.intersections {
        let curve = match int {
            SurfaceSurfaceIntersection::Curve(c) => c,
            other => panic!("Expected Curve, got {:?}", other),
        };
        assert!(!curve.is_closed);
    }
}

// ── Plane-sphere intersection (offset) ─────────────────────────────────

#[test]
fn plane_sphere_offset_circle_radius() {
    // Sphere radius 2, plane at z=1. Circle radius = sqrt(4-1) = sqrt(3).
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
        axis: crate::geometry::Vec3::Z,
    };
    let result = intersect_surfaces(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &test_tolerance(),
    )
    .unwrap();
    assert_eq!(result.intersections.len(), 1);
    let curve = match &result.intersections[0] {
        SurfaceSurfaceIntersection::Curve(c) => c,
        other => panic!("Expected Curve, got {:?}", other),
    };
    assert!(curve.is_closed);

    let expected_r = 3.0_f64.sqrt();
    use crate::geometry::curve::CurveEval;
    let (t0, t1) = curve.curve_3d.parameter_range();
    for i in 0..=8 {
        let t = t0 + (t1 - t0) * (i as f64) / 8.0;
        let p = curve.curve_3d.point_at(t);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - expected_r).abs() < 0.05, "r={r:.4}");
        assert!((p.z - 1.0).abs() < 0.01, "z={:.4}", p.z);
    }
}

#[test]
fn sphere_no_intersection_when_plane_too_far() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 5.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let result = intersect_surfaces(
        &plane,
        (-3.0, 3.0),
        (-3.0, 3.0),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &test_tolerance(),
    )
    .unwrap();
    assert!(result.intersections.is_empty());
}

// ── Coincident spheres ────────────────────────────────────────────────

#[test]
fn coincident_spheres_return_region() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let result = intersect_surfaces(
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &sphere,
        (0.0, TAU),
        (0.0, PI),
        &test_tolerance(),
    )
    .unwrap();
    assert!(
        result
            .intersections
            .iter()
            .any(|i| matches!(i, SurfaceSurfaceIntersection::CoincidentRegion { .. })),
    );
}

// ── Disjoint parameter ranges ──────────────────────────────────────────

#[test]
fn perpendicular_planes_no_overlap_in_parameter_range() {
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
    let result = intersect_surfaces(
        &s1,
        (10.0, 20.0),
        (-1.0, 1.0),
        &s2,
        (-20.0, -10.0),
        (-1.0, 1.0),
        &test_tolerance(),
    )
    .unwrap();
    assert!(result.intersections.is_empty());
}
