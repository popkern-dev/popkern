use super::*;
use std::f64::consts::{FRAC_PI_4, TAU};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// Cone: apex at origin, axis = +Z, half_angle = π/4 (45°)
fn std_cone() -> (Point3, Vec3, f64) {
    (
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_4,
    )
}

#[test]
fn secant_two_intersections() {
    // Line through the cone at two points.
    // Cone: apex=(0,0,0), axis=Z, half_angle=π/4 → r = v at height v.
    // At v=1, r=1; at v=2, r=2.
    // Line: origin=(−3, 1, 1), direction=(1, 0, 0).
    // Cone condition: x² + y² = z² (since tan(π/4)=1).
    // x² + 1 = 1 → x = 0, but we want intersections with the line y=1, z=1:
    // (x²+ 1) = 1 → x=0 only, that's one point. Let's use a better setup.
    //
    // Line: origin=(−3, 0, 2), dir=(1, 0, 0), z=2 throughout.
    // Cone at z=2: r=2. So intersection where x²+0²=4 → x=±2.
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-3.0, 0.0, 2.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    assert_eq!(
        results.len(),
        2,
        "expected 2 intersections, got {:?}",
        results
    );
    // Verify the points lie on the cone: x²+y²=z²
    for r in &results {
        if let SurfaceCurveIntersection::Point { point, .. } = r {
            let radial = (point.x * point.x + point.y * point.y).sqrt();
            assert!(
                (radial - point.z.abs()).abs() < 1e-6,
                "point not on cone: {:?}",
                point
            );
        }
    }
}

#[test]
fn tangent_one_intersection() {
    // Line tangent to cone surface.
    // Cone apex=(0,0,0), axis=Z, half_angle=π/4.
    // Line in the plane y=1, z=1: origin=(0,1,1), dir=(1,0,0).
    // Cone at z=1: r=1. y=1 is exactly the tangent point at x=0, y=1.
    // x²+1 = 1 → x=0 only (disc=0).
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(0.0, 1.0, 1.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    assert_eq!(
        results.len(),
        1,
        "expected 1 tangent intersection, got {:?}",
        results
    );
}

#[test]
fn miss_no_intersection() {
    // Line that misses the cone entirely.
    // Cone at z=1: r=1. Line at z=1, y=3 — too far from axis.
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-5.0, 3.0, 1.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    assert!(
        results.is_empty(),
        "expected no intersection, got {:?}",
        results
    );
}

#[test]
fn behind_apex_filtered() {
    // Line intersects the backward nappe (v < 0), should be filtered.
    // Line at z=−1: origin=(−2, 0, −1), dir=(1,0,0).
    // z=−1 < 0 → v = Q·Z = −1 < 0 → filtered.
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-2.0, 0.0, -1.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    assert!(
        results.is_empty(),
        "expected no intersection (behind apex), got {:?}",
        results
    );
}

#[test]
fn through_apex() {
    // Line through the apex — should produce either 0 or 1 result at the apex.
    // Line: origin=(0,0,−1), dir=(0,0,1) along axis. At apex Q·A=0 → tangential.
    // Degenerate: the line is along the axis. All points on axis satisfy Q perpendicular=0,
    // but only at apex itself does the cone equation hold (v=0, r=0).
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(0.0, 0.0, -1.0);
    let dir = Vec3::new(0.0, 0.0, 1.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    // Acceptable: 0 or 1 intersection at the apex (v=0)
    assert!(
        results.len() <= 1,
        "expected 0 or 1 apex intersections, got {:?}",
        results
    );
}

#[test]
fn surface_params_match_point_at() {
    use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-3.0, 0.0, 2.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    let cone = SurfaceGeometry::Cone {
        apex,
        axis,
        half_angle: ha,
    };
    for r in &results {
        if let SurfaceCurveIntersection::Point { point, u, v, .. } = r {
            let surface_pt = cone.point_at(*u, *v);
            assert!(
                point.distance_to(&surface_pt) < 1e-5,
                "surface params don't reproduce point: pt={:?}, surface_pt={:?}, u={}, v={}",
                point,
                surface_pt,
                u,
                v
            );
        }
    }
}

#[test]
fn degenerate_linear_case() {
    // When qa ≈ 0 (line lies on the cone surface), we get the linear case (line 56-60).
    // A line along a generator of the cone should produce a degenerate quadratic.
    let (apex, axis, ha) = std_cone();
    // Generator direction: at angle=0, the generator is (sin(ha), 0, cos(ha)) = (1/sqrt2, 0, 1/sqrt2)
    let s = ha.sin();
    let c = ha.cos();
    let dir = Vec3::new(s, 0.0, c);
    let origin = apex; // start at apex
    let results = intersect_line_cone(
        &origin,
        &dir,
        (0.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    // May return 0 or more points; the degenerate path was exercised
    let _ = results;
}

#[test]
fn t_range_filters_cone() {
    // Line intersects cone but t_range excludes the hits (line 83)
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-3.0, 0.0, 2.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (10.0, 20.0), // far beyond the intersection points
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    assert!(results.is_empty());
}

#[test]
fn v_range_filters_cone() {
    // Line intersects cone at z=2 but v_range excludes it (line 96)
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-3.0, 0.0, 2.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (0.0, TAU),
        (5.0, 10.0), // v_range excludes v=2
        &default_tol(),
    );
    assert!(results.is_empty());
}

#[test]
fn u_range_filters_cone() {
    // Line intersects cone but u_range excludes the azimuth (line 108)
    let (apex, axis, ha) = std_cone();
    let origin = Point3::new(-3.0, 0.0, 2.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_cone(
        &origin,
        &dir,
        (-10.0, 10.0),
        &apex,
        &axis,
        ha,
        (1.0, 2.0), // restricted u_range
        (0.0, f64::INFINITY),
        &default_tol(),
    );
    // Should filter out at least some intersections
    assert!(results.len() < 2);
}

#[test]
fn a_line_tangent_to_a_kilometre_cone_touches_it() {
    // A dimple's cone, apex down at (16, 10, 1) km-scaled, its rim of
    // radius 4000 at z = 5000; the plate edge x = 20000, z = 5000
    // touches the rim at (20000, 10000, 5000). The discriminant band
    // in squared units read it as a miss (PK-1143).
    let s = 1e3;
    let tol = Tolerance::default().with_linear(1e-6 * s);
    let hits = intersect_line_cone(
        &Point3::new(20.0 * s, 0.0, 5.0 * s),
        &Vec3::new(0.0, 20.0 * s, 0.0),
        (0.0, 1.0),
        &Point3::new(16.0 * s, 10.0 * s, 1.0 * s),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_4,
        (0.0, TAU),
        (0.0, 4.0 * 2f64.sqrt() * s),
        &tol,
    );
    assert_eq!(hits.len(), 1, "{hits:?}");
    let SurfaceCurveIntersection::Point { point, t, .. } = hits[0] else {
        panic!("{hits:?}");
    };
    assert!((t - 0.5).abs() < 1e-12, "{t}");
    assert!(point.distance_to(&Point3::new(20.0 * s, 10.0 * s, 5.0 * s)) < 1e-9 * s);
}

#[test]
fn a_line_just_outside_a_cone_misses_and_just_inside_touches() {
    let (apex, axis, ha) = std_cone();
    let tol = Tolerance::default();
    let at = |y: f64| {
        intersect_line_cone(
            &Point3::new(-5.0, y, 1.0),
            &Vec3::new(1.0, 0.0, 0.0),
            (-10.0, 10.0),
            &apex,
            &axis,
            ha,
            (0.0, TAU),
            (0.0, f64::INFINITY),
            &tol,
        )
        .len()
    };
    // At z = 1 the radius is 1; the cone's normal there leans 45
    // degrees, so a line at y = 1 + e is e/sqrt2 off the surface.
    assert_eq!(at(1.0 + 1e-9), 0);
    assert_eq!(at(1.0), 1);
    assert_eq!(at(1.0 - 1e-7), 1);
    assert_eq!(at(1.0 - 1e-4), 2);
}
