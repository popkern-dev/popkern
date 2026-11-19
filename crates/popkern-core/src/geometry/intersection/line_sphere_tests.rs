//! Tests for `line_sphere` (split out of the module file).
use super::*;
use std::f64::consts::{PI, TAU};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// Unit sphere at origin.
fn unit_sphere() -> (Point3, f64) {
    (Point3::new(0.0, 0.0, 0.0), 1.0)
}

#[test]
fn through_center_two_intersections() {
    // Line along X through centre of unit sphere: hits at x=−1 and x=+1.
    let (center, radius) = unit_sphere();
    let origin = Point3::new(-3.0, 0.0, 0.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_sphere(
        &origin,
        &dir,
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert_eq!(
        results.len(),
        2,
        "expected 2 intersections, got {:?}",
        results
    );
    // t values should be 2 and 4 (origin at x=−3, hits x=−1 at t=2, x=+1 at t=4)
    let ts: Vec<f64> = results
        .iter()
        .map(|r| match r {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            _ => panic!("unexpected overlap"),
        })
        .collect();
    assert!((ts[0] - 2.0).abs() < 1e-6, "t0 expected 2.0, got {}", ts[0]);
    assert!((ts[1] - 4.0).abs() < 1e-6, "t1 expected 4.0, got {}", ts[1]);
}

/// A line tangent to a ball a thousand units from the origin (an
/// ice-cream cone's seam ruling, from its apex 25 below the centre of
/// a ball of radius 15, through the rim point): one hit, at the
/// touch. Judged on the discriminant against ε², rounding gave two
/// hits about 2e-6 apart (PK-1133).
#[test]
fn a_tangent_line_far_from_the_origin_touches_once() {
    let centre = Point3::new(1234.5, -678.25, 910.125);
    for turn in [0.0f64, 0.7, 2.1] {
        let (s, c) = (turn.sin(), turn.cos());
        let apex = centre + Vec3::new(0.0, 0.0, -25.0);
        let touch = centre + Vec3::new(12.0 * c, 12.0 * s, -9.0);
        let dir = touch - apex;
        let hits = intersect_line_sphere(
            &apex,
            &dir,
            (0.0, 2.0),
            &centre,
            15.0,
            &Vec3::Z,
            (0.0, TAU),
            (0.0, PI),
            &default_tol(),
        );
        assert_eq!(hits.len(), 1, "turn {turn}: {hits:?}");
        let SurfaceCurveIntersection::Point { point, .. } = &hits[0] else {
            panic!("a point");
        };
        assert!(point.distance_to(&touch) < 1e-9, "{point:?}");
    }
}

#[test]
fn tangent_one_intersection() {
    // Line tangent to unit sphere: origin=(−2, 1, 0), dir=(1,0,0) — touches at (0,1,0).
    let (center, radius) = unit_sphere();
    let origin = Point3::new(-2.0, 1.0, 0.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_sphere(
        &origin,
        &dir,
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert_eq!(
        results.len(),
        1,
        "expected 1 tangent intersection, got {:?}",
        results
    );
    if let SurfaceCurveIntersection::Point { point, .. } = &results[0] {
        assert!(
            (point.distance_to(&Point3::new(0.0, 1.0, 0.0))) < 1e-5,
            "tangent point wrong: {:?}",
            point
        );
    }
}

#[test]
fn miss_no_intersection() {
    // Line that misses the sphere: origin=(−2, 2, 0), dir=(1,0,0) — y=2 > radius=1.
    let (center, radius) = unit_sphere();
    let origin = Point3::new(-2.0, 2.0, 0.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_sphere(
        &origin,
        &dir,
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert!(
        results.is_empty(),
        "expected no intersection, got {:?}",
        results
    );
}

#[test]
fn uv_params_match_point_at() {
    use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let (center, radius) = unit_sphere();
    let origin = Point3::new(-3.0, 0.0, 0.5);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_sphere(
        &origin,
        &dir,
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    let sphere = SurfaceGeometry::Sphere {
        center,
        radius,
        axis: crate::geometry::Vec3::Z,
    };
    assert!(!results.is_empty(), "expected intersections");
    for r in &results {
        if let SurfaceCurveIntersection::Point { point, u, v, .. } = r {
            let surface_pt = sphere.point_at(*u, *v);
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
fn north_pole_v_near_zero() {
    // Line along Z through north pole: should hit at z=+1 (v≈0) and z=−1 (v≈π).
    let (center, radius) = unit_sphere();
    let origin = Point3::new(0.0, 0.0, -3.0);
    let dir = Vec3::new(0.0, 0.0, 1.0);
    let results = intersect_line_sphere(
        &origin,
        &dir,
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert_eq!(
        results.len(),
        2,
        "expected 2 intersections, got {:?}",
        results
    );
    let vs: Vec<f64> = results
        .iter()
        .map(|r| match r {
            SurfaceCurveIntersection::Point { v, .. } => *v,
            _ => panic!("overlap unexpected"),
        })
        .collect();
    // First hit: south pole v=π, second hit: north pole v=0
    assert!(
        vs[0] > PI - 0.01,
        "first hit should be south pole, v={}",
        vs[0]
    );
    assert!(vs[1] < 0.01, "second hit should be north pole, v={}", vs[1]);
}

#[test]
fn offset_sphere_center() {
    // Sphere at (5,0,0) radius 2, line along X.
    let center = Point3::new(5.0, 0.0, 0.0);
    let radius = 2.0;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    let results = intersect_line_sphere(
        &origin,
        &dir,
        (0.0, 20.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert_eq!(
        results.len(),
        2,
        "expected 2 intersections, got {:?}",
        results
    );
    let ts: Vec<f64> = results
        .iter()
        .map(|r| match r {
            SurfaceCurveIntersection::Point { t, .. } => *t,
            _ => panic!(),
        })
        .collect();
    assert!((ts[0] - 3.0).abs() < 1e-6, "t0 expected 3.0, got {}", ts[0]);
    assert!((ts[1] - 7.0).abs() < 1e-6, "t1 expected 7.0, got {}", ts[1]);
}

#[test]
fn zero_direction_returns_empty() {
    // Zero-length direction vector should return no intersections (line 43)
    let (center, radius) = unit_sphere();
    let results = intersect_line_sphere(
        &Point3::new(1.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 0.0),
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert!(results.is_empty());
}

#[test]
fn t_range_filters_out() {
    // Line through sphere but t_range excludes the intersections (line 63)
    let (center, radius) = unit_sphere();
    let results = intersect_line_sphere(
        &Point3::new(-3.0, 0.0, 0.0),
        &Vec3::new(1.0, 0.0, 0.0),
        (5.0, 10.0), // intersections at t=2 and t=4, but range starts at 5
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, PI),
        &default_tol(),
    );
    assert!(results.is_empty());
}

#[test]
fn v_range_filters_out() {
    // Restrict v_range to exclude south pole (line 79)
    let (center, radius) = unit_sphere();
    // Line along Z axis hits north pole (v=0) and south pole (v=PI)
    let results = intersect_line_sphere(
        &Point3::new(0.0, 0.0, -3.0),
        &Vec3::new(0.0, 0.0, 1.0),
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.0, TAU),
        (0.0, 0.1), // v_range excludes most of sphere
        &default_tol(),
    );
    // Only the north pole (v=0) should pass the v_range filter
    assert!(results.len() <= 1);
}

#[test]
fn u_range_filters_out() {
    // Line through sphere with restricted u_range (line 84/102)
    let (center, radius) = unit_sphere();
    // Line along X from -3 through sphere. Hits at (-1,0,0) u=PI and (1,0,0) u=0
    let results = intersect_line_sphere(
        &Point3::new(-3.0, 0.0, 0.0),
        &Vec3::new(1.0, 0.0, 0.0),
        (-10.0, 10.0),
        &center,
        radius,
        &Vec3::Z,
        (0.5, 2.5), // restricted u_range, excludes u=0 and u=PI
        (0.0, PI),
        &default_tol(),
    );
    assert!(results.is_empty());
}
