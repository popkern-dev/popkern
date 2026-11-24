//! Tests for the overlap rescue seeding.
use crate::geometry::curve::CurveEval as _;
use crate::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use crate::geometry::surface::{SurfaceEval as _, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// `v` turned by `angle` about the unit `axis` (Rodrigues).
fn turn(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let (c, s) = (angle.pcos(), angle.psin());
    v * c + axis.cross(&v) * s + axis * (axis.dot(&v) * (1.0 - c))
}

/// Two cones touching externally at a point, their generators crossing
/// there, the second pushed `depth` into the first.
fn crossed_cones(depth: f64) -> (SurfaceGeometry, SurfaceGeometry) {
    let half: f64 = 0.3;
    let a = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: half,
    };
    // At azimuth 1 rad, off A's seam.
    let z = Vec3::new(0.0, 0.0, 1.0);
    let p = Point3::new(0.0, 0.0, 0.0)
        + turn(Vec3::new(2.0 * half.psin(), 0.0, 2.0 * half.pcos()), z, 1.0);
    let n = turn(Vec3::new(half.pcos(), 0.0, -half.psin()), z, 1.0);
    // Mirrored in the tangent plane at p, turned about the normal there.
    let mirror = |x: Point3| x - n * (2.0 * (x - p).dot(&n));
    let axis_m = Vec3::new(0.0, 0.0, 1.0) - n * (2.0 * n.z);
    let apex = p + turn(mirror(Point3::new(0.0, 0.0, 0.0)) - p, n, 0.7) - n * depth;
    let b = SurfaceGeometry::Cone {
        apex,
        axis: turn(axis_m, n, 0.7),
        half_angle: half,
    };
    (a, b)
}

/// Overlapping by ten tolerances at an external touch, two cones meet in
/// a loop about 9e-3 across that no starting-point grid held: the pair
/// emitted nothing and the common of two cone bodies was served empty
/// (W3 near-miss cells). The overlap's deepest point seeds it.
#[test]
fn cones_overlapping_by_ten_tolerances_meet_in_a_small_loop() {
    let tol = Tolerance::default();
    let (a, b) = crossed_cones(10.0 * tol.linear());
    let w = ((0.0, std::f64::consts::TAU), (0.0, 4.0));
    let r = intersect_surfaces(&a, w.0, w.1, &b, w.0, w.1, &tol).expect("intersect");
    let curves: Vec<_> = r
        .intersections
        .iter()
        .filter_map(|i| match i {
            SurfaceSurfaceIntersection::Curve(c) => Some(c),
            _ => None,
        })
        .collect();
    assert_eq!(curves.len(), 1, "one loop");
    let c = curves[0];
    assert!(c.is_closed, "closed");
    let (t0, t1) = c.curve_3d.parameter_range();
    let mut len = 0.0;
    let mut last = c.curve_3d.point_at(t0);
    for k in 0..=400 {
        let q = c.curve_3d.point_at(t0 + (t1 - t0) * k as f64 / 400.0);
        len += q.distance_to(&last);
        last = q;
        for s in [&a, &b] {
            let uv = s.closest_point(&q).expect("foot");
            let off = s.point_at(uv.0, uv.1).distance_to(&q);
            // Marched samples to the tolerance, the fit between them.
            assert!(off < 2.0 * tol.linear(), "on both cones: {off:.2e}");
        }
    }
    assert!(len < 0.1, "a small loop: {len:.3e}");
}

/// Apart by ten tolerances, the same cones meet nowhere: the rescue seeds
/// only overlaps.
#[test]
fn cones_apart_by_ten_tolerances_meet_nowhere() {
    let tol = Tolerance::default();
    let (a, b) = crossed_cones(-10.0 * tol.linear());
    let w = ((0.0, std::f64::consts::TAU), (0.0, 4.0));
    let r = intersect_surfaces(&a, w.0, w.1, &b, w.0, w.1, &tol).expect("intersect");
    assert!(r.intersections.is_empty(), "{:?}", r.intersections.len());
}

/// A ball touching a cone's wall from inside at a point (normals aligned,
/// the ball's radius under the wall's), pushed `depth` out through it.
fn ball_in_cone(depth: f64) -> (SurfaceGeometry, SurfaceGeometry) {
    let half: f64 = 0.3;
    let a = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: half,
    };
    let z = Vec3::new(0.0, 0.0, 1.0);
    let p = Point3::new(0.0, 0.0, 0.0)
        + turn(Vec3::new(2.0 * half.psin(), 0.0, 2.0 * half.pcos()), z, 1.0);
    let n = turn(Vec3::new(half.pcos(), 0.0, -half.psin()), z, 1.0);
    let r = 0.3;
    let b = SurfaceGeometry::Sphere {
        center: p - n * (r - depth),
        radius: r,
        axis: Vec3::new(0.0, 0.0, 1.0),
    };
    (a, b)
}

/// Pushed two tolerances out through the cone's wall from inside, a ball
/// meets it in a loop no starting-point grid holds: the gap's extremum there is a positive
/// maximum (one surface inside the other), which the rescue once skipped
/// for want of a negative gap (YW-I 1396: a W face poking out of a
/// cylinder met in no seed, and B - A served empty). Pulled two
/// tolerances in, it meets nothing.
#[test]
fn a_ball_poking_through_a_cone_from_inside_meets_it_in_a_small_loop() {
    let tol = Tolerance::default();
    let w = ((0.0, std::f64::consts::TAU), (0.0, 4.0));
    let ws = ((0.0, std::f64::consts::TAU), (0.0, std::f64::consts::PI));
    // The ball first, its gap read against the cone's outward normal
    // (positive where it pokes out), as the W face came first in 1396.
    let (a, b) = ball_in_cone(2.0 * tol.linear());
    let r = intersect_surfaces(&b, ws.0, ws.1, &a, w.0, w.1, &tol).expect("intersect");
    let closed: Vec<_> = r
        .intersections
        .iter()
        .filter_map(|i| match i {
            SurfaceSurfaceIntersection::Curve(c) if c.is_closed => Some(c),
            _ => None,
        })
        .collect();
    assert_eq!(closed.len(), 1, "one loop: {:?}", r.intersections.len());
    let (a, b) = ball_in_cone(-2.0 * tol.linear());
    let r = intersect_surfaces(&b, ws.0, ws.1, &a, w.0, w.1, &tol).expect("intersect");
    assert!(r.intersections.is_empty(), "{:?}", r.intersections.len());
}

/// Overlapping past the near-miss band (three hundred tolerances), the
/// same cones still meet in a loop no grid holds: the overlap's deepest
/// point is a critical point of the gap the contact search keeps beside
/// its near misses (SW-E 759: a ball sunk two thousandths into a dome
/// block met it in no seed, its common served empty).
#[test]
fn cones_overlapping_past_the_near_miss_band_meet_in_a_loop() {
    let tol = Tolerance::default();
    for depth in [300.0, 1000.0, 3000.0] {
        let (a, b) = crossed_cones(depth * tol.linear());
        let w = ((0.0, std::f64::consts::TAU), (0.0, 4.0));
        let r = intersect_surfaces(&a, w.0, w.1, &b, w.0, w.1, &tol).expect("intersect");
        let closed = r
            .intersections
            .iter()
            .filter(|i| matches!(i, SurfaceSurfaceIntersection::Curve(c) if c.is_closed))
            .count();
        assert_eq!(
            closed,
            1,
            "one loop at {depth} tolerances: {}",
            r.intersections.len()
        );
    }
}

/// A ball poking out of a cone, a thousandth the size: the chart
/// normals carry the charts' speeds (the ball's its radius squared), so
/// their cross product along the loop fell under the solver tolerance
/// though the surfaces cross at a tenth of a radian, and every march
/// step was taken for a tangency and stalled (SW-E 765: a ball against
/// a dome block at that scale fitted no loop, its common served empty).
#[test]
fn a_ball_poking_through_a_cone_a_thousandth_the_size_meets_it_in_a_loop() {
    let tol = Tolerance::default();
    let s = 1e-3;
    let (a, b) = ball_in_cone(10.0 * tol.linear() / s);
    let (a, b) = match (a, b) {
        (
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
            SurfaceGeometry::Sphere {
                center,
                radius,
                axis: sphere_axis,
            },
        ) => (
            SurfaceGeometry::Cone {
                apex: Point3::new(apex.x * s, apex.y * s, apex.z * s),
                axis,
                half_angle,
            },
            SurfaceGeometry::Sphere {
                center: Point3::new(center.x * s, center.y * s, center.z * s),
                radius: radius * s,
                axis: sphere_axis,
            },
        ),
        _ => unreachable!("a cone and a ball"),
    };
    let w = ((0.0, std::f64::consts::TAU), (0.0, 4.0 * s));
    let ws = ((0.0, std::f64::consts::TAU), (0.0, std::f64::consts::PI));
    let r = intersect_surfaces(&b, ws.0, ws.1, &a, w.0, w.1, &tol).expect("intersect");
    let closed = r
        .intersections
        .iter()
        .filter(|i| matches!(i, SurfaceSurfaceIntersection::Curve(c) if c.is_closed))
        .count();
    assert_eq!(closed, 1, "one loop: {}", r.intersections.len());
}
