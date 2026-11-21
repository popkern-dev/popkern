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
