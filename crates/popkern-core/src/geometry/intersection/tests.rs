#[cfg(test)]
use super::*;
use crate::geometry::Vec3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use std::f64::consts::{FRAC_PI_2, TAU};

fn line(ox: f64, oy: f64, dx: f64, dy: f64) -> CurveGeometry {
    CurveGeometry::Line {
        origin: Point3::new(ox, oy, 0.0),
        direction: Vec3::new(dx, dy, 0.0),
    }
}

#[test]
fn invalid_ranges_fail() {
    let l = line(0.0, 0.0, 1.0, 0.0);
    assert!(intersect_curves(&l, (1.0, 0.0), &l, (0.0, 1.0), &Tolerance::default()).is_err());
    assert!(intersect_curves(&l, (0.0, 1.0), &l, (1.0, 0.0), &Tolerance::default()).is_err());
}

#[test]
fn dispatch_line_line() {
    let r = intersect_curves(
        &line(0.0, 0.0, 1.0, 0.0),
        (0.0, 1.0),
        &line(0.5, -1.0, 0.0, 1.0),
        (0.0, 2.0),
        &Tolerance::default(),
    )
    .unwrap();
    assert_eq!(r.intersections.len(), 1);
}

#[test]
fn dispatch_line_circle_and_circle_line() {
    let l = line(-2.0, 0.0, 1.0, 0.0);
    let c = make_circle(0.0, 0.0);
    let r1 = intersect_curves(&l, (-2.0, 4.0), &c, (0.0, TAU), &Tolerance::default()).unwrap();
    assert_eq!(r1.intersections.len(), 2);
    let r2 = intersect_curves(&c, (0.0, TAU), &l, (-2.0, 4.0), &Tolerance::default()).unwrap();
    assert_eq!(r2.intersections.len(), 2);
}

#[test]
fn dispatch_arc_line_filters_range() {
    let arc = make_arc(0.0, 0.0, 0.0, FRAC_PI_2);
    let l = line(-2.0, 0.5, 1.0, 0.0);
    let r = intersect_curves(
        &arc,
        (0.0, FRAC_PI_2),
        &l,
        (-2.0, 4.0),
        &Tolerance::default(),
    )
    .unwrap();
    assert_eq!(r.intersections.len(), 1);
}

fn nurbs_line(p0: Point3, p1: Point3) -> CurveGeometry {
    CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    })
}

#[test]
fn dispatch_nurbs_fallback() {
    let c1 = CurveGeometry::Nurbs(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 3],
    });
    let c2 = nurbs_line(Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0));
    let r = intersect_curves(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default()).unwrap();
    assert!(r.intersections.len() >= 2);
}

#[test]
fn invariant_points_on_both_curves() {
    let c1 = make_circle(0.0, 0.0);
    let c2 = make_circle(1.0, 0.0);
    let result = intersect_curves(&c1, (0.0, TAU), &c2, (0.0, TAU), &Tolerance::default()).unwrap();
    for ix in &result.intersections {
        if let CurveIntersection::Point { t1, t2, .. } = ix {
            let p1 = c1.point_at(*t1);
            let p2 = c2.point_at(*t2);
            assert!(p1.distance_to(&p2) < 1e-6, "Invariant violated");
        }
    }
}

#[test]
fn angle_in_arc_range_cases() {
    assert!(angle_in_arc_range(0.5, 0.0, FRAC_PI_2));
    assert!(!angle_in_arc_range(2.0, 0.0, FRAC_PI_2));
    // Wrapping: [3pi/2, pi/2] wraps around 0
    assert!(angle_in_arc_range(0.0, 3.0 * FRAC_PI_2, FRAC_PI_2));
    assert!(!angle_in_arc_range(
        std::f64::consts::PI,
        3.0 * FRAC_PI_2,
        FRAC_PI_2
    ));
}

#[test]
fn swap_parameters_point_and_overlap() {
    let pt = swap_parameters(vec![CurveIntersection::Point {
        point: Point3::new(1.0, 2.0, 3.0),
        t1: 0.25,
        t2: 0.75,
    }]);
    if let CurveIntersection::Point { t1, t2, .. } = &pt[0] {
        assert!((*t1 - 0.75).abs() < 1e-15 && (*t2 - 0.25).abs() < 1e-15);
    }
    let ov = swap_parameters(vec![CurveIntersection::Overlap {
        t1_range: (0.1, 0.5),
        t2_range: (0.6, 0.9),
    }]);
    if let CurveIntersection::Overlap { t1_range, t2_range } = &ov[0] {
        assert!((t1_range.0 - 0.6).abs() < 1e-15 && (t2_range.0 - 0.1).abs() < 1e-15);
    }
}

fn make_arc(cx: f64, cy: f64, s: f64, e: f64) -> CurveGeometry {
    CurveGeometry::Arc {
        center: Point3::new(cx, cy, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: s,
        end_angle: e,
    }
}

fn make_circle(cx: f64, cy: f64) -> CurveGeometry {
    CurveGeometry::Circle {
        center: Point3::new(cx, cy, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    }
}

#[test]
fn dispatch_line_arc_and_arc_line() {
    let line = CurveGeometry::Line {
        origin: Point3::new(-2.0, 0.5, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let arc = make_arc(0.0, 0.0, 0.0, FRAC_PI_2);
    let r1 = intersect_curves(
        &line,
        (-2.0, 4.0),
        &arc,
        (0.0, FRAC_PI_2),
        &Tolerance::default(),
    )
    .unwrap();
    assert_eq!(r1.intersections.len(), 1);
    // Arc-Line swap
    let r2 = intersect_curves(
        &arc,
        (0.0, FRAC_PI_2),
        &line,
        (-2.0, 4.0),
        &Tolerance::default(),
    )
    .unwrap();
    assert_eq!(r2.intersections.len(), 1);
}

#[test]
fn dispatch_circle_arc_and_arc_circle() {
    let circle = make_circle(0.0, 0.0);
    let arc = make_arc(1.0, 0.0, 0.0, TAU);
    let r1 =
        intersect_curves(&circle, (0.0, TAU), &arc, (0.0, TAU), &Tolerance::default()).unwrap();
    assert!(r1.intersections.len() <= 2);
    // Arc-Circle swap
    let r2 =
        intersect_curves(&arc, (0.0, TAU), &circle, (0.0, TAU), &Tolerance::default()).unwrap();
    assert!(r2.intersections.len() <= 2);
}

#[test]
fn dispatch_arc_arc() {
    let pi = std::f64::consts::PI;
    let a1 = make_arc(0.0, 0.0, 0.0, pi);
    let a2 = make_arc(1.0, 0.0, 0.0, pi);
    let r = intersect_curves(&a1, (0.0, pi), &a2, (0.0, pi), &Tolerance::default()).unwrap();
    assert!(r.intersections.len() <= 2);
}

#[test]
fn filter_arc_range_keeps_overlap() {
    let results = vec![CurveIntersection::Overlap {
        t1_range: (0.1, 0.5),
        t2_range: (0.2, 0.6),
    }];
    assert_eq!(filter_arc_range(results, true, 0.0, 0.3).len(), 1);
}
