use super::*;
use std::f64::consts::{FRAC_PI_2, PI};

#[test]
fn line_point_at() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let p = line.point_at(0.5);
    assert!((p.x - 0.5).abs() < 1e-10);
    assert!((p.y).abs() < 1e-10);
}

#[test]
fn line_derivative() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(2.0, 0.0, 0.0),
    };
    let d = line.derivative_at(0.5, 1);
    assert!((d.x - 2.0).abs() < 1e-10);
    let d2 = line.derivative_at(0.5, 2);
    assert!((d2.length()).abs() < 1e-10);
}

#[test]
fn line_length() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!((line.length_between(0.0, 1.0) - 1.0).abs() < 1e-10);
    assert!((line.length_between(0.0, 0.5) - 0.5).abs() < 1e-10);
}

#[test]
fn circle_point_at_start() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let p = circle.point_at(0.0);
    let dist = ((p.x - 0.0).powi(2) + (p.y - 0.0).powi(2)).sqrt();
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn circle_point_at_quarter() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let p = circle.point_at(FRAC_PI_2);
    let dist = ((p.x).powi(2) + (p.y).powi(2)).sqrt();
    assert!((dist - 1.0).abs() < 1e-10);
}

#[test]
fn circle_length_full() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let len = circle.length_between(0.0, 2.0 * PI);
    assert!((len - 4.0 * PI).abs() < 1e-10);
}

#[test]
fn arc_length() {
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let len = arc.length_between(0.0, FRAC_PI_2);
    assert!((len - FRAC_PI_2).abs() < 1e-10);
}

#[test]
fn line_tangent_is_normalized() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(3.0, 4.0, 0.0),
    };
    let t = line.tangent_at(0.0);
    assert!((t.length() - 1.0).abs() < 1e-10);
}

#[test]
fn line_parameter_range() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let (t0, t1) = line.parameter_range();
    assert!(t0.is_infinite() && t0 < 0.0);
    assert!(t1.is_infinite() && t1 > 0.0);
}

#[test]
fn circle_all_curve_eval_methods() {
    let c = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    assert!((c.derivative_at(0.0, 0).length() - 2.0).abs() < 1e-10);
    assert!((c.derivative_at(0.0, 1).length() - 2.0).abs() < 1e-10);
    assert!((c.derivative_at(0.0, 2).length() - 2.0).abs() < 1e-10);
    assert!((c.derivative_at(FRAC_PI_2, 3).length() - 2.0).abs() < 1e-10);
    assert!((c.tangent_at(0.0).length() - 1.0).abs() < 1e-10);
    let (t0, t1) = c.parameter_range();
    assert!((t0).abs() < 1e-10);
    assert!((t1 - std::f64::consts::TAU).abs() < 1e-10);
    assert!((c.length_between(0.0, FRAC_PI_2) - FRAC_PI_2 * 2.0).abs() < 1e-10);
}

#[test]
fn arc_all_curve_eval_methods() {
    let a = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.5,
        end_angle: 2.0,
    };
    let p = a.point_at(0.0);
    assert!(((p.x).powi(2) + (p.y).powi(2)).sqrt() - 1.0 < 1e-10);
    assert!((a.derivative_at(0.0, 0).length() - 1.0).abs() < 1e-10);
    assert!((a.derivative_at(0.0, 1).length() - 1.0).abs() < 1e-10);
    assert!((a.derivative_at(0.0, 2).length() - 1.0).abs() < 1e-10);
    assert!((a.derivative_at(FRAC_PI_2, 3).length() - 1.0).abs() < 1e-10);
    assert!((a.tangent_at(FRAC_PI_2).length() - 1.0).abs() < 1e-10);
    let (t0, t1) = a.parameter_range();
    assert!((t0 - 0.5).abs() < 1e-10);
    assert!((t1 - 2.0).abs() < 1e-10);
    assert!((a.length_between(0.0, PI) - PI).abs() < 1e-10);
}

#[test]
fn nurbs_dispatch_all_methods() {
    use crate::geometry::nurbs_curve::NurbsCurve;
    let nurbs = CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    });
    assert!((nurbs.point_at(0.5).x - 0.5).abs() < 1e-6);
    assert!(nurbs.derivative_at(0.5, 1).length() > 0.5);
    assert!((nurbs.tangent_at(0.5).length() - 1.0).abs() < 1e-6);
    let (t0, t1) = nurbs.parameter_range();
    assert!((t0).abs() < 1e-10 && (t1 - 1.0).abs() < 1e-10);
    assert!(nurbs.length_between(0.0, 1.0) > 0.9);
}

#[test]
fn arbitrary_perpendicular_both_branches() {
    let perp1 = arbitrary_perpendicular(Vec3::new(0.0, 0.0, 1.0));
    assert!(perp1.dot(&Vec3::new(0.0, 0.0, 1.0)).abs() < 1e-10);
    let perp2 = arbitrary_perpendicular(Vec3::new(1.0, 0.0, 0.0));
    assert!(perp2.dot(&Vec3::new(1.0, 0.0, 0.0)).abs() < 1e-10);
}
