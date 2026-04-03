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

// --- closest_point tests ---

#[test]
fn line_closest_point_perpendicular() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    // Point at (3, 5, 0): closest on x-axis is t=3
    let t = line.closest_point(&Point3::new(3.0, 5.0, 0.0)).unwrap();
    assert!((t - 3.0).abs() < 1e-10, "t = {t}");
}

#[test]
fn line_closest_point_negative_t() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let t = line.closest_point(&Point3::new(-2.0, 1.0, 0.0)).unwrap();
    assert!((t - (-2.0)).abs() < 1e-10, "t = {t}");
}

#[test]
fn circle_closest_point_exterior() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    // Point at (5, 0, 0): closest on unit circle should be at t≈0 (point (1,0,0))
    let t = circle.closest_point(&Point3::new(5.0, 0.0, 0.0)).unwrap();
    let closest = circle.point_at(t);
    let dist = closest.distance_to(&Point3::new(1.0, 0.0, 0.0));
    assert!(
        dist < 0.1,
        "closest should be near (1,0,0), got {:?}",
        closest
    );
}

#[test]
fn arc_closest_point() {
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    // Point at (3, 3, 0): closest should be near t≈PI/4 (point at 45 degrees)
    let t = arc.closest_point(&Point3::new(3.0, 3.0, 0.0)).unwrap();
    let closest = arc.point_at(t);
    let r = (closest.x.powi(2) + closest.y.powi(2)).sqrt();
    assert!((r - 1.0).abs() < 0.01, "should be on unit circle, r = {r}");
}

#[test]
fn nurbs_closest_point() {
    use crate::geometry::nurbs_curve::NurbsCurve;
    let nurbs = CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    });
    // Linear NURBS from (0,0,0) to (1,0,0). Point at (0.3, 1.0, 0.0).
    let t = nurbs.closest_point(&Point3::new(0.3, 1.0, 0.0)).unwrap();
    assert!((t - 0.3).abs() < 0.05, "t = {t}, expected ~0.3");
}

// --- curve_bounding_box tests ---

#[test]
fn line_bounding_box_finite_range() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 2.0, 3.0),
    };
    let bbox = line.curve_bounding_box((0.0, 1.0));
    // At t=0: (0,0,0), at t=1: (1,2,3)
    assert!(bbox.min.x <= 0.01);
    assert!(bbox.max.x >= 0.99);
    assert!(bbox.max.z >= 2.99);
}

#[test]
fn circle_bounding_box_full() {
    use std::f64::consts::TAU;
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let bbox = circle.curve_bounding_box((0.0, TAU));
    // Full circle of radius 2: bbox should be roughly [-2, 2] x [-2, 2] x [0, 0]
    assert!(bbox.min.x < -1.9, "min.x = {}", bbox.min.x);
    assert!(bbox.max.x > 1.9, "max.x = {}", bbox.max.x);
    assert!(bbox.min.y < -1.9, "min.y = {}", bbox.min.y);
    assert!(bbox.max.y > 1.9, "max.y = {}", bbox.max.y);
    assert!(bbox.min.z.abs() < 0.01);
    assert!(bbox.max.z.abs() < 0.01);
}

#[test]
fn arc_bounding_box_quarter() {
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let bbox = arc.curve_bounding_box((0.0, FRAC_PI_2));
    // Arc basis vectors for normal=(0,0,1): u=(0,1,0), v=(-1,0,0)
    // t=0 -> (0,1,0), t=PI/2 -> (-1,0,0), so bbox spans that region
    let diag = bbox.diagonal();
    assert!(diag > 0.5, "bbox should have nonzero extent, diag={diag}");
    assert!(diag < 3.0, "bbox should not be oversized, diag={diag}");
    // z should be near zero (arc in XY plane)
    assert!(bbox.min.z.abs() < 0.01);
    assert!(bbox.max.z.abs() < 0.01);
}
