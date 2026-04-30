//! Tests for [`super::pcurve`].  Extracted to keep the implementation
//! file under the 500-line structural cap.

use super::*;

#[test]
fn linear_endpoint_evaluation_is_exact() {
    let p = Pcurve::linear((0.0, 0.5), (std::f64::consts::PI, 0.5), 0.0, 1.0);
    let p0 = p.point_at(0.0);
    let p1 = p.point_at(1.0);
    assert_eq!(p0.x, 0.0);
    assert_eq!(p0.y, 0.5);
    assert_eq!(p1.x, std::f64::consts::PI);
    assert_eq!(p1.y, 0.5);
}

#[test]
fn linear_midpoint_is_average() {
    let p = Pcurve::linear((0.0, 0.0), (2.0, 4.0), 0.0, 1.0);
    let mid = p.point_at(0.5);
    assert_eq!(mid.x, 1.0);
    assert_eq!(mid.y, 2.0);
}

#[test]
fn linear_tangent_is_constant_direction() {
    let p = Pcurve::linear((1.0, 2.0), (5.0, 4.0), 0.0, 1.0);
    let t0 = p.tangent_at(0.25);
    let t1 = p.tangent_at(0.75);
    assert_eq!(t0.x, 4.0);
    assert_eq!(t0.y, 2.0);
    assert_eq!(t1.x, 4.0);
    assert_eq!(t1.y, 2.0);
}

#[test]
fn from_nurbs_round_trip() {
    let n = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 2.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let p: Pcurve = n.clone().into();
    assert!(p.as_nurbs().is_some());
    let back = p.to_nurbs();
    assert_eq!(back.control_points, n.control_points);
}

#[test]
fn linear_to_nurbs_preserves_endpoints() {
    let p = Pcurve::linear((0.0, 1.0), (4.0, 1.0), 0.0, 1.0);
    let n = p.to_nurbs();
    assert_eq!(n.degree, 1);
    assert_eq!(n.control_points.len(), 2);
    assert_eq!(n.control_points[0], Point3::new(0.0, 1.0, 0.0));
    assert_eq!(n.control_points[1], Point3::new(4.0, 1.0, 0.0));
}

#[test]
fn is_analytic_distinguishes_variants() {
    let lin = Pcurve::linear((0.0, 0.0), (1.0, 1.0), 0.0, 1.0);
    assert!(lin.is_analytic());

    let n = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let p: Pcurve = n.into();
    assert!(!p.is_analytic());
}
