use super::*;

#[test]
fn line_to_nurbs_in_range() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let nurbs = line.to_nurbs_in_range((-2.0, 4.0));
    assert_eq!(nurbs.degree, 1);
    for t in [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 3.99] {
        let p1 = line.point_at(t);
        let p2 = nurbs.point_at(t);
        assert!((p1.x - p2.x).abs() < 1e-10, "x differs at t={t}");
    }
}

#[test]
fn circle_to_nurbs_shape() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let nurbs = circle.to_nurbs();

    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    for t in [0.0, FRAC_PI_2, PI, PI + FRAC_PI_2] {
        let p1 = circle.point_at(t);
        let p2 = nurbs.point_at(t);
        assert!(
            (p1.x - p2.x).abs() < 1e-10,
            "x differs at t={t}: {} vs {}",
            p1.x,
            p2.x
        );
        assert!(
            (p1.y - p2.y).abs() < 1e-10,
            "y differs at t={t}: {} vs {}",
            p1.y,
            p2.y
        );
    }

    for i in 0..20 {
        let t = TAU * (i as f64) / 20.0;
        let p = nurbs.point_at(t);
        let dist = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (dist - 1.0).abs() < 1e-10,
            "point not on circle at t={t}: dist={}",
            dist
        );
    }
}

#[test]
fn nurbs_to_nurbs_is_clone() {
    let nurbs = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let geom = CurveGeometry::Nurbs(nurbs.clone());
    let promoted = geom.to_nurbs();
    assert_eq!(promoted.degree, 2);
    assert_eq!(promoted.control_points.len(), 3);
}

#[test]
fn nurbs_to_nurbs_in_range_honors_edge_span() {
    let nurbs = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 2.0, 2.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let geom = CurveGeometry::Nurbs(nurbs.clone());
    let trimmed = geom.to_nurbs_in_range((0.5, 1.5));
    let (trim_start, trim_end) = trimmed.parameter_range();

    assert!((trim_start - 0.5).abs() < 1e-12);
    assert!((trim_end - 1.5).abs() < 1e-12);
    assert!(
        trimmed
            .point_at(trim_start)
            .distance_to(&nurbs.point_at(0.5))
            < 1e-12
    );
    assert!(trimmed.point_at(trim_end).distance_to(&nurbs.point_at(1.5)) < 1e-12);
}

#[test]
fn nurbs_to_nurbs_in_range_preserves_reversed_edge_span() {
    let nurbs = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 2.0, 2.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let geom = CurveGeometry::Nurbs(nurbs.clone());
    let trimmed = geom.to_nurbs_in_range((1.5, 0.5));
    let (trim_start, trim_end) = trimmed.parameter_range();

    assert!((trim_start - 0.5).abs() < 1e-12);
    assert!((trim_end - 1.5).abs() < 1e-12);
    assert!(
        trimmed
            .point_at(trim_start)
            .distance_to(&nurbs.point_at(1.5))
            < 1e-12
    );
    assert!(trimmed.point_at(trim_end).distance_to(&nurbs.point_at(0.5)) < 1e-12);
}

#[test]
fn arc_to_nurbs_quarter_circle() {
    use std::f64::consts::FRAC_PI_2;
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let nurbs = arc.to_nurbs();
    let p = nurbs.point_at(nurbs.parameter_range().0);
    let r = (p.x * p.x + p.y * p.y).sqrt();
    assert!(
        (r - 1.0).abs() < 1e-6,
        "Start point not on unit circle: r={}",
        r
    );
}

#[test]
fn arc_to_nurbs_half_circle() {
    use std::f64::consts::PI;
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: PI,
    };
    let nurbs = arc.to_nurbs();
    assert!(nurbs.degree >= 2);
    assert!(!nurbs.control_points.is_empty());
}

#[test]
fn arc_to_nurbs_from_nonzero_start() {
    use std::f64::consts::{FRAC_PI_2, PI};
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: FRAC_PI_2,
        end_angle: PI,
    };
    let nurbs = arc.to_nurbs();
    assert!(nurbs.degree >= 2);
}

#[test]
fn arc_to_nurbs_full_circle_range() {
    use std::f64::consts::TAU;
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: TAU,
    };
    let nurbs = arc.to_nurbs();
    assert!(nurbs.degree >= 2);
}

#[test]
fn to_nurbs_default_range() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    let nurbs = line.to_nurbs();
    assert_eq!(nurbs.degree, 1);
}

#[test]
fn curve_on_edge_line_segment_is_bounded_analytic_line() {
    let coe = CurveOnEdge::line_segment(
        Point3::new(1.0, 2.0, 3.0),
        Point3::new(4.0, 6.0, 3.0),
        1.0e-9,
    )
    .expect("non-degenerate line segment");
    assert!(matches!(coe.geometry, CurveGeometry::Line { .. }));
    assert_eq!(coe.t_start, 0.0);
    assert_eq!(coe.t_end, 1.0);
    assert_eq!(
        coe.geometry.point_at(coe.t_start),
        Point3::new(1.0, 2.0, 3.0)
    );
    assert_eq!(coe.geometry.point_at(coe.t_end), Point3::new(4.0, 6.0, 3.0));
}

#[test]
fn arc_to_nurbs_in_range() {
    use std::f64::consts::FRAC_PI_2;
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let nurbs = arc.to_nurbs_in_range((0.0, FRAC_PI_2));
    assert!(nurbs.degree >= 2);
}

#[test]
fn arc_to_nurbs_boundary_split_fallback() {
    use std::f64::consts::TAU;
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: TAU - 0.001,
        end_angle: TAU,
    };
    let nurbs = arc.to_nurbs();
    assert!(nurbs.degree >= 2);
}
