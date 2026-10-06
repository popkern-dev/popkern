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

/// PK-257 terminal-sliver guard: an out-of-principal arc whose end
/// carries atan2-class noise (width π/2 + 2e-8) must not strand a
/// sub-tolerance nano segment against the quarter-knot grid — the
/// full-multiplicity break poisoned downstream knot unions (PK-255).
/// The re-anchored grid keeps every interior knot clear of both ends
/// while the endpoints stay exact.
#[test]
fn arc_to_nurbs_near_quarter_end_leaves_no_sliver_span() {
    use std::f64::consts::FRAC_PI_2;
    let noise = 2.107e-8;
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.05,
        start_angle: -(FRAC_PI_2 + noise),
        end_angle: 0.0,
    };
    let nurbs = arc.to_nurbs();
    let (t_lo, t_hi) = nurbs.parameter_range();
    let domain = (t_hi - t_lo).abs().max(f64::EPSILON);
    let min_span = nurbs
        .knots
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .filter(|span| *span > 0.0)
        .fold(f64::INFINITY, f64::min);
    assert!(
        min_span > domain * 1.0e-3,
        "arc conversion strands a sliver span: min_span={min_span:.3e} of domain {domain:.3e}"
    );
    // Knots are absolute angle: parameter equals angle AT the knots
    // (the rational quadratic's interior parameterization skews, as
    // it always has); the curve must lie ON the circle everywhere.
    let mut distinct = nurbs.knots.clone();
    distinct.dedup();
    for &knot in &distinct {
        let converted = nurbs.point_at(knot);
        let exact = arc.point_at(knot);
        assert!(
            converted.distance_to(&exact) < 1.0e-12,
            "knot point deviates at t={knot}: {converted:?} vs {exact:?}"
        );
    }
    for index in 0..=32 {
        let t = t_lo + (t_hi - t_lo) * index as f64 / 32.0;
        let p = nurbs.point_at(t);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (r - 0.05).abs() < 1.0e-12,
            "converted point off the circle at t={t}: r={r}"
        );
    }
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

/// PK-462: a closed curve's closest parameter comes back inside the
/// curve's own range. Newton can settle a whole period off (−0.389 −
/// 2π on a rim circle), and a caller clamping that into the window
/// read a point ON the curve as 3.3e-3 away.
#[test]
fn closest_point_on_a_closed_curve_stays_in_its_range() {
    use crate::geometry::curve::closest_point_on_curve;
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.003,
    };
    let (t0, t1) = (&circle as &dyn CurveEval).parameter_range();
    for k in 0..64 {
        let a = std::f64::consts::TAU * k as f64 / 64.0 + 0.0123;
        let p = Point3::new(0.003 * a.cos(), 0.003 * a.sin(), 1e-12);
        let Some(t) = closest_point_on_curve(&circle, &p) else {
            panic!("no closest point at {a}");
        };
        assert!(t >= t0 && t <= t1, "t {t} outside [{t0}, {t1}]");
        assert!(
            circle.point_at(t).distance_to(&p) < 1e-11,
            "point at {a} misplaced"
        );
    }
}
