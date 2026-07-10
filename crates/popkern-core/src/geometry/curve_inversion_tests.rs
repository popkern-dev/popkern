use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn line_curve(origin: Point3, direction: Vec3, t_start: f64, t_end: f64) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Line { origin, direction },
        t_start,
        t_end,
    }
}

fn arc_curve(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    }
}

fn circle_curve(center: Point3, normal: Vec3, radius: f64) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center,
            normal,
            radius,
        },
        t_start: 0.0,
        t_end: TAU,
    }
}

// ---------- Line tier 1 ----------

#[test]
fn line_inversion_at_origin() {
    let line = line_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        -10.0,
        10.0,
    );
    let t = invert_curve_at_point(&line, Point3::new(2.5, 0.0, 0.0), &tol()).unwrap();
    assert!((t - 2.5).abs() < 1e-12);
}

#[test]
fn line_inversion_off_curve() {
    let line = line_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        -10.0,
        10.0,
    );
    let result = invert_curve_at_point(&line, Point3::new(2.5, 0.5, 0.0), &tol());
    assert!(matches!(
        result,
        Err(PointInversionError::PointOffCurve { .. })
    ));
}

#[test]
fn line_inversion_out_of_range() {
    let line = line_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        0.0,
        1.0,
    );
    let result = invert_curve_at_point(&line, Point3::new(5.0, 0.0, 0.0), &tol());
    assert!(matches!(
        result,
        Err(PointInversionError::ParameterOutOfRange { .. })
    ));
}

// ---------- Arc tier 1 ----------

/// Reconstruct the expected point at a given parameter via `point_at`,
/// then invert and check round-trip.
fn check_round_trip(curve: &CurveOnEdge, t_expected: f64) {
    let p = curve.geometry.point_at(t_expected);
    let t = invert_curve_at_point(curve, p, &tol())
        .unwrap_or_else(|e| panic!("invert failed at t={t_expected}: {e:?}"));
    assert!(
        (t - t_expected).abs() < 1e-9,
        "round-trip t={t}, expected={t_expected}, point={p:?}"
    );
}

#[test]
fn arc_inversion_at_endpoints() {
    let arc = arc_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );
    check_round_trip(&arc, 0.0);
    check_round_trip(&arc, FRAC_PI_2);
}

#[test]
fn arc_inversion_at_midpoint() {
    let arc = arc_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );
    check_round_trip(&arc, FRAC_PI_4);
}

#[test]
fn arc_inversion_off_curve() {
    let arc = arc_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );
    // Off-plane: pick a point on the arc, then nudge it along normal.
    let p_on_arc = arc.geometry.point_at(FRAC_PI_4);
    let p_off_plane = Point3::new(p_on_arc.x, p_on_arc.y, 0.5);
    let result = invert_curve_at_point(&arc, p_off_plane, &tol());
    assert!(
        matches!(result, Err(PointInversionError::PointOffCurve { .. })),
        "off-plane should be off-curve, got {result:?}"
    );
    // In-plane but wrong radius (point at half the radius from centre).
    let half_radius = Point3::new(p_on_arc.x * 0.5, p_on_arc.y * 0.5, 0.0);
    let result = invert_curve_at_point(&arc, half_radius, &tol());
    assert!(
        matches!(result, Err(PointInversionError::PointOffCurve { .. })),
        "wrong radius should be off-curve, got {result:?}"
    );
}

#[test]
fn arc_inversion_out_of_arc_range() {
    let arc = arc_curve(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        0.0,
        FRAC_PI_2,
    );
    // Pick a point on the same circle but at parameter t = π
    // (outside the [0, π/2] arc range).
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let p_on_circle_outside_arc = circle.point_at(PI);
    let result = invert_curve_at_point(&arc, p_on_circle_outside_arc, &tol());
    assert!(
        matches!(result, Err(PointInversionError::ParameterOutOfRange { .. })),
        "expected ParameterOutOfRange, got {result:?}"
    );
}

// ---------- Circle tier 1 ----------

#[test]
fn circle_inversion_round_trip_at_quarter() {
    let circle = circle_curve(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0), 1.0);
    // Round-trip at multiple parameters.
    for t_expected in [0.0_f64, FRAC_PI_2, PI, FRAC_PI_2 * 3.0] {
        let p = circle.geometry.point_at(t_expected);
        let t = invert_curve_at_point(&circle, p, &tol()).unwrap();
        // For a full circle the result is wrapped into [0, 2π).
        // 0 and 2π are the same point; allow either.
        let diff = (t - t_expected).abs().min((t - t_expected - TAU).abs());
        assert!(
            diff < 1e-9,
            "circle round-trip: got {t}, expected {t_expected}"
        );
    }
}

// ---------- as_circular_arc detector ----------

#[test]
fn as_circular_arc_recognises_arc_directly() {
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let params = arc.as_circular_arc(&tol()).unwrap();
    assert!((params.radius - 2.0).abs() < 1e-12);
    assert!((params.start_angle - 0.0).abs() < 1e-12);
    assert!((params.sweep_angle - FRAC_PI_2).abs() < 1e-12);
}

#[test]
fn as_circular_arc_returns_none_for_line() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(line.as_circular_arc(&tol()).is_none());
}

#[test]
fn as_circular_arc_recognises_nurbs_quarter_arc() {
    // A rational quadratic NURBS quarter circle in xy-plane, radius 1.
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let nurbs = CurveGeometry::Nurbs(arc.to_nurbs());
    let params = nurbs
        .as_circular_arc(&tol())
        .expect("NURBS quarter arc should be detected as circular");
    assert!(
        (params.radius - 1.0).abs() < 1e-6,
        "fitted radius should be ≈1, got {}",
        params.radius
    );
    // The fitted normal should be parallel to (0, 0, 1) (sign may flip
    // depending on sweep direction).
    let z_align = params.normal.dot(&Vec3::new(0.0, 0.0, 1.0)).abs();
    assert!(
        z_align > 0.999,
        "normal should be ±z, got {:?}",
        params.normal
    );
}

#[test]
fn as_circular_arc_returns_none_for_non_arc_nurbs() {
    // Cubic NURBS through scattered (non-circular) points.
    let nurbs = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.5, 0.0),
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(1.5, 1.0, 0.0),
            Point3::new(2.0, 0.3, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
    };
    let geom = CurveGeometry::Nurbs(nurbs);
    assert!(geom.as_circular_arc(&tol()).is_none());
}

#[test]
fn invert_nurbs_arc_falls_back_to_arc_path() {
    let arc = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
        start_angle: 0.0,
        end_angle: FRAC_PI_2,
    };
    let nurbs_curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(arc.to_nurbs()),
        t_start: 0.0,
        t_end: FRAC_PI_2,
    };
    // Pick the point at the analytic arc's midpoint and confirm the
    // NURBS-detected-as-arc tier 1b path produces the same parameter
    // (up to small fitting / detector tolerance).
    let mid_point = arc.point_at(FRAC_PI_4);
    let result = invert_curve_at_point(&nurbs_curve, mid_point, &tol());
    let t = result.expect("NURBS-stored arc should invert via Tier 1b");
    assert!(
        (t - FRAC_PI_4).abs() < 1e-3,
        "NURBS-detected arc midpoint inversion: got t={t}, expected ≈π/4"
    );
}

#[test]
fn invert_nurbs_line_uses_analytic_line_path() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    };
    let t = invert_curve_at_point(&curve, Point3::new(2.0, 0.0, 0.0), &tol()).unwrap();
    assert!((t - 0.5).abs() < 1e-12);
}

#[test]
fn invert_generic_nurbs_finds_on_curve_points() {
    let nurbs = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.5, 0.0),
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(1.5, 1.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    };
    // Generic (scan + ternary) inversion: an on-curve point inverts
    // to its parameter; a far point reports PointOffCurve. This
    // replaced the Tier-2 stub — a boolean junction pre-split needed
    // inversion on exact-arc NURBS pieces whose circular-arc
    // reconstruction answers in a phase-shifted frame.
    use crate::geometry::curve::CurveEval as _;
    let probe = curve.geometry.point_at(0.37);
    let t = invert_curve_at_point(&curve, probe, &tol()).expect("on-curve point inverts");
    assert!(
        curve.geometry.point_at(t).distance_to(&probe) < 1.0e-5,
        "inverted parameter must reproduce the probe point"
    );
    let result = invert_curve_at_point(&curve, Point3::new(5.0, -3.0, 2.0), &tol());
    assert!(matches!(
        result,
        Err(PointInversionError::PointOffCurve { .. })
    ));
}

#[test]
fn ellipse_returns_unsupported() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center: Point3::new(0.0, 0.0, 0.0),
            u_basis: Vec3::new(2.0, 0.0, 0.0),
            v_basis: Vec3::new(0.0, 1.0, 0.0),
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        },
        t_start: 0.0,
        t_end: FRAC_PI_2,
    };
    let result = invert_curve_at_point(&curve, Point3::new(2.0, 0.0, 0.0), &tol());
    assert!(matches!(
        result,
        Err(PointInversionError::UnsupportedCurve { .. })
    ));
}

// ---------- circle_through_three_points ----------

#[test]
fn fit_circle_through_three_points_xy_plane() {
    let p0 = Point3::new(1.0, 0.0, 0.0);
    let p1 = Point3::new(0.0, 1.0, 0.0);
    let p2 = Point3::new(-1.0, 0.0, 0.0);
    let (center, radius, normal) = circle_through_three_points(p0, p1, p2).unwrap();
    assert!(center.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-12);
    assert!((radius - 1.0).abs() < 1e-12);
    assert!((normal.dot(&Vec3::new(0.0, 0.0, 1.0)) - 1.0).abs() < 1e-12);
}

#[test]
fn fit_circle_returns_none_for_collinear() {
    let p0 = Point3::new(0.0, 0.0, 0.0);
    let p1 = Point3::new(1.0, 0.0, 0.0);
    let p2 = Point3::new(2.0, 0.0, 0.0);
    assert!(circle_through_three_points(p0, p1, p2).is_none());
}
