use super::*;
use crate::geometry::NurbsCurve;
use crate::geometry::curve::CurveGeometry;

fn circle_nurbs(center: Point3, normal: Vec3, radius: f64) -> NurbsCurve {
    CurveGeometry::Circle {
        center,
        normal,
        radius,
    }
    .to_nurbs()
}

#[test]
fn a_rational_circle_is_recognised_with_its_radius_and_orientation() {
    let center = Point3::new(1.0, -2.0, 0.5);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let nurbs = circle_nurbs(center, normal, 0.18);
    let range = nurbs.parameter_range();
    let fit = recognize_circular(&nurbs, range, 1e-9).expect("circle");
    assert!(fit.closed);
    assert!((fit.radius - 0.18).abs() < 1e-12, "{}", fit.radius);
    assert!(fit.center.distance_to(&center) < 1e-12);
    assert!((fit.normal.dot(&normal).abs() - 1.0).abs() < 1e-12);
    assert!(fit.max_deviation < 1e-12);
    let arc = arc_on_edge(&nurbs, range, &fit);
    let CurveGeometry::Arc {
        start_angle,
        end_angle,
        ..
    } = arc.geometry
    else {
        panic!("arc");
    };
    assert!((end_angle - start_angle - std::f64::consts::TAU).abs() < 1e-12);
    // The arc starts where the curve starts and runs the same way.
    assert!(
        arc.geometry
            .point_at(arc.t_start)
            .distance_to(&nurbs.point_at(range.0))
            < 1e-12
    );
    let quarter = range.0 + (range.1 - range.0) * 0.25;
    let ahead = nurbs.point_at(quarter) - fit.center;
    let d0 = nurbs.point_at(range.0) - fit.center;
    assert!(d0.cross(&ahead).dot(&fit.normal) > 0.0, "right-hand sense");
}

#[test]
fn an_open_arc_sub_range_is_recognised_with_its_end_angles() {
    let nurbs = circle_nurbs(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), 2.0);
    let (t0, t1) = nurbs.parameter_range();
    let range = (t0, t0 + (t1 - t0) * 0.4);
    let fit = recognize_circular(&nurbs, range, 1e-9).expect("arc");
    assert!(!fit.closed);
    assert!((fit.radius - 2.0).abs() < 1e-9);
    let arc = arc_on_edge(&nurbs, range, &fit);
    assert!(
        arc.geometry
            .point_at(arc.t_end)
            .distance_to(&nurbs.point_at(range.1))
            < 1e-9
    );
    assert!(arc.t_end > arc.t_start && arc.t_end - arc.t_start < std::f64::consts::TAU);
}

#[test]
fn an_ellipse_and_a_line_are_not_circles() {
    let ellipse = CurveGeometry::Ellipse {
        center: Point3::new(0.0, 0.0, 0.0),
        u_basis: Vec3::new(2.0, 0.0, 0.0),
        v_basis: Vec3::new(0.0, 1.0, 0.0),
        start_angle: 0.0,
        end_angle: std::f64::consts::TAU,
    }
    .to_nurbs();
    assert!(recognize_circular(&ellipse, ellipse.parameter_range(), 1e-6).is_none());
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(recognize_circular(&line, (0.0, 1.0), 1e-6).is_none());
}

/// A straight run at a coarse gate is not a circle either: the
/// least-squares fit of collinear samples answers an enormous radius
/// whose sagitta sits under the gate (a bore's seam line read as a
/// circle of radius 8078 with a zero-degree sweep at 1e-3, PK-418).
#[test]
fn a_straight_nurbs_seam_is_not_a_circle_at_a_coarse_gate() {
    let line = CurveGeometry::Line {
        origin: Point3::new(0.2, 0.0, -0.566),
        direction: Vec3::new(0.0, 0.0, 1.0),
    }
    .to_nurbs();
    for gate in [1e-9, 1e-6, 1e-3, 1e-2] {
        assert!(
            recognize_circular(&line, line.parameter_range(), gate).is_none(),
            "gate {gate:e}"
        );
    }
    // A genuine shallow arc whose sagitta clears the gate still
    // recognises: radius 10, chord ~1.6 → sagitta 0.032 > 1e-3.
    let shallow = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 10.0,
        start_angle: 0.0,
        end_angle: 0.16,
    }
    .to_nurbs();
    assert!(recognize_circular(&shallow, shallow.parameter_range(), 1e-3).is_some());
}

/// PK-779: at a coarse gate a quarter arc whose chord is under 100
/// gates still determines its circle by its bulge (sagitta 0.21 against
/// a 1e-2 gate), while a flat arc with neither a long chord nor a bulge
/// over 10 gates stays unrecognised.
#[test]
fn a_quarter_rim_is_recognised_at_a_coarse_gate() {
    let quarter = crate::geometry::curve::CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.1),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.717,
        start_angle: 0.0,
        end_angle: std::f64::consts::FRAC_PI_2,
    }
    .to_nurbs();
    let range = quarter.parameter_range();
    let fit = recognize_circular(&quarter, range, 1e-2).expect("recognised");
    assert!((fit.radius - 0.717).abs() < 1e-6, "radius {}", fit.radius);
    // A 0.5-long arc of radius 50: sagitta 6.25e-4, under 10 gates.
    let flat = crate::geometry::curve::CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 50.0,
        start_angle: 0.0,
        end_angle: 0.01,
    }
    .to_nurbs();
    let range = flat.parameter_range();
    assert!(recognize_circular(&flat, range, 1e-2).is_none());
}
