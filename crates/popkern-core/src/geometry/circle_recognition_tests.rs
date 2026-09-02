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
