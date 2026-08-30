use super::*;
use popkern_core::geometry::curve::CurveEval;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn radial_error(curve: &NurbsCurve, cx: f64, cy: f64, r: f64) -> f64 {
    let (t0, t1) = curve.parameter_range();
    (0..=400)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / 400.0;
            let p = curve.point_at(t);
            (((p.x - cx).powi(2) + (p.y - cy).powi(2)).sqrt() - r).abs()
        })
        .fold(0.0, f64::max)
}

#[test]
fn arc_spans_are_exactly_circular_and_angle_labelled() {
    // 250° arc → three sub-spans; every sample lies on the circle and
    // the wire's parameter is the angle at each span joint.
    let a0 = -0.7;
    let a1 = a0 + 250.0_f64.to_radians();
    let spans = arc_spans(1.0, -2.0, 3.0, a0, a1, 0.5);
    assert_eq!(spans.len(), 3);
    let wire = assemble_quadratic(&spans);
    assert!(radial_error(&wire, 1.0, -2.0, 3.0) < 1e-12);
    assert_eq!(wire.parameter_range(), (a0, a1));
    for s in &spans {
        for (t, p) in [(s.t0, s.p0), (s.t1, s.p2)] {
            let q = wire.point_at(t);
            assert!(q.distance_to(&p) < 1e-12, "knot {t} must land on {p:?}");
        }
    }
    assert!((wire.point_at(a0).x - (1.0 + 3.0 * a0.cos())).abs() < 1e-12);
}

#[test]
fn quarter_arc_is_one_span_and_full_turn_is_four() {
    assert_eq!(arc_spans(0.0, 0.0, 1.0, 0.0, FRAC_PI_2, 0.0).len(), 1);
    assert_eq!(arc_spans(0.0, 0.0, 1.0, 0.0, PI, 0.0).len(), 2);
    assert_eq!(arc_spans(0.0, 0.0, 1.0, 0.0, TAU, 0.0).len(), 4);
}

#[test]
fn mixed_line_and_arc_spans_assemble_with_full_multiplicity_joints() {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let mut spans = vec![line_span(a, b, 0.0, 1.0)];
    spans.extend(
        arc_spans(0.5, 0.0, 0.5, 0.0, PI, 0.0)
            .into_iter()
            .map(|mut s| {
                s.t0 += 1.0;
                s.t1 += 1.0;
                s
            }),
    );
    let wire = assemble_quadratic(&spans);
    assert_eq!(wire.degree, 2);
    assert_eq!(wire.control_points.len(), 2 * spans.len() + 1);
    // Interior knots doubled: 3 + 2*(n-1) + 3.
    assert_eq!(wire.knots.len(), 2 * spans.len() + 4);
    let (t0, t1) = wire.parameter_range();
    assert!(wire.point_at(t0).distance_to(&a) < 1e-12);
    assert!(
        wire.point_at(t1).distance_to(&a) < 1e-12,
        "closed back to the start"
    );
}

#[test]
fn circle_and_polygon_wires_are_closed_and_exact() {
    let c = circle_wire(2.0, 3.0, 0.75, 1.0);
    assert!(radial_error(&c, 2.0, 3.0, 0.75) < 1e-12);
    let p = polygon_wire(&[(0.0, 0.0), (2.0, 0.0), (1.0, 1.0)], 0.0);
    assert_eq!(p.degree, 1);
    assert!(p.point_at(0.0).distance_to(&p.point_at(1.0)) < 1e-12);
}
