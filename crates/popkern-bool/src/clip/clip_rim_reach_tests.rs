use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{PortableMath, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn segment(p: (f64, f64), q: (f64, f64)) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(p.0, p.1, 0.0), Point3::new(q.0, q.1, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// The distance from `(u, v)` to the segment `p`–`q`.
fn to_segment(u: f64, v: f64, p: (f64, f64), q: (f64, f64)) -> f64 {
    let (dx, dy) = (q.0 - p.0, q.1 - p.1);
    let len2 = dx * dx + dy * dy;
    let s = if len2 > 0.0 {
        (((u - p.0) * dx + (v - p.1) * dy) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (u - p.0 - s * dx).phypot(v - p.1 - s * dy)
}

/// One rim of each form with a closed-form box, and the form without.
fn rims() -> Vec<Pcurve> {
    vec![
        Pcurve::linear((-1.0, -0.5), (1.5, 0.8), 0.0, 1.0),
        Pcurve::elliptical((0.2, -0.1), 1.4, 0.6, 0.5, 0.3, 4.0, 0.0, 1.0),
        Pcurve::elliptical(
            (0.0, 0.0),
            1.2,
            1.2,
            0.0,
            0.0,
            std::f64::consts::TAU,
            0.0,
            1.0,
        ),
        Pcurve::sinusoidal(-0.5, 0.2, 0.4, -0.3, 0.0, 2.5),
        Pcurve::from_nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(0.0, 1.5, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 0.6, 1.0],
        }),
    ]
}

/// Whatever rim is cleared of a segment lies farther than the pad from
/// it at every one of 4000 points of the rim; each rim is both cleared
/// of some segments and kept for others, the full circle among them
/// for segments inside it.
#[test]
fn a_cleared_rim_never_comes_within_the_pad() {
    let pad = 0.05;
    for rim in rims() {
        let (t0, t1) = rim.parameter_range();
        let (mut cleared, mut kept) = (0, 0);
        for i in 0..14 {
            for j in 0..14 {
                let c = (-3.5 + 0.5 * f64::from(i), -3.5 + 0.5 * f64::from(j));
                for (du, dv) in [(0.3, 0.1), (0.05, -0.4), (1.2, 0.9)] {
                    let (p, q) = ((c.0 - du, c.1 - dv), (c.0 + du, c.1 + dv));
                    if !rim_clear_of_query(&rim, &segment(p, q), pad) {
                        kept += 1;
                        continue;
                    }
                    cleared += 1;
                    for k in 0..=4000 {
                        let r = rim.point_at(t0 + (t1 - t0) * f64::from(k) / 4000.0);
                        let d = to_segment(r.x, r.y, p, q);
                        assert!(d > pad, "{rim:?}: {p:?}-{q:?} at {r:?}, {d}");
                    }
                }
            }
        }
        assert!(
            cleared > 50 && kept > 20,
            "{rim:?}: cleared {cleared}, kept {kept}"
        );
    }
    // Inside the full circle, by its disc: a segment wholly inside it
    // is cleared though their boxes overlap, one within the pad is not.
    let circle = &rims()[2];
    assert!(rim_clear_of_query(
        circle,
        &segment((-0.5, 0.0), (0.5, 0.3)),
        pad
    ));
    assert!(!rim_clear_of_query(
        circle,
        &segment((-0.5, 0.0), (1.2 - 0.5 * pad, 0.0)),
        pad
    ));
}

/// A form with no closed-form box here, and a record of a weight not
/// positive, are never cleared.
#[test]
fn only_known_boxes_and_positive_weights_clear() {
    let far = segment((50.0, 50.0), (51.0, 50.0));
    let arc = Pcurve {
        geometry: PcurveGeometry::CylinderArcViaRationalQuad {
            theta_a: 0.0,
            theta_b: 1.0,
            v_const: 0.0,
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    assert!(!rim_clear_of_query(&arc, &far, 1e-4));
    let mut negative = far.clone();
    negative.weights[1] = -1.0;
    assert!(!rim_clear_of_query(&rims()[0], &negative, 1e-4));
    assert!(rim_clear_of_query(&rims()[0], &far, 1e-4));
}

/// A plate with a round window: on its top face, a record inside the
/// face crosses neither the rim nor the hole, one through the hole
/// crosses the hole, and one across the rim crosses it twice; a chord
/// cutting the hole's rim 1e-5 below its top still crosses it.
#[test]
fn the_crossing_searches_keep_their_verdicts() {
    let tol = Tolerance::default();
    let plate = BRep::prism(
        &[
            Point3::new(-10.0, -10.0, 0.0),
            Point3::new(10.0, -10.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(-10.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        4.0,
        tol,
    )
    .expect("plate");
    let window = BRep::cylinder(3.0, 6.0, tol)
        .expect("window")
        .transform(&Transform3::translation(4.0, 4.0, -1.0))
        .expect("through");
    let plate =
        crate::boolean(crate::BoolOperation::Subtraction, plate, window, tol, None).expect("cut");
    let (top, origin, u_axis, v_axis) = plate
        .faces
        .iter()
        .find_map(|(h, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::surface::SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            }) if (origin.z - 4.0).abs() < 1e-9
                && normal.z.abs() > 0.5
                && !f.inner_loops.is_empty() =>
            {
                Some((h, *origin, *u_axis, *v_axis))
            }
            _ => None,
        })
        .expect("the holed top");
    let uv = |x: f64, y: f64| {
        let d = Point3::new(x, y, 4.0) - origin;
        (d.dot(&u_axis), d.dot(&v_axis))
    };
    let record = |a: (f64, f64), b: (f64, f64)| segment(uv(a.0, a.1), uv(b.0, b.1));
    let crossings = |c: &NurbsCurve| {
        super::super::clip_uv_crossings::find_uv_crossings(&plate, top, c, &tol)
            .expect("pcurves")
            .len()
    };
    let holed = |c: &NurbsCurve| {
        super::super::clip_inner_crossing::crosses_inner_loop(&plate, top, c, &tol)
    };
    let inside = record((-6.0, -6.0), (-2.0, -5.0));
    assert_eq!((crossings(&inside), holed(&inside)), (0, false));
    let through_hole = record((0.0, 4.0), (8.0, 4.0));
    assert_eq!((crossings(&through_hole), holed(&through_hole)), (0, true));
    let across_rim = record((-12.0, 0.0), (12.0, 0.0));
    assert_eq!(crossings(&across_rim), 2);
    let grazing = record((4.0 - 2.0, 7.0 - 1e-5), (4.0 + 2.0, 7.0 - 1e-5));
    assert!(holed(&grazing));
}
