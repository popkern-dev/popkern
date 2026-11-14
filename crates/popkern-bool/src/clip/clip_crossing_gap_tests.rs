use super::*;
use popkern_core::geometry::{Point3, Vec3};

fn sphere(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius,
        axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn nurbs(points: &[(f64, f64)]) -> NurbsCurve {
    let degree = points.len() - 1;
    let mut knots = vec![0.0; degree + 1];
    knots.extend(vec![1.0; degree + 1]);
    NurbsCurve {
        degree,
        knots,
        control_points: points
            .iter()
            .map(|&(u, v)| Point3::new(u, v, 0.0))
            .collect(),
        weights: vec![1.0; points.len()],
    }
}

/// The rim latitude v0 of a sphere's cap, and an arc of the chart rising
/// to `delta` short of it at its apex (t = 1/2), as the circle two balls
/// overlap in rises to a ball's rim on a box's top (PK-1180).
fn rim_and_arc(delta: f64, lift: f64) -> (NurbsCurve, NurbsCurve) {
    let v0 = 2.0 * std::f64::consts::FRAC_PI_3;
    let a = 1e-3;
    (
        nurbs(&[(-0.5, v0), (0.5, v0)]),
        nurbs(&[
            (lift - 0.1, v0 - delta - a),
            (lift, v0 - delta + a),
            (lift + 0.1, v0 - delta - a),
        ]),
    )
}

/// A near miss under the tolerance in the chart (7.7e-7 rad) is 1.54e-6
/// apart on a sphere of radius 2, past the tolerance in space, however
/// near its apex the intersector reported it, and on a lift by a period;
/// a shallow crossing refines to no gap.
#[test]
fn a_near_miss_on_a_stretched_chart_is_no_crossing() {
    let tol = Tolerance::default();
    let s = sphere(2.0);
    let tau = std::f64::consts::TAU;
    for lift in [0.0, tau] {
        let (rim, arc) = rim_and_arc(7.7e-7, lift);
        for t in [0.5, 0.5 - 1e-4, 0.5 + 1e-4] {
            let gap = crossing_gap(&s, (&rim, &arc), (0.5, t), 0.0);
            // Never under the least distance; near a tangency the
            // refinement creeps along the curves.
            assert!(
                gap > 1.54e-6 - 2e-8 && gap < 1.6e-6,
                "lift {lift} t {t}: gap {gap}"
            );
            assert!(afar(Some(&s), (&rim, &arc), (0.5, t), 0.0, &tol));
        }
        // Rising 2e-7 past the rim, the arc crosses it at t = 1/2 -
        // sqrt(5e-5), u = -sqrt(5e-5) / 5 (the rim's t = 1/2 + u): a
        // report near there refines onto it; at the apex the curves run
        // parallel and the refinement stops 4e-7 apart, within reach.
        let (rim, arc) = rim_and_arc(-2e-7, lift);
        let (t, u) = (0.5 - 5e-5f64.sqrt(), -(5e-5f64.sqrt()) / 5.0);
        let gap = crossing_gap(&s, (&rim, &arc), (0.5 + u + 1e-5, t + 1e-4), 0.0);
        assert!(gap < 1e-9, "lift {lift}: crossing gap {gap}");
        assert!(!afar(Some(&s), (&rim, &arc), (0.5, 0.5), 0.0, &tol));
    }
}

/// A real crossing reported off its place by the chart's tolerance lies
/// past the tolerance in space on a wide chart (a sphere of radius 100),
/// and refines to its true place; with no surface nothing is judged.
#[test]
fn a_crossing_reported_off_its_place_refines_onto_it() {
    let tol = Tolerance::default();
    let s = sphere(100.0);
    let v0 = 1.2;
    let rim = nurbs(&[(-0.5, v0), (0.5, v0)]);
    let chord = nurbs(&[(-0.01, v0 - 0.01), (0.01, v0 + 0.01)]);
    let off = (0.5, 0.5 + 1e-6);
    let p = |c: &NurbsCurve, t: f64| {
        let q = c.point_at(t);
        popkern_core::geometry::surface::SurfaceEval::point_at(&s, q.x, q.y)
    };
    assert!(p(&rim, off.0).distance_to(&p(&chord, off.1)) > 2.0 * tol.linear());
    assert!(crossing_gap(&s, (&rim, &chord), off, 0.0) < 1e-9);
    assert!(!afar(Some(&s), (&rim, &chord), off, 0.0, &tol));
    let (rim, arc) = rim_and_arc(7.7e-7, 0.0);
    assert!(!afar(None, (&rim, &arc), (0.5, 0.5), 0.0, &tol));
}

/// Only the crossings from `from` on are judged: those before are the
/// passes that judged their own.
#[test]
fn only_the_lifted_crossings_are_judged_again() {
    let tol = Tolerance::default();
    let s = sphere(2.0);
    let (rim, arc) = rim_and_arc(7.7e-7, 0.0);
    let he = HalfEdgeHandle::default();
    let at = |t: f64| UvCrossing {
        t_curve: t,
        boundary_he: he,
        t_boundary: 0.5,
    };
    let boundary = vec![(he, rim)];
    let mut crossings = vec![at(0.5), at(0.49), at(0.51)];
    keep_one_place(Some(&s), &boundary, &arc, &mut crossings, 1, &tol);
    assert_eq!(crossings.len(), 1);
    assert_eq!(crossings[0].t_curve, 0.5);
}
