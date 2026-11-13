use super::*;
use popkern_core::geometry::{Point3, Vec3};

fn sphere() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
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

const V0: f64 = 2.0 * std::f64::consts::FRAC_PI_3;

/// The rim latitude of a sphere of radius 2's cap, as one boundary edge.
fn rim() -> Vec<(HalfEdgeHandle, NurbsCurve)> {
    vec![(HalfEdgeHandle::default(), nurbs(&[(-0.5, V0), (0.5, V0)]))]
}

/// A closed record of the chart: an arc rising to `delta` short of the
/// rim at t = 1/4 (`a` below it at its ends, t = 0 and 1/2), then a
/// straight return.
fn arc(delta: f64, a: f64) -> NurbsCurve {
    let low = V0 - delta - a;
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        control_points: [
            (-0.1, low),
            (0.0, V0 - delta + a),
            (0.1, low),
            (0.0, low),
            (-0.1, low),
        ]
        .iter()
        .map(|&(u, v)| Point3::new(u, v, 0.0))
        .collect(),
        weights: vec![1.0; 5],
    }
}

/// Reports at arc parameters `s` (the arc's own, 0..1): t = s / 2.
fn on_arc(s: f64) -> UvCrossing {
    at(0.5 * s)
}

fn at(t: f64) -> UvCrossing {
    UvCrossing {
        t_curve: t,
        boundary_he: HalfEdgeHandle::default(),
        t_boundary: 0.5,
    }
}

/// A closed record passing 4e-7 inside the rim in space (2e-7 rad), its
/// stretch within the reach reported as its only two crossings, grazes
/// it: both reports go. Rising 2e-7 rad past the rim between two true crossings it is
/// still one place with it there (a sub-tolerance excursion), and goes
/// too.
#[test]
fn a_pass_within_the_reach_is_no_crossing() {
    let tol = Tolerance::default();
    let s = sphere();
    for delta in [2e-7f64, -2e-7] {
        let half = 0.5 * ((5e-7 - delta) / 1e-3).sqrt();
        let crossings = vec![on_arc(0.5 - half), on_arc(0.5 + half)];
        let kept = drop_grazes((Some(&s), &rim()), &arc(delta, 1e-3), crossings, &tol);
        assert!(kept.is_empty(), "delta {delta}: {} kept", kept.len());
    }
}

/// A line crossing the rim shallowly, its stretch within the reach
/// reported at both ends, leaves to the other side: a crossing, its
/// reports kept (the run resolution makes them one). A pair between
/// which the record leaves the rim by more than the reach is two
/// crossings, kept; an open record's graze, or one among other
/// crossings, keeps its reports for the clip's arms; with no surface
/// nothing is judged.
#[test]
fn a_crossing_or_a_parting_keeps_its_reports() {
    let tol = Tolerance::default();
    let s = sphere();
    let line = nurbs(&[(-0.1, V0 - 1e-4), (0.1, V0 + 1e-4)]);
    let crossings = vec![at(0.5 - 2.4e-3), at(0.5 + 2.4e-3)];
    let kept = drop_grazes((Some(&s), &rim()), &line, crossings, &tol);
    assert_eq!(kept.len(), 2);
    let crossings = vec![on_arc(0.4), on_arc(0.6)];
    assert_eq!(
        drop_grazes((Some(&s), &rim()), &arc(1e-4, 1e-3), crossings, &tol).len(),
        2
    );
    let half = 0.5 * (3e-7f64 / 1e-3).sqrt();
    let open = nurbs(&[
        (-0.1, V0 - 2e-7 - 1e-3),
        (0.0, V0 - 2e-7 + 1e-3),
        (0.1, V0 - 2e-7 - 1e-3),
    ]);
    let crossings = vec![at(0.5 - half), at(0.5 + half)];
    assert_eq!(
        drop_grazes((Some(&s), &rim()), &open, crossings, &tol).len(),
        2
    );
    let crossings = vec![on_arc(0.1), on_arc(0.5 - half), on_arc(0.5 + half)];
    assert_eq!(
        drop_grazes((Some(&s), &rim()), &arc(2e-7, 1e-3), crossings, &tol).len(),
        3
    );
    let crossings = vec![on_arc(0.5 - half), on_arc(0.5 + half)];
    assert_eq!(
        drop_grazes((None, &rim()), &arc(2e-7, 1e-3), crossings, &tol).len(),
        2
    );
}
