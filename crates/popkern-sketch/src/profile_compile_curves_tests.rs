use super::*;
use crate::plane::SketchPlane;
use crate::profile::{AnalyticProfile, PathSegment, ProfileShape};
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;

fn samples(curve: &NurbsCurve, n: usize) -> Vec<Point3> {
    let (t0, t1) = curve.parameter_range();
    (0..=n)
        .map(|i| curve.point_at(t0 + (t1 - t0) * i as f64 / n as f64))
        .collect()
}

/// The area a closed wire encloses (shoelace over dense samples).
fn enclosed_area(curve: &NurbsCurve) -> f64 {
    let pts = samples(curve, 20_000);
    pts.windows(2)
        .map(|w| w[0].x * w[1].y - w[1].x * w[0].y)
        .sum::<f64>()
        .abs()
        / 2.0
}

/// The worst distance of `curve`'s points from the level set `f = 0`
/// of an implicit equation, over its whole range.
fn worst_residual(curve: &NurbsCurve, f: impl Fn(f64, f64) -> f64) -> f64 {
    samples(curve, 2000)
        .iter()
        .map(|p| f(p.x, p.y).abs())
        .fold(0.0, f64::max)
}

/// A profile's face: validates, every edge exact on `f = 0` (or on a
/// line), and extruded to a body that validates.
fn face_serves(
    profile: &AnalyticProfile<f64>,
    on_curve: impl Fn(f64, f64) -> f64,
    line_edges: usize,
) {
    let result = profile
        .to_brep(SketchPlane::xy(), Tolerance::default(), None)
        .expect("face compiles");
    let mut brep = result.brep;
    brep.validate().expect("lamina validates");
    // A lamina carries each curve on both of its faces' loops: count
    // each place once, by its midpoint.
    let mut seen: Vec<Point3> = Vec::new();
    let mut lines = 0;
    for e in brep.edges.values() {
        let c = e.curve.as_ref().expect("edge curve");
        let worst = (0..=50)
            .map(|i| c.t_start + (c.t_end - c.t_start) * i as f64 / 50.0)
            .map(|t| {
                let p = c.geometry.point_at(t);
                on_curve(p.x, p.y).abs()
            })
            .fold(0.0, f64::max);
        let mid = c.geometry.point_at((c.t_start + c.t_end) / 2.0);
        if seen.iter().any(|q| q.distance_to(&mid) < 1e-9) {
            continue;
        }
        seen.push(mid);
        if worst > 1e-9 {
            lines += 1;
        }
    }
    assert_eq!(lines, line_edges, "edges off the curve");
    let face = result.regions[0].face;
    brep.extrude(face, Vec3::new(0.0, 0.0, 1.0), 2.0, None)
        .expect("extrudes");
    brep.validate().expect("extruded body validates");
}

#[test]
fn an_ellipse_is_exact_in_both_outputs() {
    let profile = AnalyticProfile {
        center: [1.0, -2.0],
        ..AnalyticProfile::new(ProfileShape::Ellipse { rx: 5.0, ry: 3.0 })
    };
    let on = |x: f64, y: f64| ((x - 1.0) / 5.0).powi(2) + ((y + 2.0) / 3.0).powi(2) - 1.0;
    let wire = profile.placed_wire(&SketchPlane::xy()).expect("wire");
    assert_eq!(wire.degree, 2);
    assert!(worst_residual(&wire, on) < 1e-12);
    let area = enclosed_area(&wire);
    assert!((area - 15.0 * std::f64::consts::PI).abs() < 1e-5, "{area}");
    face_serves(&profile, on, 0);
}

#[test]
fn a_conic_of_weight_one_is_a_parabola_and_closes_with_its_chord() {
    // y = x² from (-1, 1) to (1, 1): shoulder (0, -1).
    let profile = AnalyticProfile::new(ProfileShape::Path {
        start: [-1.0, 1.0],
        segments: vec![PathSegment::Conic {
            control: [0.0, -1.0],
            to: [1.0, 1.0],
            weight: 1.0,
        }],
    });
    let wire = profile.local_wire(&SketchPlane::xy()).expect("wire");
    // The region between y = x² and y = 1 has area 4/3.
    assert!((enclosed_area(&wire) - 4.0 / 3.0).abs() < 1e-6);
    face_serves(&profile, |x, y| y - x * x, 1);
}

#[test]
fn a_conic_above_weight_one_is_a_hyperbola() {
    // The arc of x y = 1 from (1/2, 2) to (2, 1/2): its end tangents
    // y = -4x + 4 and y = -x/4 + 1 meet at (4/5, 4/5), and the weight
    // that puts the curve through the vertex (1, 1) is the hyperbola's.
    let (a, b, m) = ((0.5, 2.0), (2.0, 0.5), (0.8, 0.8));
    // For a conic through q at t = 1/2: q = (a + 2 w m + b) / (2 + 2 w).
    let w = ((a.0 + b.0) / 2.0 - 1.0) / (1.0 - m.0);
    assert!(w > 1.0);
    let Piece::Quadratic(span) = conic_piece(a, m, b, w).expect("conic") else {
        panic!("one span");
    };
    let curve = assemble_quadratic(&[span]);
    assert!(worst_residual(&curve, |x, y| x * y - 1.0) < 1e-12);
}

#[test]
fn a_cubic_segment_joins_the_lines_and_arcs_round_it() {
    // A slot: a cubic across the top, a line down, a half circle
    // round the bottom and a line back up.
    let profile = AnalyticProfile::new(ProfileShape::Path {
        start: [-2.0, 0.0],
        segments: vec![
            PathSegment::Spline {
                control: vec![[-1.0, 2.0], [1.0, -1.0]],
                to: [2.0, 0.0],
                degree: 3,
                weights: vec![],
                knots: vec![],
            },
            PathSegment::Line { to: [2.0, -2.0] },
            PathSegment::Arc {
                via: [0.0, -4.0],
                to: [-2.0, -2.0],
            },
        ],
    });
    let wire = profile.local_wire(&SketchPlane::xy()).expect("wire");
    assert_eq!(wire.degree, 3);
    let pts = samples(&wire, 4);
    assert!(pts[0].distance_to(&pts[4]) < 1e-12, "closed");
    // The half disc below, the 4 × 2 box, and the cubic's signed
    // area over its chord y = 0 (by quadrature of the cubic alone).
    let cubic = {
        let c = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: [(-2.0, 0.0), (-1.0, 2.0), (1.0, -1.0), (2.0, 0.0)]
                .into_iter()
                .map(point)
                .collect(),
            weights: vec![1.0; 4],
        };
        // ∫ y dx along the cubic: above the chord less below it.
        let pts = samples(&c, 20_000);
        pts.windows(2)
            .map(|w| (w[1].x - w[0].x) * (w[0].y + w[1].y) / 2.0)
            .sum::<f64>()
    };
    let exact = 2.0 * std::f64::consts::PI + 8.0 + cubic;
    let got = enclosed_area(&wire);
    assert!((got - exact).abs() < 1e-6, "{got} vs {exact}");
    // Every wire sample lies on one of the pieces: the cubic's points
    // are its own (checked by the face path's edge curves).
    face_serves(
        &profile,
        |x, y| {
            let on_arc = (x * x + (y + 2.0).powi(2)).sqrt() - 2.0;
            on_arc.abs().min((x.abs() - 2.0).abs())
        },
        1,
    );
}

#[test]
fn a_closed_spline_clamps_to_the_same_curve_smooth_at_its_seam() {
    let points = vec![(0.0, 0.0), (4.0, -1.0), (6.0, 2.0), (3.0, 5.0), (-1.0, 3.0)];
    let weights = vec![1.0, 2.0, 1.0, 0.5, 1.0];
    let open = periodic_unclamped(&points, 3, &weights);
    let wire = periodic_wire(&points, 3, &weights).expect("wire");
    assert_eq!(wire.knots[..4], [0.0; 4]);
    for i in 0..=200 {
        let t = 5.0 * i as f64 / 200.0;
        let d = open.point_at(t).distance_to(&wire.point_at(t));
        assert!(d < 1e-12, "t {t}: {d}");
    }
    // C2 across the seam of the periodic form.
    let (start, end) = (open.derivatives(0.0, 2), open.derivatives(5.0, 2));
    for k in 0..3 {
        assert!((start[k] - end[k]).length() < 1e-9, "derivative {k}");
    }
    let profile = AnalyticProfile::new(ProfileShape::Spline {
        points: points.iter().map(|p| [p.0, p.1]).collect(),
        degree: 3,
        weights,
    });
    let result = profile
        .to_brep(SketchPlane::xy(), Tolerance::default(), None)
        .expect("face");
    let mut brep = result.brep;
    brep.validate().expect("lamina validates");
    brep.extrude(result.regions[0].face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .expect("extrudes");
    brep.validate().expect("prism validates");
}

#[test]
fn malformed_curves_are_named() {
    let named = |r: Result<Piece, ProfileError>, needle: &str| match r {
        Err(ProfileError::InvalidSpline { reason, .. }) => {
            assert!(reason.contains(needle), "{reason}")
        }
        Err(e) => panic!("{e}"),
        Ok(_) => panic!("{needle}: accepted"),
    };
    named(
        conic_piece((0.0, 0.0), (1.0, 1.0), (2.0, 0.0), 0.0),
        "not positive",
    );
    let to = (3.0, 0.0);
    let ctl = [(1.0, 1.0), (2.0, 1.0)];
    named(spline_piece((0.0, 0.0), &ctl, to, 4, &[], &[]), "degree 4");
    named(
        spline_piece((0.0, 0.0), &ctl, to, 3, &[1.0], &[]),
        "1 weights",
    );
    named(
        spline_piece(
            (0.0, 0.0),
            &ctl,
            to,
            2,
            &[],
            &[0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
        ),
        "clamped",
    );
    named(
        spline_piece(
            (0.0, 0.0),
            &ctl,
            to,
            3,
            &[],
            &[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        ),
        "knots for",
    );
    named(
        periodic_piece(vec![(0.0, 0.0), (1.0, 0.0)], 1, &[]),
        "at least 3",
    );
}
