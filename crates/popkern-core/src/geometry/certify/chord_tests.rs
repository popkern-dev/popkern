use super::*;
use crate::geometry::Vec3;
use crate::geometry::certify::surface_jet::SurfaceJetSource;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

fn cubic() -> CurveGeometry {
    CurveGeometry::Nurbs(NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(2.0, -1.0, 0.5),
            Point3::new(3.0, 1.5, 0.0),
            Point3::new(4.0, 0.0, 1.0),
        ],
        weights: vec![1.0; 5],
    })
}

/// A quarter of the unit circle about the origin, rational quadratic.
fn quarter() -> CurveGeometry {
    let h = std::f64::consts::FRAC_1_SQRT_2;
    CurveGeometry::Nurbs(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        weights: vec![1.0, h, 1.0],
    })
}

/// The largest sampled distance from the curve over `[a, b]` to `p0 p1`.
fn sampled(curve: &CurveGeometry, a: f64, b: f64, p0: Point3, p1: Point3) -> f64 {
    (0..=2000)
        .map(|k| {
            let t = a + (b - a) * k as f64 / 2000.0;
            dist_to_segment(curve.point_at(t), p0, p1)
        })
        .fold(0.0, f64::max)
}

#[test]
fn curve_bounds_cover_the_sampled_deviation_and_tighten_to_the_target() {
    for (curve, a, b) in [
        (cubic(), 0.1, 0.4),
        (cubic(), 0.3, 0.8),
        (quarter(), 0.0, 1.0),
    ] {
        let (p0, p1) = (curve.point_at(a), curve.point_at(b));
        let truth = sampled(&curve, a, b, p0, p1);
        let loose = curve_chord_bound(&curve, a, b, p0, p1, f64::INFINITY).expect("bound");
        let tight = curve_chord_bound(&curve, a, b, p0, p1, 1.01 * truth).expect("bound");
        assert!(
            loose >= truth && tight >= truth,
            "{truth} vs {loose}, {tight}"
        );
        assert!(tight <= loose, "{tight} > {loose}");
        assert!(tight <= 1.05 * truth + 1e-12, "tight {tight} vs {truth}");
    }
    // The rational quarter's sagitta is exact: 1 − cos(π/4).
    let q = quarter();
    let (p0, p1) = (q.point_at(0.0), q.point_at(1.0));
    let sagitta = 1.0 - std::f64::consts::FRAC_1_SQRT_2;
    let b = curve_chord_bound(&q, 0.0, 1.0, p0, p1, 1.001 * sagitta).expect("bound");
    assert!(b >= sagitta && b <= 1.01 * sagitta, "{b} vs {sagitta}");
}

#[test]
fn analytic_curves_bound_their_chords() {
    let circle = CurveGeometry::Circle {
        center: Point3::new(1.0, 2.0, 3.0),
        normal: Vec3::new(0.0, 0.6, 0.8),
        radius: 2.5,
    };
    let ellipse = CurveGeometry::Ellipse {
        center: Point3::new(0.0, 0.0, 0.0),
        u_basis: Vec3::new(3.0, 0.0, 0.0),
        v_basis: Vec3::new(0.0, 1.0, 0.0),
        start_angle: 0.0,
        end_angle: 6.0,
    };
    let parabola = CurveGeometry::Parabola {
        vertex: Point3::new(0.0, 0.0, 0.0),
        linear: Vec3::new(1.0, 0.0, 0.0),
        quadratic: Vec3::new(0.0, 0.5, 0.0),
        t_start: -2.0,
        t_end: 2.0,
    };
    let hyperbola = CurveGeometry::Hyperbola {
        center: Point3::new(0.0, 0.0, 0.0),
        semi_major: Vec3::new(2.0, 0.0, 0.0),
        semi_minor: Vec3::new(0.0, 1.0, 0.0),
        t_start: -1.5,
        t_end: 1.5,
    };
    for (curve, a, b) in [
        (&circle, 0.3, 1.1),
        (&ellipse, 1.0, 1.6),
        (&parabola, -0.4, 0.7),
        (&hyperbola, 0.6, 1.2),
    ] {
        let (p0, p1) = (curve.point_at(a), curve.point_at(b));
        let truth = sampled(curve, a, b, p0, p1);
        let bound = curve_chord_bound(curve, a, b, p0, p1, 0.0).expect("bound");
        assert!(bound >= truth, "{curve:?}: {bound} < {truth}");
        assert!(bound <= 4.0 * truth + 1e-9, "{curve:?}: {bound} vs {truth}");
    }
    // Moved endpoints add their distance from the curve.
    let (p0, p1) = (circle.point_at(0.3), circle.point_at(1.1));
    let off = p1 + Vec3::new(0.0, 0.0, 0.1);
    let truth = sampled(&circle, 0.3, 1.1, p0, off);
    let bound = curve_chord_bound(&circle, 0.3, 1.1, p0, off, 0.0).expect("bound");
    assert!(bound >= truth && bound >= 0.1, "{bound} vs {truth}");
}

#[test]
fn flattening_meets_the_tolerance_on_every_segment() {
    let curve = cubic();
    let segments = certified_flatten(&curve, 0.0, 1.0, 1e-3, 12).expect("flattened");
    assert!(segments.len() > 4);
    assert_eq!(segments[0].t0, 0.0);
    assert_eq!(segments.last().map(|s| s.t1), Some(1.0));
    for pair in segments.windows(2) {
        assert_eq!(pair[0].t1, pair[1].t0);
    }
    for s in &segments {
        assert!(s.bound <= 1e-3);
        let (p0, p1) = (curve.point_at(s.t0), curve.point_at(s.t1));
        assert!(sampled(&curve, s.t0, s.t1, p0, p1) <= s.bound);
    }
    assert!(matches!(
        certified_flatten(&curve, 0.0, 1.0, 1e-9, 2),
        Err(FlattenError::DepthLimit { .. })
    ));
}

/// The largest sampled `|L(x) − S(x)|` over the triangle.
fn sampled_triangle(s: &SurfaceGeometry, w: [(f64, f64); 3], p: [Point3; 3]) -> f64 {
    let n = 60;
    let mut worst = 0.0_f64;
    for i in 0..=n {
        for j in 0..=n - i {
            let (l1, l2) = (i as f64 / n as f64, j as f64 / n as f64);
            let l0 = 1.0 - l1 - l2;
            let u = l0 * w[0].0 + l1 * w[1].0 + l2 * w[2].0;
            let v = l0 * w[0].1 + l1 * w[1].1 + l2 * w[2].1;
            let l = Point3::new(
                l0 * p[0].x + l1 * p[1].x + l2 * p[2].x,
                l0 * p[0].y + l1 * p[1].y + l2 * p[2].y,
                l0 * p[0].z + l1 * p[1].z + l2 * p[2].z,
            );
            worst = worst.max(l.distance_to(&s.point_at(u, v)));
        }
    }
    worst
}

/// A bicubic bump: flat at the corners, peaked in the middle.
fn bump() -> SurfaceGeometry {
    let mut control_points = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            let z = if (1..3).contains(&i) && (1..3).contains(&j) {
                2.0
            } else {
                0.0
            };
            control_points.push(Point3::new(i as f64, j as f64, z));
        }
    }
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: knots.clone(),
        knots_v: knots,
        control_points,
        num_u: 4,
        num_v: 4,
        weights: vec![1.0; 16],
    })
}

fn surfaces() -> Vec<(SurfaceGeometry, [(f64, f64); 3])> {
    let axis = Vec3::new(0.3, -0.4, 0.866).normalized().expect("axis");
    vec![
        (
            SurfaceGeometry::Cylinder {
                origin: Point3::new(0.0, 0.0, 0.0),
                axis: Vec3::Z,
                radius: 2.0,
            },
            [(0.1, 0.0), (0.6, 0.2), (0.3, 1.0)],
        ),
        (
            SurfaceGeometry::Sphere {
                center: Point3::new(1.0, 0.0, 0.0),
                radius: 1.5,
                axis,
            },
            [(0.5, 0.8), (1.0, 0.9), (0.7, 1.3)],
        ),
        (
            SurfaceGeometry::Torus {
                center: Point3::new(0.0, 0.0, 0.0),
                axis: Vec3::Z,
                major_radius: 3.0,
                minor_radius: 1.0,
            },
            [(0.2, 0.1), (0.5, 0.6), (0.1, 0.9)],
        ),
        (bump(), [(0.1, 0.1), (0.9, 0.2), (0.4, 0.9)]),
    ]
}

#[test]
fn triangle_bounds_cover_the_sampled_deviation() {
    for (s, w) in surfaces() {
        let p = w.map(|(u, v)| s.point_at(u, v));
        let mut source = SurfaceJetSource::new(&s).expect("source");
        let truth = sampled_triangle(&s, w, p);
        let none = C0Breaks::default();
        let loose = triangle_chord_bound(&mut source, &none, w, p, f64::INFINITY).expect("loose");
        let tight = triangle_chord_bound(&mut source, &none, w, p, 1.2 * truth).expect("tight");
        assert!(truth > 1e-3, "{s:?}: a flat fixture proves nothing");
        assert!(
            loose >= truth && tight >= truth,
            "{s:?}: {truth} vs {loose}, {tight}"
        );
        assert!(
            tight <= loose && tight <= 3.0 * truth,
            "{s:?}: {tight} vs {truth}"
        );
    }
}

/// A cubic B-spline with one raised control point in the middle of a
/// 9 × 9 net: flat outside `[1/6, 5/6]²`, a bump inside.
fn local_bump() -> SurfaceGeometry {
    let mut knots = vec![0.0; 4];
    knots.extend((1..6).map(|k| k as f64 / 6.0));
    knots.extend([1.0; 4]);
    // Greville abscissae: x(u) = u and y(v) = v exactly, so the flat
    // part has no second derivative at all.
    let greville = |i: usize| (knots[i + 1] + knots[i + 2] + knots[i + 3]) / 3.0;
    let mut control_points = Vec::new();
    for i in 0..9 {
        for j in 0..9 {
            let z = if (i, j) == (4, 4) { 2.0 } else { 0.0 };
            control_points.push(Point3::new(greville(i), greville(j), z));
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: knots.clone(),
        knots_v: knots,
        control_points,
        num_u: 9,
        num_v: 9,
        weights: vec![1.0; 81],
    })
}

#[test]
fn a_bound_from_vertex_curvature_alone_undercuts_the_bump() {
    // Negative control: M taken from the second partials at the three
    // vertices (sampling) misses a bump inside the triangle whose
    // vertices sit where the surface is flat; the certified bound does
    // not.
    let s = local_bump();
    let w = [(0.05, 0.05), (0.98, 0.05), (0.05, 0.98)];
    let p = w.map(|(u, v)| s.point_at(u, v));
    let truth = sampled_triangle(&s, w, p);
    assert!(truth > 0.1, "the bump is inside: {truth}");
    let mut source = SurfaceJetSource::new(&s).expect("source");
    let mut sampled_m = 0.0_f64;
    for (u, v) in w {
        let jet = source
            .jet_piecewise(Iv::point(u), Iv::point(v))
            .expect("jet")
            .0;
        let [a, b, c] = [(2, 0), (1, 1), (0, 2)].map(|(i, j)| jet.partial(i, j).mid().length());
        sampled_m = sampled_m.max(0.5 * (a + c) + (0.25 * (a - c).powi(2) + b * b).sqrt());
    }
    let sampled_bound = 0.5 * sampled_m * 0.5;
    assert!(sampled_bound < truth, "{sampled_bound} vs {truth}");
    let certified =
        triangle_chord_bound(&mut source, &C0Breaks::default(), w, p, 0.0).expect("bound");
    assert!(certified >= truth, "{certified} vs {truth}");
}

#[test]
fn a_crease_line_inside_a_triangle_is_cut_along() {
    // Degree 2 in u with a double knot at 0.5: C⁰ there (a crease).
    let mut control_points = Vec::new();
    for i in 0..5 {
        for j in 0..3 {
            let z = if i == 2 { 1.0 } else { 0.0 };
            control_points.push(Point3::new(i as f64 * 0.5, j as f64, z));
        }
    }
    let s = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points,
        num_u: 5,
        num_v: 3,
        weights: vec![1.0; 15],
    });
    let w = [(0.2, 0.1), (0.8, 0.2), (0.45, 0.9)];
    let p = w.map(|(u, v)| s.point_at(u, v));
    let truth = sampled_triangle(&s, w, p);
    let mut source = SurfaceJetSource::new(&s).expect("source");
    let breaks = C0Breaks {
        u: vec![0.5],
        ..C0Breaks::default()
    };
    let bound = triangle_chord_bound(&mut source, &breaks, w, p, 0.0).expect("bound");
    assert!(bound >= truth, "{bound} vs {truth}");
    assert_eq!(
        chord_triangle::split(w, 0, 0.5).len(),
        3,
        "a triangle and a quad"
    );
}

#[test]
fn the_sphere_jet_follows_a_tilted_pole() {
    let s = SurfaceGeometry::Sphere {
        center: Point3::new(0.5, -1.0, 2.0),
        radius: 1.25,
        axis: Vec3::new(1.0, 1.0, 0.0).normalized().expect("axis"),
    };
    let mut source = SurfaceJetSource::new(&s).expect("source");
    for (u, v) in [(0.3, 0.4), (2.0, 1.2), (5.0, 2.9)] {
        let jet = source.jet(Iv::point(u), Iv::point(v)).expect("jet");
        let p = s.point_at(u, v);
        assert!(point_to_box(p, jet.partial(0, 0)) < 1e-12, "({u}, {v})");
        let du = s.partial_u(u, v);
        assert!(
            (jet.partial(1, 0).mid() - du).length() < 1e-12,
            "S_u at ({u}, {v})"
        );
    }
}
