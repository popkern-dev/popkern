//! Tests for `cylinder_cylinder` (split out for the 500-line cap).
use super::*;
use crate::geometry::curve::CurveEval as _;

/// Two overlapping PARALLEL cylinders meet in exactly two lines
/// parallel to the axis; each line's samples lie on BOTH
/// surfaces and the pcurves agree with the 3D curve. The
/// marching path solved this class numerically — and on the
/// full-chart wall pair it found only ONE of the two lines, so a
/// plain two-cylinder union SERVED the tool body alone (volume
/// 0.2356 for a true union of 0.4797, measured PK-427).
#[test]
fn parallel_overlapping_cylinders_yield_two_exact_lines() {
    let a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.3,
    };
    let b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.4, 0.0, -0.1),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.25,
    };
    let tol = Tolerance::default();
    let curves = try_cylinder_cylinder(&a, (0.0, 1.0), &b, (0.0, 1.2), &tol)
        .expect("parallel pair takes the analytic arm");
    assert_eq!(curves.len(), 2, "two lines");
    for c in &curves {
        let SurfaceSurfaceIntersection::Curve(cos) = c else {
            panic!("line curve expected");
        };
        let (t0, t1) = cos.curve_3d.parameter_range();
        for k in 0..=4 {
            let t = t0 + (t1 - t0) * k as f64 / 4.0;
            let p = cos.curve_3d.point_at(t);
            let da = ((p.x * p.x + p.y * p.y).sqrt() - 0.3).abs();
            let db = (((p.x - 0.4) * (p.x - 0.4) + p.y * p.y).sqrt() - 0.25).abs();
            assert!(
                da < 1e-12 && db < 1e-12,
                "on both surfaces: {da:.2e} {db:.2e}"
            );
            use crate::geometry::SurfaceEval as _;
            let uv1 = cos.curve_uv1.point_at(t);
            let uv2 = cos.curve_uv2.point_at(t);
            assert!(a.point_at(uv1.x, uv1.y).distance_to(&p) < 1e-12);
            assert!(b.point_at(uv2.x, uv2.y).distance_to(&p) < 1e-12);
        }
    }
}

/// Separated and nested parallel pairs are exactly empty;
/// near-coincident pairs decline to the coincident detector.
#[test]
fn parallel_disjoint_nested_and_coincident_regimes() {
    let tol = Tolerance::default();
    let at = |x: f64, r: f64| SurfaceGeometry::Cylinder {
        origin: Point3::new(x, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    };
    let far = try_cylinder_cylinder(&at(0.0, 0.3), (0.0, 1.0), &at(1.0, 0.25), (0.0, 1.0), &tol)
        .expect("analytic arm");
    assert!(far.is_empty(), "separated pair is exactly empty");
    let nested = try_cylinder_cylinder(&at(0.0, 0.3), (0.0, 1.0), &at(0.01, 0.1), (0.0, 1.0), &tol)
        .expect("analytic arm");
    assert!(nested.is_empty(), "nested pair is exactly empty");
    assert!(
        try_cylinder_cylinder(&at(0.0, 0.3), (0.0, 1.0), &at(0.0, 0.3), (0.0, 1.0), &tol).is_none(),
        "coincident pair stays with the coincident detector"
    );
}

#[test]
fn cross_drill_pierce_loops_exact() {
    let body = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 16.0,
    };
    let drill = SurfaceGeometry::Cylinder {
        origin: Point3::new(-40.0, 0.0, 20.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 6.0,
    };
    let tol = Tolerance::default();
    let curves = try_cylinder_cylinder(&body, (0.0, 40.0), &drill, (0.0, 80.0), &tol)
        .expect("cylinder pair must take the analytic path");
    assert_eq!(curves.len(), 2, "full pierce must yield two loops");
    for c in &curves {
        let SurfaceSurfaceIntersection::Curve(cos) = c else {
            panic!("expected curve");
        };
        let (t0, t1) = cos.curve_3d.parameter_range();
        let mut worst = 0.0_f64;
        for i in 0..=200 {
            let t = t0 + (t1 - t0) * f64::from(i) / 200.0;
            let p = cos.curve_3d.point_at(t);
            let d_body = ((p.x * p.x + p.y * p.y).sqrt() - 16.0).abs();
            let d_drill = ((p.y * p.y + (p.z - 20.0) * (p.z - 20.0)).sqrt() - 6.0).abs();
            worst = worst.max(d_body.max(d_drill));
        }
        assert!(
            worst < 1e-6,
            "pierce loop deviates {worst:.3e} from the surfaces"
        );
    }
}

/// The worst distance of a fitted curve from two cylinders about `z`
/// (radius `ra`) and about `x` through `(0, 0, zc)` (radius `rb`).
fn worst_off(
    cos: &crate::geometry::curve_on_surface::CurveOnSurface,
    ra: f64,
    rb: f64,
    zc: f64,
) -> f64 {
    let (t0, t1) = cos.curve_3d.parameter_range();
    (0..=2000)
        .map(|i| {
            let p = cos
                .curve_3d
                .point_at(t0 + (t1 - t0) * f64::from(i) / 2000.0);
            let da = ((p.x * p.x + p.y * p.y).sqrt() - ra).abs();
            let db = ((p.y * p.y + (p.z - zc) * (p.z - zc)).sqrt() - rb).abs();
            da.max(db)
        })
        .fold(0.0, f64::max)
}

fn tee_pair(rb: f64) -> (SurfaceGeometry, SurfaceGeometry) {
    (
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 40.0,
        },
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 60.0),
            axis: Vec3::new(1.0, 0.0, 0.0),
            radius: rb,
        },
    )
}

#[test]
fn equal_radius_tee_splits_at_the_singular_crossings() {
    // PK-727: equal radii, intersecting axes — two ellipses crossing at
    // (0, ±40, 60). The scan must return four OPEN branches, each from
    // one crossing to the other, not one self-touching loop.
    let (main, branch) = tee_pair(40.0);
    let tol = Tolerance::default();
    let curves = try_cylinder_cylinder(&main, (-60.0, 180.0), &branch, (-50.0, 150.0), &tol)
        .expect("analytic path");
    assert_eq!(curves.len(), 4, "two crossing ellipses = four branches");
    for c in &curves {
        let SurfaceSurfaceIntersection::Curve(cos) = c else {
            panic!("expected curve");
        };
        assert!(!cos.is_closed, "branches between crossings are open");
        let (t0, t1) = cos.curve_3d.parameter_range();
        for end in [cos.curve_3d.point_at(t0), cos.curve_3d.point_at(t1)] {
            let at_crossing = (end.x.abs() < 1e-6)
                && ((end.y.abs() - 40.0).abs() < 1e-6)
                && ((end.z - 60.0).abs() < 1e-6);
            assert!(at_crossing, "branch end {end:?} is not a crossing point");
        }
        let off = worst_off(cos, 40.0, 40.0, 60.0);
        assert!(off < 1e-6, "branch strays {off:.3e} from the cylinders");
    }
}

#[test]
fn millimetre_pierce_loops_are_fitted_to_tolerance() {
    // PK-728: 512 fixed samples left an r = 30 × r = 40 loop 2.4e-6 off
    // its surfaces at millimetre scale; the fit is densified until it
    // lies within a tenth of the linear tolerance.
    let (main, branch) = tee_pair(30.0);
    let tol = Tolerance::default();
    let curves = try_cylinder_cylinder(&main, (-60.0, 180.0), &branch, (-100.0, 300.0), &tol)
        .expect("analytic path");
    assert_eq!(curves.len(), 2, "full pierce: two loops");
    for c in &curves {
        let SurfaceSurfaceIntersection::Curve(cos) = c else {
            panic!("expected curve");
        };
        let off = worst_off(cos, 40.0, 30.0, 60.0);
        assert!(
            off <= 0.1 * tol.linear() * 1.5,
            "loop strays {off:.3e} (budget {:.1e})",
            0.1 * tol.linear()
        );
    }
}

/// PK-782: two perpendicular r 5 cylinders whose axes are 9.53 apart
/// meet in one small loop (0.47 of overlap). Sampled uniformly in the
/// ruling angle, its tangency ends (roots moving as √(θ − θ_end))
/// never fitted to tolerance and the loop densified to the cap, 32770
/// points, each later crossing test against it taking seconds. Cosine
/// spacing makes the loop smooth in its sampling parameter: one loop,
/// a few hundred points, on both cylinders within the fit budget.
#[test]
fn a_small_near_tangent_loop_fits_in_hundreds_of_points() {
    let a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let b = SurfaceGeometry::Cylinder {
        origin: Point3::new(9.53, 0.0, 60.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        radius: 5.0,
    };
    let tol = Tolerance::default();
    let loops = try_cylinder_cylinder(&a, (0.0, 120.0), &b, (-20.0, 20.0), &tol).expect("analytic");
    assert_eq!(loops.len(), 1, "one loop");
    let SurfaceSurfaceIntersection::Curve(cos) = &loops[0] else {
        panic!("a curve");
    };
    let points = cos.curve_3d.control_points.len();
    assert!(points <= 500, "{points} control points");
    let (t0, t1) = cos.curve_3d.parameter_range();
    let worst = (0..=4000)
        .map(|i| {
            let p = cos
                .curve_3d
                .point_at(t0 + (t1 - t0) * f64::from(i) / 4000.0);
            let da = ((p.x * p.x + p.y * p.y).sqrt() - 5.0).abs();
            let db = (((p.x - 9.53).powi(2) + (p.z - 60.0).powi(2)).sqrt() - 5.0).abs();
            da.max(db)
        })
        .fold(0.0, f64::max);
    assert!(
        worst <= 0.1 * tol.linear(),
        "off the cylinders by {worst:e}"
    );
}

/// The curves of `try_cylinder_cylinder` for unit cylinders A (axis z
/// through the origin) and B (axis `axis` through `origin`): each closed
/// curve's length and its worst distance off either cylinder, sampled
/// densely.
fn unit_pair_loops(origin: Point3, axis: Vec3) -> Vec<(f64, f64)> {
    let a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let b = SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius: 1.0,
    };
    let tol = Tolerance::default();
    let curves = try_cylinder_cylinder(&a, (-2.0, 2.0), &b, (-2.0, 2.0), &tol).expect("arm");
    let an = axis.normalized().expect("axis");
    curves
        .iter()
        .map(|c| {
            let SurfaceSurfaceIntersection::Curve(cos) = c else {
                panic!("a curve");
            };
            assert!(cos.is_closed, "a closed loop");
            let (t0, t1) = cos.curve_3d.parameter_range();
            let (mut len, mut off) = (0.0_f64, 0.0_f64);
            let mut last = cos.curve_3d.point_at(t0);
            for k in 0..=4000 {
                let p = cos.curve_3d.point_at(t0 + (t1 - t0) * k as f64 / 4000.0);
                len += p.distance_to(&last);
                last = p;
                let rel = p - origin;
                let db = ((rel - an * rel.dot(&an)).length() - 1.0).abs();
                off = off.max(((p.x * p.x + p.y * p.y).sqrt() - 1.0).abs().max(db));
            }
            (len, off)
        })
        .collect()
}

/// Square cylinders overlapping by 1e-5 at an external touch meet in a
/// loop about 9e-3 rad wide round A: narrower than the discriminant
/// scan's step, so a touch between two samples answered "no
/// intersection" (W3 near-miss cells). Turned round A's axis so the
/// touch falls at every phase of the scan.
#[test]
fn square_cylinders_overlapping_by_a_hair_meet_in_one_small_loop() {
    for k in 0..9 {
        let phi = 0.0123 * k as f64 / 9.0 + 0.4;
        let (c, s) = (phi.pcos(), phi.psin());
        let origin = Point3::new((2.0 - 1e-5) * c, (2.0 - 1e-5) * s, 0.0);
        let loops = unit_pair_loops(origin, Vec3::new(-s, c, 0.0));
        assert_eq!(loops.len(), 1, "phase {k}: one loop, got {loops:?}");
        let (len, off) = loops[0];
        assert!(len < 0.05, "phase {k}: a small loop, length {len:.3e}");
        assert!(off < 1e-7, "phase {k}: on both cylinders, off {off:.2e}");
    }
}

/// Nearly parallel cylinders overlapping by 1e-5 meet in a long thin loop
/// (0.25 by 6e-3) through θ = 0 of A's scan: the entry bracket [θ_511, 0]
/// sent the bisection across the whole turn and the loop swept nearly a
/// full turn of negative discriminants, fitted 160 long (W3 near-miss
/// cells).
#[test]
fn nearly_parallel_cylinders_overlapping_by_a_hair_meet_in_one_thin_loop() {
    let alpha: f64 = 0.05;
    let loops = unit_pair_loops(
        Point3::new(2.0 - 1e-5, 0.0, 0.0),
        Vec3::new(0.0, alpha.psin(), alpha.pcos()),
    );
    assert_eq!(loops.len(), 1, "one loop, got {loops:?}");
    let (len, off) = loops[0];
    // The relative curvature's eigenvalues ≈ 2 and sin²α / 2: a thin
    // ellipse of semi-axes √(2δ/λ), its perimeter ≈ four long semi-axes.
    let long = (2.0 * 1e-5 / (alpha.psin() * alpha.psin() / 2.0)).sqrt();
    assert!(
        (len - 4.0 * long).abs() < 0.05 * 4.0 * long,
        "a thin loop, length {len:.4} for {:.4}",
        4.0 * long
    );
    assert!(off < 1e-7, "on both cylinders, off {off:.2e}");
}

/// Equal cylinders with crossing axes (two ellipses crossing at two tangent
/// points) shifted apart by 2e-6 across one axis: at one tangency the
/// rulings meet B in roots 4e-3 apart, at the other in none over a dip of
/// the discriminant about 4e-3 rad wide, between two samples of the scan.
/// Read as a full pierce, the two root loops touched at the dip (clamped to
/// its double root), the boolean met them in a pinch vertex and declined a
/// valid difference (scene 1512 of the generated corpus, W3 V cells); the
/// curve is one loop turning on both sides of the dip.
#[test]
fn shifted_crossing_cylinders_meet_in_one_loop_past_a_sub_sample_dip() {
    for k in 0..7 {
        let phi = 0.0123 * k as f64 / 7.0 + 0.3;
        let (c, s) = (phi.pcos(), phi.psin());
        // B's axis along x turned by φ round z, shifted 2e-6 along y turned.
        let origin = Point3::new(-2e-6 * s, 2e-6 * c, 0.0);
        let loops = unit_pair_loops(origin, Vec3::new(c, s, 0.0));
        assert_eq!(loops.len(), 1, "phase {k}: one loop, got {loops:?}");
        let (len, off) = loops[0];
        // Two ellipses of semi-axes 1 and √2: 7.64 each.
        assert!((len - 15.28).abs() < 0.05, "phase {k}: length {len:.4}");
        assert!(off < 1e-7, "phase {k}: on both cylinders, off {off:.2e}");
    }
}
