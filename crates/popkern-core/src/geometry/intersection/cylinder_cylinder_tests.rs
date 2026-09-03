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
