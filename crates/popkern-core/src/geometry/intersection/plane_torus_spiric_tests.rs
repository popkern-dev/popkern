//! Spiric-section sampler tests: loop-count structure, on-surface and
//! co-parameterization residuals, a line×torus oracle for branch
//! completeness, and dispatch integration (PK-279).
use super::*;
use crate::geometry::Vec3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::intersection::line_torus::intersect_line_torus;
use crate::geometry::intersection::plane_torus::try_plane_torus;
use crate::geometry::intersection::surface_surface::intersect_surfaces;
use crate::geometry::surface::SurfaceEval;

fn torus() -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    }
}

fn full() -> (f64, f64) {
    (0.0, TAU)
}

fn torus_residual(p: &Point3) -> f64 {
    let rho = (p.x * p.x + p.y * p.y).sqrt();
    (((rho - 2.0) * (rho - 2.0) + p.z * p.z).sqrt() - 0.5).abs()
}

fn curves_of(
    result: Vec<SurfaceSurfaceIntersection>,
) -> Vec<crate::geometry::curve_on_surface::CurveOnSurface> {
    result
        .into_iter()
        .filter_map(|ix| match ix {
            SurfaceSurfaceIntersection::Curve(c) => Some(c),
            SurfaceSurfaceIntersection::CoincidentRegion { .. } => None,
        })
        .collect()
}

fn verify(
    curves: &[crate::geometry::curve_on_surface::CurveOnSurface],
    plane: &SurfaceGeometry,
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
) {
    let SurfaceGeometry::Plane { origin, normal, .. } = plane else {
        panic!("plane expected");
    };
    for c in curves {
        assert!(c.is_closed, "spiric loops are closed");
        for i in 0..=256 {
            let t = i as f64 / 256.0;
            let p = c.curve_3d.point_at(t);
            let plane_err = (p - *origin).dot(normal).abs();
            assert!(plane_err < 1e-6, "off plane by {plane_err}");
            assert!(
                torus_residual(&p) < 1e-6,
                "off torus by {}",
                torus_residual(&p)
            );
            let uv1 = c.curve_uv1.point_at(t);
            let uv2 = c.curve_uv2.point_at(t);
            assert!(s1.point_at(uv1.x, uv1.y).distance_to(&p) < 1e-5);
            assert!(s2.point_at(uv2.x, uv2.y).distance_to(&p) < 1e-5);
        }
    }
}

/// Oracle: shoot a fan of in-plane lines through the configuration;
/// every line×torus hit must lie within tolerance of some emitted
/// curve. This checks branch COMPLETENESS independently of the
/// solver's own machinery.
fn oracle_check(
    plane_origin: Point3,
    axis1: Vec3,
    axis2: Vec3,
    curves: &[crate::geometry::curve_on_surface::CurveOnSurface],
) {
    let tol = Tolerance::default();
    let mut checked = 0;
    for k in 0..48 {
        let ang = TAU * k as f64 / 48.0;
        let dir = axis1 * ang.cos() + axis2 * ang.sin();
        let hits = intersect_line_torus(
            &plane_origin,
            &dir,
            (-10.0, 10.0),
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            2.0,
            0.5,
            (0.0, TAU),
            (0.0, TAU),
            &tol,
        );
        for hit in &hits {
            let SurfaceCurveIntersection::Point { point: p, .. } = hit else {
                continue;
            };
            let p = *p;
            // Point-to-SEGMENT distance against a dense polyline of
            // each curve (point-to-sample would gate on the sampling
            // density instead of the curve).
            let best = curves
                .iter()
                .map(|c| {
                    let pts: Vec<Point3> = (0..=1024)
                        .map(|i| c.curve_3d.point_at(i as f64 / 1024.0))
                        .collect();
                    pts.windows(2)
                        .map(|seg| {
                            let ab = seg[1] - seg[0];
                            let len_sq = ab.dot(&ab);
                            let t = if len_sq > 0.0 {
                                ((p - seg[0]).dot(&ab) / len_sq).clamp(0.0, 1.0)
                            } else {
                                0.0
                            };
                            (seg[0] + ab * t).distance_to(&p)
                        })
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(f64::INFINITY, f64::min);
            assert!(
                best < 1e-3,
                "oracle point {p:?} is {best} from every emitted curve — a branch is missing"
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "oracle produced no intersection points");
}

#[test]
fn oblique_plane_produces_verified_spiric_loops() {
    let tol = Tolerance::default();
    let torus = torus();
    let n = Vec3::new(0.5, 0.0, 0.75_f64.sqrt()); // 30° off the axis
    let u_axis = Vec3::new(0.75_f64.sqrt(), 0.0, -0.5);
    let v_axis = Vec3::new(0.0, 1.0, 0.0);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.1, 0.0, 0.1),
        normal: n,
        u_axis,
        v_axis,
    };
    let out = try_plane_torus(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("plane×torus must be claimed");
    let curves = curves_of(out);
    assert!(!curves.is_empty(), "oblique cut must produce loops");
    verify(&curves, &plane, &plane, &torus);
    oracle_check(Point3::new(0.1, 0.0, 0.1), u_axis, v_axis, &curves);
}

#[test]
fn axis_parallel_plane_missing_the_hole_gives_one_loop() {
    let tol = Tolerance::default();
    let torus = torus();
    // x = 2 passes through the tube centre ring: a single loop.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(2.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let out = try_plane_torus(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 1, "tube-centre cut is one loop");
    verify(&curves, &plane, &plane, &torus);
    oracle_check(
        Point3::new(2.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        &curves,
    );
}

#[test]
fn axis_parallel_plane_through_the_hole_gives_two_loops() {
    let tol = Tolerance::default();
    let torus = torus();
    // x = 1 crosses the ring on both sides: two loops.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let out = try_plane_torus(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2, "hole-crossing cut is two loops");
    verify(&curves, &plane, &plane, &torus);
    oracle_check(
        Point3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        &curves,
    );
}

#[test]
fn near_perpendicular_plane_crosses_over_to_two_full_loops() {
    let tol = Tolerance::default();
    let torus = torus();
    // Just outside the perpendicular angular gate: the spiric arm's
    // root-free two-loop case (the analytic continuation of the
    // perpendicular circle pair).
    let tilt = 2.0e-3_f64;
    let n = Vec3::new(tilt.sin(), 0.0, tilt.cos());
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.1),
        normal: n,
        u_axis: Vec3::new(tilt.cos(), 0.0, -tilt.sin()),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let out = try_plane_torus(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2, "near-perpendicular offset cut: two loops");
    verify(&curves, &plane, &plane, &torus);
}

#[test]
fn separated_plane_is_claimed_empty() {
    let tol = Tolerance::default();
    let torus = torus();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(4.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let out = try_plane_torus(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("claimed");
    assert!(curves_of(out).is_empty());
}

/// End-to-end dispatch: `intersect_surfaces` must serve the general
/// spiric case through the analytic arm (fitted cubic loops, correct
/// count) in both operand orders.
#[test]
fn dispatch_serves_spiric_sections_in_both_operand_orders() {
    let tol = Tolerance::default();
    let torus = torus();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let pf = intersect_surfaces(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        full(),
        full(),
        &tol,
    )
    .expect("plane-first");
    let tf = intersect_surfaces(
        &torus,
        full(),
        full(),
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &tol,
    )
    .expect("torus-first");
    let pf_curves = curves_of(pf.intersections);
    let tf_curves = curves_of(tf.intersections);
    assert_eq!(pf_curves.len(), 2);
    assert_eq!(tf_curves.len(), 2);
    for c in &tf_curves {
        // Operand order swaps the pcurve assignment.
        for i in 0..=64 {
            let t = i as f64 / 64.0;
            let p = c.curve_3d.point_at(t);
            let uv1 = c.curve_uv1.point_at(t);
            assert!(torus.point_at(uv1.x, uv1.y).distance_to(&p) < 1e-5);
        }
    }
}
