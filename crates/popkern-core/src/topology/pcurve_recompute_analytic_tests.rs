//! Tests for exact analytic pcurve candidates (PK-210 arcs, PK-287
//! circles and NURBS curves in plane faces).

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// PK-287: the exact NURBS-in-plane candidate must reproduce a
/// rational NURBS circle to machine precision at MATCHING parameters
/// — trace AND parameterization, including the rational weights.
#[test]
fn exact_plane_nurbs_candidate_is_parameter_exact() {
    use super::exact_plane_pcurve;
    use crate::geometry::Vec3;
    use crate::geometry::curve::{CurveEval, CurveGeometry};
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI, TAU};

    let r = 0.3;
    let w = FRAC_1_SQRT_2;
    let circle = crate::geometry::nurbs_curve::NurbsCurve {
        degree: 2,
        knots: vec![
            0.0,
            0.0,
            0.0,
            FRAC_PI_2,
            FRAC_PI_2,
            PI,
            PI,
            PI + FRAC_PI_2,
            PI + FRAC_PI_2,
            TAU,
            TAU,
            TAU,
        ],
        control_points: vec![
            Point3::new(0.5 + r, 0.5, 1.0),
            Point3::new(0.5 + r, 0.5 + r, 1.0),
            Point3::new(0.5, 0.5 + r, 1.0),
            Point3::new(0.5 - r, 0.5 + r, 1.0),
            Point3::new(0.5 - r, 0.5, 1.0),
            Point3::new(0.5 - r, 0.5 - r, 1.0),
            Point3::new(0.5, 0.5 - r, 1.0),
            Point3::new(0.5 + r, 0.5 - r, 1.0),
            Point3::new(0.5 + r, 0.5, 1.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    };
    let curve = CurveGeometry::Nurbs(circle);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let pcurve = exact_plane_pcurve(&plane, &curve, 0.0, TAU).expect("nurbs-in-plane candidate");
    for i in 0..=64 {
        let t = TAU * i as f64 / 64.0;
        let on_curve = curve.point_at(t);
        let uv = pcurve.point_at(t);
        let on_trace = plane.point_at(uv.x, uv.y);
        assert!(
            on_curve.distance_to(&on_trace) < 1e-12,
            "t={t}: parameter-matched deviation {:.3e}",
            on_curve.distance_to(&on_trace)
        );
    }
}

/// PK-287: full Circle curves in a plane get the same exact
/// angle-parameterized elliptical candidate as Arcs.
#[test]
fn exact_plane_circle_candidate_is_parameter_exact() {
    use super::exact_plane_pcurve;
    use crate::geometry::Vec3;
    use crate::geometry::curve::{CurveEval, CurveGeometry};
    use std::f64::consts::TAU;

    let curve = CurveGeometry::Circle {
        center: Point3::new(0.5, 0.5, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.3,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let pcurve = exact_plane_pcurve(&plane, &curve, 0.0, TAU).expect("circle-in-plane candidate");
    for i in 0..=32 {
        let t = TAU * i as f64 / 32.0;
        let on_curve = curve.point_at(t);
        let uv = pcurve.point_at(t);
        let on_trace = plane.point_at(uv.x, uv.y);
        assert!(
            on_curve.distance_to(&on_trace) < 1e-12,
            "t={t}: parameter-matched deviation {:.3e}",
            on_curve.distance_to(&on_trace)
        );
    }
}
