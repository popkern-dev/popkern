//! Pins for the measured recompute interpolant (PK-555).

use super::*;
use crate::geometry::Vec3;

fn plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// A quadratic B-spline with a double interior knot at 0.5: two arcs
/// meeting with a parametric kink (C0), like a conic rim written as
/// rational quadratic segments.
fn kinked_curve() -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.2, -3.0, 0.0),
                Point3::new(5.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 5],
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// Plane uv of the edge curve at each fraction (the exact projection).
fn sampler(curve: &CurveOnEdge) -> impl Fn(&[f64]) -> Vec<Point3> + '_ {
    move |fractions: &[f64]| {
        fractions
            .iter()
            .map(|&f| {
                let p = curve.geometry.point_at(f);
                Point3::new(p.x, p.y, 0.0)
            })
            .collect()
    }
}

#[test]
fn breakpoints_are_the_curves_c0_knots_inside_the_range() {
    let curve = kinked_curve();
    assert_eq!(breakpoint_fractions(&curve, (0.0, 1.0)), vec![0.5]);
    // Reversed range: fraction measured from t_lo.
    assert_eq!(breakpoint_fractions(&curve, (1.0, 0.0)), vec![0.5]);
    // A sub-range that excludes the knot has none.
    assert!(breakpoint_fractions(&curve, (0.0, 0.4)).is_empty());
}

#[test]
fn the_interpolant_follows_a_kinked_edge_within_the_target() {
    let curve = kinked_curve();
    let surface = plane();
    let linear = 1e-6;
    let fit =
        measured_projection_interpolant(&surface, &curve, (0.0, 1.0), linear, &sampler(&curve))
            .expect("fit");
    // Co-parameterized: the image at each fraction is the edge point.
    let worst = (0..=400)
        .map(|k| {
            let f = k as f64 / 400.0;
            let uv = fit.point_at(f);
            surface
                .point_at(uv.x, uv.y)
                .distance_to(&curve.geometry.point_at(f))
        })
        .fold(0.0_f64, f64::max);
    assert!(worst <= 0.1 * linear, "kinked edge followed to {worst:e}");

    // The historical single smooth interpolant through 23 uniform
    // samples overshoots at the kink by orders of magnitude more.
    let uniform = uniform_interpolant(&sampler(&curve)(&segment_fractions(
        &[],
        &[PROJECTION_SAMPLES - 1],
    )))
    .expect("uniform");
    let historical = (0..=400)
        .map(|k| {
            let f = k as f64 / 400.0;
            let uv = uniform.point_at(f);
            Point3::new(uv.x, uv.y, 0.0).distance_to(&curve.geometry.point_at(f))
        })
        .fold(0.0_f64, f64::max);
    assert!(historical > 1e3 * linear, "historical error {historical:e}");
}
