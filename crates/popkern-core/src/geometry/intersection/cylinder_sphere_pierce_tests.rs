// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

//! The off-axis cylinder×sphere pierce: the branch count of each regime
//! and every curve, and both its pcurves' images, on both surfaces.

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::cylinder_sphere::try_cylinder_sphere;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::{PI, TAU};

fn tol() -> Tolerance {
    Tolerance::default().with_linear(1e-6).with_angular(1e-6)
}

/// A cylinder of radius 6 along y through (x, 0, 0) against a ball of
/// radius 20 on the origin (its pole along `pole`).
fn pierce(x: f64, pole: Vec3) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(x, -30.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        radius: 6.0,
    };
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 20.0,
        axis: pole,
    };
    try_cylinder_sphere(
        &cyl,
        ((0.0, TAU), (0.0, 60.0)),
        &ball,
        ((0.0, TAU), (0.0, PI)),
        &tol(),
    )
}

/// Every curve lies on both surfaces, and so do its pcurves' images.
fn assert_on_both(curves: &[SurfaceSurfaceIntersection], x: f64, pole: Vec3) {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(x, -30.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        radius: 6.0,
    };
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 20.0,
        axis: pole,
    };
    for c in curves {
        let SurfaceSurfaceIntersection::Curve(c) = c else {
            panic!("a point where a curve belongs");
        };
        let (t0, t1) = c.curve_3d.parameter_range();
        let (a0, a1) = c.curve_uv1.parameter_range();
        let (b0, b1) = c.curve_uv2.parameter_range();
        for k in 0..=400 {
            let f = k as f64 / 400.0;
            let p = c.curve_3d.point_at(t0 + (t1 - t0) * f);
            let rho = ((p.x - x).powi(2) + p.z * p.z).sqrt();
            assert!(
                (rho - 6.0).abs() < 1e-6,
                "off the cylinder by {}",
                rho - 6.0
            );
            let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!((r - 20.0).abs() < 1e-6, "off the ball by {}", r - 20.0);
            let q1 = c.curve_uv1.point_at(a0 + (a1 - a0) * f);
            let q2 = c.curve_uv2.point_at(b0 + (b1 - b0) * f);
            assert!(cyl.point_at(q1.x, q1.y).distance_to(&p) < 1e-6);
            assert!(ball.point_at(q2.x, q2.y).distance_to(&p) < 1e-6);
        }
    }
}

#[test]
fn a_bore_through_the_ball_off_its_centre_meets_it_in_two_loops() {
    let pole = Vec3::new(0.0, 0.0, 1.0);
    for x in [0.5, 5.0, 8.0] {
        let curves = pierce(x, pole).expect("answered across the pole axis");
        assert_eq!(curves.len(), 2, "x {x}");
        assert_on_both(&curves, x, pole);
    }
}

/// Along the pole axis a bore holding the pole winds its loops round
/// the chart's poles: the scan leaves those to the marching path. One
/// beside the pole (x 8, the pole 2 outside its wall) stays in the
/// chart and is answered on both surfaces.
#[test]
fn a_bore_round_the_pole_is_left_to_the_marching_path() {
    let pole = Vec3::new(0.0, 1.0, 0.0);
    for x in [0.5, 5.0] {
        assert!(pierce(x, pole).is_none(), "x {x}");
    }
    let curves = pierce(8.0, pole).expect("beside the pole");
    assert_eq!(curves.len(), 2);
    assert_on_both(&curves, 8.0, pole);
}

/// Breaking out of the ball's side (|x| + r > R): one loop round the
/// break-out.
#[test]
fn a_bore_breaking_out_of_the_side_meets_it_in_one_loop() {
    let curves = pierce(16.0, Vec3::new(0.0, 0.0, 1.0)).expect("answered");
    assert_eq!(curves.len(), 1);
    assert_on_both(&curves, 16.0, Vec3::new(0.0, 0.0, 1.0));
}

#[test]
fn a_bore_beside_the_ball_misses_it() {
    assert!(
        pierce(27.0, Vec3::new(0.0, 0.0, 1.0))
            .expect("answered")
            .is_empty()
    );
}
