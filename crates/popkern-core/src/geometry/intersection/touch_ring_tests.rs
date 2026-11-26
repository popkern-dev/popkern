use super::is_touch_band;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

const TOL: f64 = 1e-6;

fn floor() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(-5.0, -5.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn ball(sunk: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 1.0 - sunk),
        radius: 1.0,
        axis: Vec3::new(1.0, 0.0, 0.0),
    }
}

/// The circle of radius `rho` in the floor about the origin.
fn ring(rho: f64) -> CurveOnSurface {
    let c = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: rho,
    }
    .to_nurbs_in_range((0.0, std::f64::consts::TAU));
    CurveOnSurface {
        curve_3d: c.clone(),
        curve_uv1: c.clone(),
        curve_uv2: c,
        is_closed: true,
    }
}

fn band(sunk: f64, rho: f64) -> bool {
    let w_floor = ((0.0, 10.0), (0.0, 10.0));
    let b = ball(sunk);
    let w_ball = {
        use crate::geometry::surface::SurfaceEval as _;
        b.parameter_range()
    };
    is_touch_band(
        &ring(rho),
        (&floor(), w_floor),
        (&b, w_ball),
        (TOL, TOL * (1.0 + 1e-9), 10.0),
    )
}

/// A ball resting on a floor: the circle where their gap is the
/// tolerance, √(2·r·ε) about the touch, is the touch.
#[test]
fn the_tolerance_ring_round_a_resting_ball_is_its_touch() {
    assert!(band(0.0, (2.0 * TOL).sqrt()));
}

/// Negative control: sunk ten tolerances, the ball cuts a real lens,
/// and its rim (√(2·r·δ) across) is a crossing.
#[test]
fn the_rim_of_a_lens_ten_tolerances_deep_is_a_crossing() {
    assert!(!band(10.0 * TOL, (2.0 * 10.0 * TOL).sqrt()));
}

/// A circle far wider than the touch's band is no part of it.
#[test]
fn a_wide_circle_about_a_touch_is_no_band() {
    assert!(!band(0.0, 0.5));
}
