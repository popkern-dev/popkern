#[cfg(test)]
use super::*;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[test]
fn rotate_90_degrees() {
    let p = Point3::new(1.0, 0.0, 0.0);
    let r = rotate_point(
        &p,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    );
    assert!((r.x).abs() < 1e-10);
    assert!((r.y - 1.0).abs() < 1e-10);
}

#[test]
fn rotate_360_returns_to_start() {
    let p = Point3::new(1.0, 0.0, 0.0);
    let r = rotate_point(
        &p,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
    );
    assert!((r.x - 1.0).abs() < 1e-10);
    assert!((r.y).abs() < 1e-10);
}

#[test]
fn rotate_around_offset_axis() {
    let p = Point3::new(2.0, 0.0, 0.0);
    let r = rotate_point(
        &p,
        &Point3::new(1.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
    );
    assert!((r.x).abs() < 1e-10);
    assert!((r.y).abs() < 1e-10);
}

#[test]
fn revolution_curve_full_is_circle() {
    let c = revolution_curve(
        &Point3::new(1.0, 0.0, 0.0),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
    );
    assert!(matches!(c, CurveGeometry::Circle { .. }));
}

#[test]
fn revolution_curve_partial_is_arc() {
    let c = revolution_curve(
        &Point3::new(1.0, 0.0, 0.0),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    );
    assert!(matches!(c, CurveGeometry::Arc { .. }));
}

#[test]
fn revolution_curve_radius_correct() {
    let c = revolution_curve(
        &Point3::new(3.0, 0.0, 0.0),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
    );
    if let CurveGeometry::Circle { radius, .. } = c {
        assert!((radius - 3.0).abs() < 1e-10);
    } else {
        panic!("Expected Circle");
    }
}

#[test]
fn revolution_curve_offset_point() {
    let c = revolution_curve(
        &Point3::new(2.0, 0.0, 3.0),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
    );
    if let CurveGeometry::Circle { center, radius, .. } = c {
        assert!((center.z - 3.0).abs() < 1e-10);
        assert!((radius - 2.0).abs() < 1e-10);
    } else {
        panic!("Expected Circle");
    }
}

#[test]
fn revolution_surface_cylinder() {
    use crate::geometry::surface::SurfaceEval;
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    let surface = revolution_surface(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
    )
    .expect("should construct cylinder");

    for u in [0.0, 0.5, 0.99] {
        for v_frac in [0.0, 0.25, 0.5, 0.75] {
            let (v_min, v_max) = (
                surface.knots_v[surface.degree_v],
                surface.knots_v[surface.knots_v.len() - surface.degree_v - 1],
            );
            let v = v_min + v_frac * (v_max - v_min);
            let p = surface.point_at(u, v);
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (r - 1.0).abs() < 0.15,
                "Not on cylinder: r={r} at u={u}, v={v}"
            );
        }
    }
}

#[test]
fn revolution_surface_partial_90() {
    use crate::geometry::surface::SurfaceEval;
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    let surface = revolution_surface(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        FRAC_PI_2,
    )
    .expect("should construct quarter-cylinder");

    let p_start = surface.point_at(0.0, surface.knots_v[surface.degree_v]);
    assert!(
        (p_start.x - 1.0).abs() < 0.15,
        "Start should be at x=1, got x={}",
        p_start.x
    );
}

#[test]
fn revolution_surface_invalid_angle() {
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(
        revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            0.0
        )
        .is_err()
    );
}

#[test]
fn revolution_surface_negative_angle_fails() {
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(
        revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            -1.0
        )
        .is_err()
    );
}

#[test]
fn revolution_surface_zero_axis_fails() {
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(
        revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 0.0),
            PI
        )
        .is_err()
    );
}

#[test]
fn revolution_surface_180_degrees() {
    use crate::geometry::surface::SurfaceEval;
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    let surface = revolution_surface(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        PI,
    )
    .expect("should construct half-cylinder");
    let p = surface.point_at(0.5, PI / 2.0);
    let r = (p.x * p.x + p.y * p.y).sqrt();
    assert!((r - 1.0).abs() < 0.15, "Not on cylinder: r={r}");
}

#[test]
fn revolution_surface_270_degrees() {
    use crate::geometry::surface::SurfaceEval;
    let angle = 3.0 * FRAC_PI_2;
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
        weights: vec![1.0, 1.0],
    };
    let surface = revolution_surface(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        angle,
    )
    .expect("should construct 3/4 cylinder");
    let p = surface.point_at(0.5, angle / 2.0);
    let r = (p.x * p.x + p.y * p.y).sqrt();
    assert!((r - 1.0).abs() < 0.15, "Not on cylinder: r={r}");
}

#[test]
fn revolution_surface_point_on_axis_degenerate() {
    let profile = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0), // on axis
            Point3::new(1.0, 0.0, 1.0),
        ],
        weights: vec![1.0, 1.0],
    };
    let surface = revolution_surface(
        &profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        TAU,
    );
    assert!(surface.is_ok());
}
