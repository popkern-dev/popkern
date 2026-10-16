//! Coaxial torus analytic-arm tests: exact circles, tangency-silent
//! band, coincident deferral, fallback gates, window honouring, and
//! operand-order equivalence (PK-279).
use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::intersect_surfaces;
use crate::geometry::surface::SurfaceEval;
use std::f64::consts::TAU;

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

/// Max residual of `curve` against both surfaces plus its
/// co-parameterization error, over a dense parameter sweep.
fn curve_residuals(
    curve: &crate::geometry::curve_on_surface::CurveOnSurface,
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
) -> (f64, f64) {
    let mut worst_surface = 0.0f64;
    let mut worst_sync = 0.0f64;
    for i in 0..=256 {
        let t = i as f64 / 256.0;
        let p = curve.curve_3d.point_at(t);
        let uv1 = curve.curve_uv1.point_at(t);
        let uv2 = curve.curve_uv2.point_at(t);
        let p1 = s1.point_at(uv1.x, uv1.y);
        let p2 = s2.point_at(uv2.x, uv2.y);
        worst_sync = worst_sync.max(p1.distance_to(&p)).max(p2.distance_to(&p));
        worst_surface = worst_surface
            .max(implicit_residual(s1, &p))
            .max(implicit_residual(s2, &p));
    }
    (worst_surface, worst_sync)
}

/// Closed-form distance-to-surface residual for the analytic types
/// used in these tests.
fn implicit_residual(surface: &SurfaceGeometry, p: &Point3) -> f64 {
    match surface {
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let d = *p - *center;
            let h = d.dot(axis);
            let rho = (d - *axis * h).length();
            (((rho - major_radius) * (rho - major_radius) + h * h).sqrt() - minor_radius).abs()
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let d = *p - *origin;
            let h = d.dot(axis);
            ((d - *axis * h).length() - radius).abs()
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let d = *p - *apex;
            let h = d.dot(axis);
            let rho = (d - *axis * h).length();
            (rho - half_angle.tan().abs() * h.abs()).abs() * half_angle.cos()
        }
        SurfaceGeometry::Sphere { center, radius, .. } => ((*p - *center).length() - radius).abs(),
        SurfaceGeometry::Plane { origin, normal, .. } => (*p - *origin).dot(normal).abs(),
        _ => f64::INFINITY,
    }
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

#[test]
fn coaxial_cylinder_meets_tube_in_two_exact_circles() {
    let tol = Tolerance::default();
    let torus = torus();
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let out = try_torus_coaxial(&torus, full(), full(), &cyl, full(), (-6.0, 6.0), &tol)
        .expect("coaxial pair must be claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2, "cylinder through the tube: two circles");
    for c in &curves {
        assert!(c.is_closed);
        assert_eq!(c.curve_3d.degree, 2, "exact rational circle expected");
        let (surf, sync) = curve_residuals(c, &torus, &cyl);
        assert!(surf < 1e-9, "on-surface residual {surf}");
        assert!(sync < 5e-6, "pcurve sync residual {sync}");
    }
}

#[test]
fn tangent_and_separated_cylinders_are_silent() {
    let tol = Tolerance::default();
    let torus = torus();
    for rc in [2.5, 1.5, 2.6, 1.2] {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: rc,
        };
        let out = try_torus_coaxial(&torus, full(), full(), &cyl, full(), (-6.0, 6.0), &tol)
            .expect("coaxial pair must be claimed");
        assert!(
            curves_of(out).is_empty(),
            "rc={rc}: tangent/separated must be silent"
        );
    }
}

#[test]
fn skew_or_offset_axes_fall_back_to_marching() {
    let tol = Tolerance::default();
    let torus = torus();
    let skew = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 1.0).normalized().unwrap(),
        radius: 2.0,
    };
    assert!(try_torus_coaxial(&torus, full(), full(), &skew, full(), (-6.0, 6.0), &tol).is_none());
    let offset = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.5, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    assert!(
        try_torus_coaxial(&torus, full(), full(), &offset, full(), (-6.0, 6.0), &tol).is_none()
    );
}

#[test]
fn stacked_coaxial_tori_meet_in_two_circles() {
    let tol = Tolerance::default();
    let t1 = torus();
    let t2 = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.6),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    let out = try_torus_coaxial(&t1, full(), full(), &t2, full(), full(), &tol)
        .expect("coaxial tori must be claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2);
    for c in &curves {
        let (surf, sync) = curve_residuals(c, &t1, &t2);
        assert!(surf < 1e-9, "on-surface residual {surf}");
        assert!(sync < 5e-6, "pcurve sync residual {sync}");
    }
}

#[test]
fn coincident_tori_defer_to_the_coincident_detector() {
    let tol = Tolerance::default();
    let t1 = torus();
    let t2 = torus();
    assert!(
        try_torus_coaxial(&t1, full(), full(), &t2, full(), full(), &tol).is_none(),
        "identical tori must decline so the coincident detector owns the verdict"
    );
    // End-to-end: the dispatcher must still report the coincidence.
    let result = intersect_surfaces(&t1, full(), full(), &t2, full(), full(), &tol)
        .expect("intersect_surfaces");
    assert!(
        result
            .intersections
            .iter()
            .any(|ix| matches!(ix, SurfaceSurfaceIntersection::CoincidentRegion { .. })),
        "coincident coaxial tori must yield a CoincidentRegion"
    );
}

#[test]
fn coaxial_cone_countersink_meets_rim_fillet_in_two_circles() {
    let tol = Tolerance::default();
    let torus = torus();
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 1.5),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    let out = try_torus_coaxial(&torus, full(), full(), &cone, full(), (-6.0, 6.0), &tol)
        .expect("coaxial cone must be claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2, "45° coaxial cone: two circles");
    for c in &curves {
        let (surf, sync) = curve_residuals(c, &torus, &cone);
        assert!(surf < 1e-9, "on-surface residual {surf}");
        assert!(sync < 5e-6, "pcurve sync residual {sync}");
    }
}

#[test]
fn axis_centred_sphere_meets_torus_in_two_circles() {
    let tol = Tolerance::default();
    let torus = torus();
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.3),
        radius: 2.3,
        axis: crate::geometry::Vec3::Z,
    };
    let out = try_torus_coaxial(
        &torus,
        full(),
        full(),
        &sphere,
        full(),
        (0.0, std::f64::consts::PI),
        &tol,
    )
    .expect("axis-centred sphere must be claimed");
    let curves = curves_of(out);
    assert_eq!(curves.len(), 2);
    for c in &curves {
        let (surf, sync) = curve_residuals(c, &torus, &sphere);
        assert!(surf < 1e-9, "on-surface residual {surf}");
        assert!(sync < 5e-6, "pcurve sync residual {sync}");
    }
}

#[test]
fn tube_window_filters_stations_outside_it() {
    let tol = Tolerance::default();
    let torus = torus();
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    // Stations for rc = R sit at tube angles π/2 and 3π/2. A window
    // around π/2 only must keep exactly one circle.
    let window = (1.0, 2.0);
    let out = try_torus_coaxial(&torus, full(), window, &cyl, full(), (-6.0, 6.0), &tol)
        .expect("claimed");
    assert_eq!(curves_of(out).len(), 1, "only the in-window station serves");
}

#[test]
fn operand_order_is_equivalent() {
    let tol = Tolerance::default();
    let torus = torus();
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let ab = curves_of(
        try_torus_coaxial(&torus, full(), full(), &cyl, full(), (-6.0, 6.0), &tol).unwrap(),
    );
    let ba = curves_of(
        try_torus_coaxial(&cyl, full(), (-6.0, 6.0), &torus, full(), full(), &tol).unwrap(),
    );
    assert_eq!(ab.len(), ba.len());
    // Same 3D circle set; pcurves attach to the right operand.
    for c in &ba {
        let (surf, sync) = curve_residuals(c, &cyl, &torus);
        assert!(surf < 1e-9);
        assert!(sync < 5e-6);
    }
    for a in &ab {
        let pa = a.curve_3d.point_at(0.37);
        let best = ba
            .iter()
            .map(|b| {
                (0..=4096)
                    .map(|i| b.curve_3d.point_at(i as f64 / 4096.0).distance_to(&pa))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(f64::INFINITY, f64::min);
        // Gate at the sampling spacing of the densely-sampled circle.
        assert!(
            best < 5e-3,
            "curve sets differ between operand orders: {best}"
        );
    }
}
