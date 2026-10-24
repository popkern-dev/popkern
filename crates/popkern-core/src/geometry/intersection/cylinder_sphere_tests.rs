//! Tests for the analytic coaxial cylinder–sphere arm (PK-390).

use super::*;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::intersect_surfaces;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::{PI, TAU};

fn tol() -> Tolerance {
    Tolerance::default().with_linear(1e-6).with_angular(1e-6)
}

fn z_cylinder(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

fn sphere(center: Point3, radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center,
        radius,
        axis: crate::geometry::Vec3::Z,
    }
}

fn must_some<T>(opt: Option<T>, msg: &str) -> T {
    match opt {
        Some(v) => v,
        None => panic!("{msg}"),
    }
}

#[test]
fn equal_radii_on_axis_is_tangent_and_silent() {
    // The fillet-strip × corner-sphere junction: same rolling-ball
    // radius, centre on the strip's axis.
    let cyl = z_cylinder(0.08);
    let sph = sphere(Point3::new(0.0, 0.0, 0.32), 0.08);
    let r = must_some(
        try_cylinder_sphere(
            &cyl,
            ((0.0, TAU), (-1.0, 1.0)),
            &sph,
            ((0.0, TAU), (0.0, PI)),
            &tol(),
        ),
        "coaxial pair is analytic",
    );
    assert!(r.is_empty(), "tangent contact emits nothing");
    let r = must_some(
        try_cylinder_sphere(
            &sph,
            ((0.0, TAU), (0.0, PI)),
            &cyl,
            ((0.0, TAU), (-1.0, 1.0)),
            &tol(),
        ),
        "either operand order",
    );
    assert!(r.is_empty());
}

#[test]
fn smaller_sphere_inside_the_cylinder_is_empty() {
    let cyl = z_cylinder(1.0);
    let sph = sphere(Point3::new(0.0, 0.0, 0.5), 0.4);
    let r = must_some(
        try_cylinder_sphere(
            &cyl,
            ((0.0, TAU), (-2.0, 2.0)),
            &sph,
            ((0.0, TAU), (0.0, PI)),
            &tol(),
        ),
        "coaxial",
    );
    assert!(r.is_empty());
}

#[test]
fn larger_sphere_cuts_two_exact_circles() {
    let cyl = z_cylinder(0.6);
    let sph = sphere(Point3::new(0.0, 0.0, 0.5), 1.0);
    let r = must_some(
        try_cylinder_sphere(
            &cyl,
            ((0.0, TAU), (-2.0, 3.0)),
            &sph,
            ((0.0, TAU), (0.0, PI)),
            &tol(),
        ),
        "coaxial",
    );
    assert_eq!(r.len(), 2);
    let h = (1.0f64 - 0.36).sqrt();
    let mut stations: Vec<f64> = Vec::new();
    for curve in &r {
        let SurfaceSurfaceIntersection::Curve(c) = curve else {
            panic!("expected a curve");
        };
        assert!(c.is_closed);
        for i in 0..=16 {
            let t = i as f64 / 16.0;
            let p = c.curve_3d.point_at(t);
            let radial = (p.x * p.x + p.y * p.y).sqrt();
            assert!((radial - 0.6).abs() < 1e-12, "on the cylinder: {radial}");
            let d = p.distance_to(&Point3::new(0.0, 0.0, 0.5));
            assert!((d - 1.0).abs() < 1e-12, "on the sphere: {d}");
            // UV images trace the same circle (image equality; the
            // rational-quadratic 3-D parameter is not the chart angle,
            // as for every analytic circle arm).
            let on_circle = |q: Point3| {
                let radial = (q.x * q.x + q.y * q.y).sqrt();
                (radial - 0.6).abs() < 1e-9 && (q.z - p.z).abs() < 1e-9
            };
            let uv1 = c.curve_uv1.point_at(t);
            assert!(on_circle(cyl.point_at(uv1.x, uv1.y)), "cylinder pcurve");
            let uv2 = c.curve_uv2.point_at(t);
            assert!(on_circle(sph.point_at(uv2.x, uv2.y)), "sphere pcurve");
        }
        stations.push(c.curve_3d.point_at(0.0).z);
    }
    stations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    assert!((stations[0] - (0.5 - h)).abs() < 1e-12);
    assert!((stations[1] - (0.5 + h)).abs() < 1e-12);
    // Window clip: only the upper circle when the cylinder starts above the lower one.
    let r = must_some(
        try_cylinder_sphere(
            &cyl,
            ((0.0, TAU), (0.0, 3.0)),
            &sph,
            ((0.0, TAU), (0.0, PI)),
            &tol(),
        ),
        "coaxial",
    );
    assert_eq!(r.len(), 1);
}

#[test]
fn off_axis_pairs_keep_the_marching_path() {
    let cyl = z_cylinder(0.6);
    let off = sphere(Point3::new(0.1, 0.0, 0.5), 1.0);
    assert!(
        try_cylinder_sphere(
            &cyl,
            ((0.0, TAU), (-2.0, 3.0)),
            &off,
            ((0.0, TAU), (0.0, PI)),
            &tol()
        )
        .is_none()
    );
    let tilted = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 0.6,
    };
    // Tangent contact is silent regardless of the chart orientation.
    let tangent = sphere(Point3::new(0.5, 0.0, 0.0), 0.6);
    let r = must_some(
        try_cylinder_sphere(
            &tilted,
            ((0.0, TAU), (-2.0, 3.0)),
            &tangent,
            ((0.0, TAU), (0.0, PI)),
            &tol(),
        ),
        "tangent is analytic in any orientation",
    );
    assert!(r.is_empty());
}

#[test]
fn intersect_surfaces_takes_the_arm_for_the_fillet_junction_pair() {
    // The PK-390 pair: strip cylinder r=0.08 with a 20-wide padded
    // window against its corner sphere. Marching seeded ~2100 tangent
    // points and took ~55 s; the arm answers in microseconds.
    let cyl = z_cylinder(0.08);
    let sph = sphere(Point3::new(0.0, 0.0, 0.32), 0.08);
    let t = std::time::Instant::now();
    let r = intersect_surfaces(
        &cyl,
        (0.0, TAU),
        (-6.48, 13.52),
        &sph,
        (0.0, TAU),
        (0.0, PI),
        &tol(),
    );
    let r = match r {
        Ok(r) => r,
        Err(e) => panic!("{e:?}"),
    };
    assert!(r.intersections.is_empty(), "{:?}", r.intersections.len());
    assert!(t.elapsed().as_millis() < 500, "took {:?}", t.elapsed());
}

/// PK-937: a sphere drilled across its pole axis. The circles are as
/// exact as the pole-aligned case: the 3D curves lie on both surfaces
/// to rounding, and the sphere pcurve (a chart fit) maps onto them.
#[test]
fn a_cylinder_across_the_sphere_pole_emits_exact_circles() {
    let x_cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(-1.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 0.2,
    };
    let ball = sphere(Point3::new(0.0, 0.0, 0.0), 0.5);
    let station = 0.21_f64.sqrt();
    for swapped in [false, true] {
        let (cyl_w, sph_w) = (((0.0, TAU), (0.0, 2.0)), ((0.0, TAU), (0.0, PI)));
        let r = must_some(
            if swapped {
                try_cylinder_sphere(&ball, sph_w, &x_cylinder, cyl_w, &tol())
            } else {
                try_cylinder_sphere(&x_cylinder, cyl_w, &ball, sph_w, &tol())
            },
            "coaxial across the pole is analytic",
        );
        assert_eq!(r.len(), 2);
        let mut xs = Vec::new();
        for c in &r {
            let SurfaceSurfaceIntersection::Curve(c) = c else {
                panic!("expected curves");
            };
            let (sph_uv, cyl_uv) = if swapped {
                (&c.curve_uv1, &c.curve_uv2)
            } else {
                (&c.curve_uv2, &c.curve_uv1)
            };
            for k in 0..=64 {
                let t = k as f64 / 64.0;
                let p = c.curve_3d.point_at(t);
                let on_sphere = (p - Point3::new(0.0, 0.0, 0.0)).length() - 0.5;
                let on_cyl = (p.y * p.y + p.z * p.z).sqrt() - 0.2;
                assert!(on_sphere.abs() < 1e-14 && on_cyl.abs() < 1e-14, "t {t}");
                let s = sph_uv.point_at(t);
                let q = ball.point_at(s.x, s.y);
                assert!((q - p).length() < 1e-6, "sphere pcurve at t {t}");
                let s = cyl_uv.point_at(t);
                let q = x_cylinder.point_at(s.x, s.y);
                assert!((q - p).length() < 1e-6, "cylinder pcurve at t {t}");
            }
            xs.push(c.curve_3d.point_at(0.0).x);
        }
        xs.sort_by(f64::total_cmp);
        assert!((xs[0] + station).abs() < 1e-14 && (xs[1] - station).abs() < 1e-14);
    }
}
