//! Plane-through-apex generator lines (PK-931).

use super::super::try_plane_cone;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn cone() -> SurfaceGeometry {
    // BRep::cone(1, 2): apex (0, 0, 2), axis -z, base radius 1.
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.5_f64.atan(),
    }
}

fn plane(normal: Vec3, u_axis: Vec3) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal,
        u_axis,
        v_axis: normal.cross(&u_axis),
    }
}

/// Both generators of `plane` x the cone: every sample of each lies on
/// both surfaces, its cone pcurve has one constant u, and the 3D line
/// and both pcurves agree at matching fractions.
fn assert_generators(plane: &SurfaceGeometry) {
    let tol = Tolerance::default();
    let cone = cone();
    let Some(result) = try_plane_cone(
        plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cone,
        (0.0, std::f64::consts::TAU),
        (0.0, 2.0),
        &tol,
    ) else {
        panic!("plane x cone is analytic");
    };
    assert_eq!(result.len(), 2, "two generators");
    for r in &result {
        let SurfaceSurfaceIntersection::Curve(c) = r else {
            panic!("a curve");
        };
        let u0 = c.curve_uv2.point_at(0.0).x;
        for k in 0..=8 {
            let f = k as f64 / 8.0;
            let p = c.curve_3d.point_at(f);
            let q1 = c.curve_uv1.point_at(f);
            let q2 = c.curve_uv2.point_at(f);
            assert!(plane.point_at(q1.x, q1.y).distance_to(&p) < 1e-12);
            assert!(cone.point_at(q2.x, q2.y).distance_to(&p) < 1e-12);
            assert!(
                (q2.x - u0).abs() < 1e-15,
                "cone u wanders: {} vs {u0}",
                q2.x
            );
        }
        assert!(
            c.curve_3d
                .point_at(0.0)
                .distance_to(&Point3::new(0.0, 0.0, 2.0))
                < 1e-15
        );
    }
}

#[test]
fn a_plane_containing_the_axis_cuts_two_generator_lines() {
    assert_generators(&plane(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)));
    assert_generators(&plane(Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0)));
}

#[test]
fn a_tilted_plane_through_the_apex_cuts_two_generators() {
    // The plane x = z/4 - 1/2 passes through the apex (0, 0, 2) at a
    // slant, cutting the cone in two generators off the axis plane.
    let n = Vec3::new(1.0, 0.0, -0.25);
    let n = n * (1.0 / n.length());
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 2.0),
        normal: n,
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: n.cross(&Vec3::new(0.0, 1.0, 0.0)),
    };
    assert_generators(&plane);
}

#[test]
fn a_plane_tangent_along_a_generator_cuts_nothing() {
    // The cone's half-angle plane: normal at the half angle off the
    // radial direction touches the cone along one generator.
    let a = 0.5_f64.atan();
    let n = Vec3::new(a.cos(), 0.0, a.sin());
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 2.0),
        normal: n,
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: n.cross(&Vec3::new(0.0, 1.0, 0.0)),
    };
    let tol = Tolerance::default();
    let result = try_plane_cone(
        &plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        &cone(),
        (0.0, std::f64::consts::TAU),
        (0.0, 2.0),
        &tol,
    );
    assert!(result.is_some_and(|r| r.is_empty()));
}
