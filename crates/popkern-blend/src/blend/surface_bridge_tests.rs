//! Tests for the surface-aware bridge builder. Validates that the
//! resulting cubic Bezier:
//! 1. Interpolates `pa` at t=0 and `pb` at t=1.
//! 2. Samples in between lie on the supporting surface to a
//!    reasonable tolerance.
//! 3. Periodic unwrap chooses the short arc (cylinder, sphere,
//!    torus).

use super::build_straight_line_bridge;
use super::build_surface_bridge;
use super::point_on_surface_distance;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

fn plane_xy() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn cylinder_z(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

fn sphere_at_origin(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius,
        axis: popkern_core::geometry::Vec3::Z,
    }
}

fn cone_from_apex() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 2.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: std::f64::consts::FRAC_PI_6,
    }
}

fn torus_z(major: f64, minor: f64) -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: major,
        minor_radius: minor,
    }
}

fn sample_and_check_on_surface(
    curve: &popkern_core::geometry::NurbsCurve,
    surface: &SurfaceGeometry,
    tol: f64,
) {
    let (t0, t1) = curve.parameter_range();
    for i in 0..=8 {
        let t = t0 + (t1 - t0) * (i as f64) / 8.0;
        let p = curve.point_at(t);
        let d = point_on_surface_distance(&p, surface);
        assert!(
            d < tol,
            "bridge sample at t={t:.3} lies {d:.4e} from surface; tol={tol:e}"
        );
    }
}

#[test]
fn bridge_on_plane_is_straight_line() {
    let plane = plane_xy();
    let pa = Point3::new(0.0, 0.0, 0.0);
    let pb = Point3::new(1.0, 2.0, 0.0);
    let bridge = build_surface_bridge(&pa, &pb, &plane, 1e-9).expect("plane bridge");
    let (t0, t1) = bridge.parameter_range();
    let start = bridge.point_at(t0);
    let end = bridge.point_at(t1);
    assert!(pa.distance_to(&start) < 1e-10, "starts at pa");
    assert!(pb.distance_to(&end) < 1e-10, "ends at pb");
    sample_and_check_on_surface(&bridge, &plane, 1e-10);
}

#[test]
fn bridge_on_cylinder_hugs_surface() {
    // Two points on a cylinder of radius 1, separated by ~90°
    // around the axis — the bridge should follow the cylindrical
    // arc, not a straight chord.
    let cyl = cylinder_z(1.0);
    let pa = Point3::new(1.0, 0.0, 0.5);
    let pb = Point3::new(0.0, 1.0, 0.7);
    let bridge = build_surface_bridge(&pa, &pb, &cyl, 1e-9).expect("cylinder bridge");
    sample_and_check_on_surface(&bridge, &cyl, 1e-2);
}

#[test]
fn bridge_on_cylinder_picks_short_arc() {
    // Two points on the same circle, near-diametric opposites
    // in angle but close in (wrapped) UV: theta ~ -0.1 vs +0.1
    // via the +x direction. Short-arc unwrap must pick the 0.2-rad
    // span rather than the 6.1-rad wraparound.
    let cyl = cylinder_z(1.0);
    let eps = 0.1_f64;
    let pa = Point3::new(eps.cos(), -eps.sin(), 0.0);
    let pb = Point3::new(eps.cos(), eps.sin(), 0.0);
    let bridge = build_surface_bridge(&pa, &pb, &cyl, 1e-9).expect("cylinder bridge");
    sample_and_check_on_surface(&bridge, &cyl, 1e-2);
    // The bridge midpoint should be near (1, 0, 0) — the short-arc
    // midpoint. Going the long way would put it near (-1, 0, 0).
    let (t0, t1) = bridge.parameter_range();
    let mid = bridge.point_at(0.5 * (t0 + t1));
    assert!(
        mid.x > 0.99,
        "short-arc midpoint must be near +x side of the cylinder; got {mid:?}"
    );
}

#[test]
fn bridge_on_sphere_hugs_surface() {
    let sph = sphere_at_origin(2.0);
    // Two points on the sphere.
    let pa = Point3::new(2.0, 0.0, 0.0);
    let pb = Point3::new(0.0, 2.0, 0.0);
    let bridge = build_surface_bridge(&pa, &pb, &sph, 1e-9).expect("sphere bridge");
    // Tolerance is looser on sphere because UV-linear in the
    // spherical parameterisation isn't quite a great circle —
    // but it stays very close to the surface over a ~90° arc.
    sample_and_check_on_surface(&bridge, &sph, 5e-2);
}

#[test]
fn bridge_on_cone_hugs_surface() {
    let cone = cone_from_apex();
    // Two points on the cone at the same axial distance, about
    // 60° apart in angle.
    let half_angle = std::f64::consts::FRAC_PI_6;
    let v = 1.0; // axial distance from apex
    let r = v * half_angle.tan();
    let pa = Point3::new(r, 0.0, 2.0 - v);
    let pb = Point3::new(
        r * (std::f64::consts::FRAC_PI_3).cos(),
        r * (std::f64::consts::FRAC_PI_3).sin(),
        2.0 - v,
    );
    let bridge = build_surface_bridge(&pa, &pb, &cone, 1e-9).expect("cone bridge");
    sample_and_check_on_surface(&bridge, &cone, 5e-2);
}

#[test]
fn bridge_on_torus_hugs_surface() {
    let torus = torus_z(2.0, 0.5);
    // Two points on the torus's outer equator (minor angle = 0,
    // varying major angle).
    let pa = Point3::new(2.5, 0.0, 0.0);
    let pb = Point3::new(
        2.5 * (std::f64::consts::FRAC_PI_4).cos(),
        2.5 * (std::f64::consts::FRAC_PI_4).sin(),
        0.0,
    );
    let bridge = build_surface_bridge(&pa, &pb, &torus, 1e-9).expect("torus bridge");
    sample_and_check_on_surface(&bridge, &torus, 5e-2);
}

#[test]
fn bridge_off_surface_returns_none() {
    let plane = plane_xy();
    let pa = Point3::new(0.0, 0.0, 0.0);
    let pb = Point3::new(1.0, 0.0, 1.0); // 1 unit above the plane.
    let bridge = build_surface_bridge(&pa, &pb, &plane, 1e-3);
    assert!(bridge.is_none(), "off-surface endpoint should return None");
}

#[test]
fn straight_line_fallback_constructs_valid_bridge() {
    let pa = Point3::new(0.0, 0.0, 0.0);
    let pb = Point3::new(1.0, 2.0, 3.0);
    let bridge = build_straight_line_bridge(&pa, &pb).expect("straight-line fallback");
    let (t0, t1) = bridge.parameter_range();
    assert!(pa.distance_to(&bridge.point_at(t0)) < 1e-10);
    assert!(pb.distance_to(&bridge.point_at(t1)) < 1e-10);
}

#[test]
fn point_on_surface_distance_plane() {
    let plane = plane_xy();
    assert!(point_on_surface_distance(&Point3::new(1.0, 2.0, 0.0), &plane) < 1e-10);
    assert!((point_on_surface_distance(&Point3::new(1.0, 2.0, 3.0), &plane) - 3.0).abs() < 1e-10);
}

#[test]
fn point_on_surface_distance_cylinder() {
    let cyl = cylinder_z(1.0);
    // Point on cylinder.
    assert!(point_on_surface_distance(&Point3::new(1.0, 0.0, 0.5), &cyl) < 1e-10);
    // Point outside cylinder (radial = 2 vs radius 1).
    let d = point_on_surface_distance(&Point3::new(2.0, 0.0, 0.5), &cyl);
    assert!((d - 1.0).abs() < 1e-10);
}
