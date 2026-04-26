//! Verify each surface variant produces a symmetric quadric form
//! that vanishes on sample surface points.

use super::quadric::{Mat4Sym, surface_to_quadric};
use crate::geometry::{Point3, Vec3, surface::SurfaceGeometry};

const TOL: f64 = 1e-12;

fn assert_vanishes_on(q: &Mat4Sym, p: Point3, ctx: &str) {
    let v = q.eval_affine(p.x, p.y, p.z);
    assert!(v.abs() < TOL, "{ctx}: pᵀAp = {v} at {p:?} — expected 0",);
}

#[test]
fn plane_quadric_symmetric_and_vanishing() {
    // xy-plane through origin.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let q = surface_to_quadric(&plane).expect("plane → quadric");
    assert!(q.symmetry_residual() < TOL);
    // Points on xy-plane (z=0) should vanish.
    assert_vanishes_on(&q, Point3::new(0.0, 0.0, 0.0), "plane origin");
    assert_vanishes_on(&q, Point3::new(2.5, -1.7, 0.0), "plane interior");
    // Off-plane point should NOT vanish.
    let off = q.eval_affine(0.0, 0.0, 1.0);
    assert!(
        (off - 1.0).abs() < TOL,
        "plane z=1: expected n²=1, got {off}"
    );
}

#[test]
fn plane_quadric_with_offset_origin() {
    // z = 5 plane.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 5.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let q = surface_to_quadric(&plane).unwrap();
    assert_vanishes_on(&q, Point3::new(3.0, -2.0, 5.0), "z=5 plane");
    // Off-plane: pᵀAp = (z − 5)² = 4 at z=7.
    let v = q.eval_affine(0.0, 0.0, 7.0);
    assert!((v - 4.0).abs() < TOL);
}

#[test]
fn sphere_quadric_symmetric_and_vanishing() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 2.0, 3.0),
        radius: 5.0,
    };
    let q = surface_to_quadric(&sphere).unwrap();
    assert!(q.symmetry_residual() < TOL);
    // Point on sphere: (1+5, 2, 3) is distance 5 from centre.
    assert_vanishes_on(&q, Point3::new(6.0, 2.0, 3.0), "sphere +x");
    assert_vanishes_on(&q, Point3::new(1.0, -3.0, 3.0), "sphere -y");
    assert_vanishes_on(&q, Point3::new(1.0, 2.0, 8.0), "sphere +z");
    // Centre point: pᵀAp = −ρ² = −25.
    let v = q.eval_affine(1.0, 2.0, 3.0);
    assert!((v + 25.0).abs() < TOL);
}

#[test]
fn cylinder_quadric_symmetric_and_vanishing() {
    // z-axis cylinder, origin (1, 0, 0), radius 2.
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    let q = surface_to_quadric(&cyl).unwrap();
    assert!(q.symmetry_residual() < TOL);
    // Surface points satisfy (x-1)² + y² = 4 (any z).
    assert_vanishes_on(&q, Point3::new(3.0, 0.0, 0.0), "cyl +x");
    assert_vanishes_on(&q, Point3::new(-1.0, 0.0, 7.5), "cyl −x at z=7.5");
    assert_vanishes_on(&q, Point3::new(1.0, 2.0, -3.0), "cyl +y at z=−3");
    // Axis point: pᵀAp = −r² = −4.
    let v = q.eval_affine(1.0, 0.0, 4.5);
    assert!((v + 4.0).abs() < TOL);
}

#[test]
fn cylinder_quadric_oblique_axis() {
    // Axis along (1, 1, 1)/√3, origin at origin, radius 1.
    let inv = 1.0 / 3.0_f64.sqrt();
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(inv, inv, inv),
        radius: 1.0,
    };
    let q = surface_to_quadric(&cyl).unwrap();
    assert!(q.symmetry_residual() < TOL);
    // Pick a point at unit perpendicular distance from the axis line
    // through origin in direction (1,1,1)/√3. Vector (1, -1, 0)/√2 is
    // perpendicular to (1,1,1)/√3 and has unit length.
    let inv2 = 1.0 / 2.0_f64.sqrt();
    assert_vanishes_on(&q, Point3::new(inv2, -inv2, 0.0), "oblique cyl r=1");
    // Move along axis: still on cylinder.
    let p_along = Point3::new(inv2 + inv, -inv2 + inv, inv);
    assert_vanishes_on(&q, p_along, "oblique cyl + axis step");
}

#[test]
fn cone_quadric_symmetric_and_vanishing() {
    // z-axis cone, apex at origin, half-angle 45° (so tan θ = 1).
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    let q = surface_to_quadric(&cone).unwrap();
    assert!(q.symmetry_residual() < TOL);
    // Surface: r = z·tan θ = z. Points (x, y, z) with x² + y² = z².
    assert_vanishes_on(&q, Point3::new(3.0, 0.0, 3.0), "45° cone +x");
    assert_vanishes_on(&q, Point3::new(0.0, -2.0, 2.0), "45° cone −y");
    assert_vanishes_on(&q, Point3::new(1.0, 1.0, 2.0_f64.sqrt()), "45° cone diag");
    // Apex is on the surface.
    assert_vanishes_on(&q, Point3::new(0.0, 0.0, 0.0), "45° cone apex");
}

#[test]
fn cone_quadric_offset_apex_oblique_axis() {
    // Apex at (1, 1, 0), axis along (1, 0, 1)/√2, half-angle = arctan(0.5).
    let inv = 1.0 / 2.0_f64.sqrt();
    let theta = (0.5_f64).atan();
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(1.0, 1.0, 0.0),
        axis: Vec3::new(inv, 0.0, inv),
        half_angle: theta,
    };
    let q = surface_to_quadric(&cone).unwrap();
    assert!(q.symmetry_residual() < TOL);
    // Apex must be on the surface.
    assert_vanishes_on(&q, Point3::new(1.0, 1.0, 0.0), "oblique cone apex");
    // Walk distance d=√2 along axis → axial coord √2 from apex; cone
    // radius at this axial offset = √2 · tan θ = √2/2.
    // Axis direction (inv, 0, inv); apex + √2·axis = (2, 1, 1).
    // Perp direction: y is perp to axis (axis has 0 y-component).
    // So (2, 1 + √2/2, 1) is on the surface.
    let half_root2 = 2.0_f64.sqrt() / 2.0;
    assert_vanishes_on(
        &q,
        Point3::new(2.0, 1.0 + half_root2, 1.0),
        "oblique cone +y",
    );
    assert_vanishes_on(
        &q,
        Point3::new(2.0, 1.0 - half_root2, 1.0),
        "oblique cone −y",
    );
}

#[test]
fn torus_returns_none() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    assert!(surface_to_quadric(&torus).is_none());
}
