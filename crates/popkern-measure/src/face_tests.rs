//! Tests for face measurement functions.

use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn make_cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("cube")
}

#[test]
fn unit_cube_face_area_is_one() {
    let cube = make_cube();
    let tol = Tolerance::default();
    let face = cube.faces.keys().next().unwrap();
    let area = face_area(&cube, face, &tol).unwrap();
    assert!((area - 1.0).abs() < 1e-10, "expected 1.0, got {area}");
}

#[test]
fn unit_cube_face_perimeter_is_four() {
    let cube = make_cube();
    let face = cube.faces.keys().next().unwrap();
    let perim = face_perimeter(&cube, face).unwrap();
    assert!((perim - 4.0).abs() < 1e-10, "expected 4.0, got {perim}");
}

#[test]
fn unit_cube_face_normal_is_unit() {
    let cube = make_cube();
    let face = cube.faces.keys().next().unwrap();
    let n = face_normal(&cube, face).unwrap();
    assert!((n.length() - 1.0).abs() < 1e-10);
}

#[test]
fn planar_face_radius_is_incompatible() {
    let cube = make_cube();
    let face = cube.faces.keys().next().unwrap();
    assert!(matches!(
        face_radius(&cube, face),
        Err(MeasureError::IncompatibleSurface { .. })
    ));
}

#[test]
fn planar_face_curvature_is_zero() {
    let cube = make_cube();
    let face = cube.faces.keys().next().unwrap();
    let (k_min, k_max) = face_curvature(&cube, face).unwrap();
    assert!(k_min.abs() < 1e-10);
    assert!(k_max.abs() < 1e-10);
}

/// Find the face with the given surface type in a BRep.
fn find_face_by_surface<F>(brep: &BRep, pred: F) -> Option<popkern_core::topology::FaceHandle>
where
    F: Fn(&SurfaceGeometry) -> bool,
{
    brep.faces.iter().find_map(|(fh, fd)| {
        fd.surface
            .as_ref()
            .and_then(|s| if pred(&s.geometry) { Some(fh) } else { None })
    })
}

#[test]
fn cylinder_lateral_face_area_exact() {
    // r=1, h=2 → lateral area = 2π * r * h = 4π
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = find_face_by_surface(&cyl, |g| matches!(g, SurfaceGeometry::Cylinder { .. }))
        .expect("cylinder face");
    let area = face_area(&cyl, face, &tol).unwrap();
    let expected = std::f64::consts::TAU * 1.0 * 2.0; // 4π
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}

#[test]
fn cylinder_lateral_face_area_custom_radius() {
    // r=2.5, h=3.0 → lateral area = 2π * 2.5 * 3.0 = 15π
    let cyl = BRep::cylinder(2.5, 3.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = find_face_by_surface(&cyl, |g| matches!(g, SurfaceGeometry::Cylinder { .. }))
        .expect("cylinder face");
    let area = face_area(&cyl, face, &tol).unwrap();
    let expected = std::f64::consts::TAU * 2.5 * 3.0;
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}

#[test]
fn sphere_face_area_exact() {
    // r=1 → area = 4π
    let sph = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = find_face_by_surface(&sph, |g| matches!(g, SurfaceGeometry::Sphere { .. }))
        .expect("sphere face");
    let area = face_area(&sph, face, &tol).unwrap();
    let expected = 4.0 * std::f64::consts::PI;
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}

#[test]
fn sphere_face_area_custom_radius() {
    // r=2.5 → area = 4π * 2.5² = 25π
    let sph = BRep::sphere(2.5, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = find_face_by_surface(&sph, |g| matches!(g, SurfaceGeometry::Sphere { .. }))
        .expect("sphere face");
    let area = face_area(&sph, face, &tol).unwrap();
    let expected = 4.0 * std::f64::consts::PI * 2.5 * 2.5;
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}

#[test]
fn cone_lateral_face_area_exact() {
    // r=1, h=1 => lateral area = pi * r * sqrt(r^2 + h^2) = pi * sqrt(2)
    let cone = BRep::cone(1.0, 1.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = find_face_by_surface(&cone, |g| matches!(g, SurfaceGeometry::Cone { .. }))
        .expect("cone lateral face");
    let area = face_area(&cone, face, &tol).unwrap();
    let expected = std::f64::consts::PI * (2.0_f64).sqrt();
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}

#[test]
fn cone_lateral_face_area_custom_params() {
    // r=2, h=3 => lateral area = pi * r * sqrt(r^2 + h^2) = pi * 2 * sqrt(13)
    let cone = BRep::cone(2.0, 3.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = find_face_by_surface(&cone, |g| matches!(g, SurfaceGeometry::Cone { .. }))
        .expect("cone lateral face");
    let area = face_area(&cone, face, &tol).unwrap();
    let expected = std::f64::consts::PI * 2.0 * 13.0_f64.sqrt();
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}

#[test]
fn torus_rectangular_patch_area_exact() {
    // A rectangular torus patch at (u=[0,pi/2], v=[0,pi/2]) with R=3, r=1.
    // Exact: r * du * |R*dv + r*(sin(v1)-sin(v0))|
    //       = 1 * pi/2 * |3*(pi/2) + 1*(1-0)|
    //       = pi/2 * (3*pi/2 + 1)
    let tol = Tolerance::default();
    let torus = BRep::torus(3.0, 1.0, tol).unwrap();
    // Find a face with both delta_u and delta_v = pi/2 (rectangular patch).
    use crate::face_helpers::pcurve_uv_range;
    use std::f64::consts::PI;
    let half = PI / 2.0;
    // Find a face at (u=[0,pi/2], v=[0,pi/2]) — the first v-band.
    let face = torus
        .faces
        .iter()
        .find_map(|(fh, fd)| {
            fd.surface.as_ref().and_then(|s| {
                if !matches!(s.geometry, SurfaceGeometry::Torus { .. }) {
                    return None;
                }
                let ((u0, u1), (v0, v1)) = pcurve_uv_range(&torus, fh)?;
                let du = (u1 - u0).abs();
                let dv = (v1 - v0).abs();
                if (du - half).abs() < 1e-9 && (dv - half).abs() < 1e-9 && v0.abs() < 1e-9 {
                    Some(fh)
                } else {
                    None
                }
            })
        })
        .expect("rectangular torus patch at v=[0,pi/2]");
    let area = face_area(&torus, face, &tol).unwrap();
    // Expected: pi/2 * (3*pi/2 + 1)  (v0=0, v1=pi/2)
    let expected = (PI / 2.0) * (3.0 * PI / 2.0 + 1.0);
    assert!(
        (area - expected).abs() < 1e-10,
        "expected {expected}, got {area}"
    );
}
