//! Tests for pencil operations: characteristic polynomial,
//! degenerate-fiber search.

use super::pencil::{char_poly_in_lambda, det4, find_degenerate_fiber, linear_comb};
use super::quadric::{Mat4Sym, surface_to_quadric};
use crate::geometry::{Point3, Vec3, surface::SurfaceGeometry};

const TOL: f64 = 1e-9;

#[test]
fn det_of_diagonal_matrix() {
    let m = [
        [2.0, 0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0, 0.0],
        [0.0, 0.0, 5.0, 0.0],
        [0.0, 0.0, 0.0, 7.0],
    ];
    assert!((det4(&m) - 210.0).abs() < TOL);
}

#[test]
fn det_of_upper_triangular() {
    let m = [
        [1.0, 4.0, 5.0, 6.0],
        [0.0, 2.0, 3.0, 1.0],
        [0.0, 0.0, 4.0, 2.0],
        [0.0, 0.0, 0.0, 3.0],
    ];
    // Triangular det = product of diagonal.
    assert!((det4(&m) - 24.0).abs() < TOL);
}

#[test]
fn linear_comb_is_linear() {
    let a = Mat4Sym([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 3.0, 0.0],
        [0.0, 0.0, 0.0, 4.0],
    ]);
    let b = Mat4Sym([
        [0.5, 0.0, 0.0, 0.0],
        [0.0, 1.5, 0.0, 0.0],
        [0.0, 0.0, 2.5, 0.0],
        [0.0, 0.0, 0.0, 3.5],
    ]);
    let m = linear_comb(&a, &b, 2.0, 3.0);
    assert!((m.0[0][0] - (2.0 * 1.0 + 3.0 * 0.5)).abs() < TOL);
    assert!((m.0[3][3] - (2.0 * 4.0 + 3.0 * 3.5)).abs() < TOL);
}

#[test]
fn char_poly_recovers_known_polynomial() {
    // A = diag(1, 1, 1, -1), B = diag(1, 1, 1, -4).
    // λA + B = diag(λ+1, λ+1, λ+1, -λ-4).
    // det = (λ+1)³ · (-λ-4) = -(λ+1)³(λ+4)
    //     = -(λ⁴ + 7λ³ + 15λ² + 13λ + 4)
    //     = -λ⁴ - 7λ³ - 15λ² - 13λ - 4
    let s1 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let s2 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
    };
    let a = surface_to_quadric(&s1).unwrap();
    let b = surface_to_quadric(&s2).unwrap();
    let coefs = char_poly_in_lambda(&a, &b);
    let expected = [-1.0, -7.0, -15.0, -13.0, -4.0];
    for i in 0..5 {
        assert!(
            (coefs[i] - expected[i]).abs() < TOL,
            "coef[{i}]: {} vs {}",
            coefs[i],
            expected[i]
        );
    }
}

#[test]
fn degenerate_fiber_concentric_spheres() {
    // Two concentric spheres ρ=1, ρ=2. λA + B is rank-deficient at
    // λ = -1 (the fiber collapses to the "spheres at radius 0" — a
    // doubled point at the centre, in projective form a rank-1 quadric).
    let s1 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let s2 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
    };
    let a = surface_to_quadric(&s1).unwrap();
    let b = surface_to_quadric(&s2).unwrap();
    let fiber = find_degenerate_fiber(&a, &b).expect("should find a fiber");
    let det = det4(&fiber.quadric.0);
    assert!(det.abs() < 1e-6, "fiber det should be ~0, got {det}");
}

#[test]
fn degenerate_fiber_disjoint_spheres() {
    // Two disjoint unit spheres. Their pencil includes the radical
    // plane (squared linear form), which is degenerate.
    let s1 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let s2 = SurfaceGeometry::Sphere {
        center: Point3::new(3.0, 0.0, 0.0),
        radius: 1.0,
    };
    let a = surface_to_quadric(&s1).unwrap();
    let b = surface_to_quadric(&s2).unwrap();
    let fiber = find_degenerate_fiber(&a, &b).expect("should find a fiber");
    let det = det4(&fiber.quadric.0);
    assert!(det.abs() < 1e-6);
}

#[test]
fn degenerate_fiber_two_orthogonal_cylinders() {
    // x² + y² = 1 and y² + z² = 1 (the Steinmetz solid surfaces).
    // Intersection is two ellipses lying in the planes x = ±z.
    // Pencil should contain a plane-pair degenerate fiber.
    let cyl_z = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let cyl_x = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 1.0,
    };
    let a = surface_to_quadric(&cyl_z).unwrap();
    let b = surface_to_quadric(&cyl_x).unwrap();
    let fiber = find_degenerate_fiber(&a, &b).expect("should find a fiber");
    let det = det4(&fiber.quadric.0);
    assert!(det.abs() < 1e-6, "Steinmetz fiber det={det}");
}

#[test]
fn degenerate_input_a_returned_as_fiber() {
    // A is a plane (rank 1, det = 0); B is a sphere. Should immediately
    // return A as the degenerate fiber.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let a = surface_to_quadric(&plane).unwrap();
    let b = surface_to_quadric(&sphere).unwrap();
    let fiber = find_degenerate_fiber(&a, &b).expect("plane is degenerate");
    let det = det4(&fiber.quadric.0);
    assert!(det.abs() < 1e-6);
}
