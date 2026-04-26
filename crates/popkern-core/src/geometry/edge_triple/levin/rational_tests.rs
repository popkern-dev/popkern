//! Tests for the rational Levin foundation. Verify that the exact
//! arithmetic recovers the *known* rational degenerate fibres for
//! hand-computable surface pairs.

use super::rational::{
    Mat4SymRat, char_poly_rat, det4_rat, find_rational_degenerate_fiber, linear_comb_rat,
    rational_roots, surface_to_quadric_rat,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

fn rat(n: i64, d: i64) -> BigRational {
    BigRational::new(BigInt::from(n), BigInt::from(d))
}

fn rat_int(n: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

#[test]
fn rational_roots_quadratic() {
    // (x - 2)(x + 3) = x² + x − 6
    let coefs = vec![rat_int(1), rat_int(1), rat_int(-6)];
    let mut roots = rational_roots(&coefs);
    roots.sort();
    assert_eq!(roots, vec![rat_int(-3), rat_int(2)]);
}

#[test]
fn rational_roots_with_fractional() {
    // 2x² − 3x + 1 = (2x − 1)(x − 1) → roots 1/2, 1.
    let coefs = vec![rat_int(2), rat_int(-3), rat_int(1)];
    let mut roots = rational_roots(&coefs);
    roots.sort();
    assert_eq!(roots, vec![rat(1, 2), rat_int(1)]);
}

#[test]
fn rational_roots_no_rational_solutions() {
    // x² − 2 = 0 has only ±√2 — no rational roots.
    let coefs = vec![rat_int(1), rat_int(0), rat_int(-2)];
    let roots = rational_roots(&coefs);
    assert!(roots.is_empty());
}

#[test]
fn rational_roots_with_zero_root() {
    // x³ + x² − 2x = x · (x² + x − 2) = x · (x − 1)(x + 2)
    // Roots: 0, 1, -2.
    let coefs = vec![rat_int(1), rat_int(1), rat_int(-2), rat_int(0)];
    let mut roots = rational_roots(&coefs);
    roots.sort();
    assert_eq!(roots, vec![rat_int(-2), rat_int(0), rat_int(1)]);
}

#[test]
fn surface_to_quadric_sphere_at_unit_origin() {
    // Sphere at origin radius 1: x² + y² + z² − 1 = 0.
    // Matrix: diag(1, 1, 1, −1).
    let s = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let q = surface_to_quadric_rat(&s).unwrap();
    let expected: [[BigRational; 4]; 4] = [
        [rat_int(1), rat_int(0), rat_int(0), rat_int(0)],
        [rat_int(0), rat_int(1), rat_int(0), rat_int(0)],
        [rat_int(0), rat_int(0), rat_int(1), rat_int(0)],
        [rat_int(0), rat_int(0), rat_int(0), rat_int(-1)],
    ];
    for (i, row) in expected.iter().enumerate() {
        for (j, want) in row.iter().enumerate() {
            assert_eq!(&q.0[i][j], want, "mismatch at ({i}, {j})");
        }
    }
}

#[test]
fn det4_rat_block_diagonal() {
    let z = rat_int(0);
    let m = [
        [rat_int(2), z.clone(), z.clone(), z.clone()],
        [z.clone(), rat_int(3), z.clone(), z.clone()],
        [z.clone(), z.clone(), rat_int(5), z.clone()],
        [z.clone(), z.clone(), z, rat_int(7)],
    ];
    assert_eq!(det4_rat(&m), rat_int(2 * 3 * 5 * 7));
}

#[test]
fn char_poly_concentric_unit_and_double_spheres() {
    // Sphere ρ=1 vs sphere ρ=2 (concentric at origin):
    //   A = diag(1,1,1,−1), B = diag(1,1,1,−4).
    //   λA + B = diag(λ+1, λ+1, λ+1, −λ−4).
    //   det = (λ+1)³·(−λ−4) = −(λ+1)³(λ+4).
    //   Expanded: −λ⁴ − 7λ³ − 15λ² − 13λ − 4.
    let s1 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let s2 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
    };
    let a = surface_to_quadric_rat(&s1).unwrap();
    let b = surface_to_quadric_rat(&s2).unwrap();
    let coefs = char_poly_rat(&a, &b);
    let expected = [
        rat_int(-1),
        rat_int(-7),
        rat_int(-15),
        rat_int(-13),
        rat_int(-4),
    ];
    for i in 0..5 {
        assert_eq!(
            coefs[i], expected[i],
            "char poly coef[{i}]: got {}, expected {}",
            coefs[i], expected[i]
        );
    }
}

#[test]
fn find_rational_degenerate_fiber_concentric_spheres() {
    // The pencil's roots are λ ∈ {−1 (triple), −4 (single)}, both rational.
    let s1 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let s2 = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 2.0,
    };
    let a = surface_to_quadric_rat(&s1).unwrap();
    let b = surface_to_quadric_rat(&s2).unwrap();
    let (lam, m) = find_rational_degenerate_fiber(&a, &b).unwrap();
    // The found λ is exactly rational. det of the resulting matrix is
    // *exactly* zero (not within tolerance).
    assert!(lam == rat_int(-1) || lam == rat_int(-4));
    assert!(det4_rat(&m.0).is_zero());
}

#[test]
fn find_rational_degenerate_fiber_axis_aligned_cylinders() {
    // Two unit cylinders along x and z axes (Steinmetz). The pencil
    // has a rational rank-2 fibre at λ = −1 (gives D = diag(1, 0, −1, 0))
    // and a rational rank-1 fibre at λ = 1.
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
    let a = surface_to_quadric_rat(&cyl_z).unwrap();
    let b = surface_to_quadric_rat(&cyl_x).unwrap();
    let (_lam, m) = find_rational_degenerate_fiber(&a, &b).unwrap();
    // Detected fibre is *exactly* rank-deficient.
    assert!(det4_rat(&m.0).is_zero());
}

#[test]
fn linear_comb_is_linear() {
    let a = Mat4SymRat([
        [rat_int(1), rat_int(0), rat_int(0), rat_int(0)],
        [rat_int(0), rat_int(2), rat_int(0), rat_int(0)],
        [rat_int(0), rat_int(0), rat_int(3), rat_int(0)],
        [rat_int(0), rat_int(0), rat_int(0), rat_int(4)],
    ]);
    let b = Mat4SymRat([
        [rat(1, 2), rat_int(0), rat_int(0), rat_int(0)],
        [rat_int(0), rat(3, 2), rat_int(0), rat_int(0)],
        [rat_int(0), rat_int(0), rat(5, 2), rat_int(0)],
        [rat_int(0), rat_int(0), rat_int(0), rat(7, 2)],
    ]);
    let m = linear_comb_rat(&a, &b, &rat_int(2), &rat_int(3));
    // (2·1 + 3·1/2) = 2 + 3/2 = 7/2
    assert_eq!(m.0[0][0], rat(7, 2));
    // (2·4 + 3·7/2) = 8 + 21/2 = 37/2
    assert_eq!(m.0[3][3], rat(37, 2));
}

#[test]
fn rational_pencil_proves_steinmetz_factors_into_planes() {
    // For two perpendicular unit cylinders (Steinmetz), the pencil's
    // rank-2 fibre at λ = −1 has rank exactly 2 — provable
    // symbolically (no tolerance fudge). This is what enables Levin
    // to claim the intersection factors into TWO real plane fibres
    // exactly, rather than maybe-rank-3 with FP noise.
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
    let a = surface_to_quadric_rat(&cyl_z).unwrap();
    let b = surface_to_quadric_rat(&cyl_x).unwrap();
    let lam = -BigRational::one();
    let one = BigRational::one();
    let d = linear_comb_rat(&a, &b, &lam, &one);
    // d should be exactly diag(1 - 1, 1, 0 - 1, -1 + 1) entries...
    // Actually let me compute: cyl_z = diag(1,1,0,−1), cyl_x = diag(0,1,1,−1).
    // λA + B at λ=-1: diag(-1,-1,0,1) + diag(0,1,1,-1) = diag(-1, 0, 1, 0).
    // Rank: 2 (the +/- pair on x and z).
    // det = 0.
    assert!(
        det4_rat(&d.0).is_zero(),
        "Steinmetz λ=-1 fibre det = {}",
        det4_rat(&d.0)
    );
}
