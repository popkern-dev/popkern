//! Stage 2b.40 tests: symbolic tangent-intersection detection.
//!
//! Cases:
//! - Cylinder with inscribed sphere — tangent along the equator
//!   (canonical 1-D-tangent case; rank-1 fibre at λ=−1).
//! - Two transverse spheres — intersection circle, not tangent.
//! - Two tangent spheres (kissing at a point) — *not* a tangent-curve
//!   case (intersection is 0-D), so no rank-1 fibre exists.
//! - Plane input directly — rank-1 by construction; helper recovers
//!   the plane's own normal.

use super::rational::{Mat4SymRat, Rat, surface_to_quadric_rat};
use super::tangent::{extract_plane_from_rank1, is_rank1_or_less, surfaces_are_tangent};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use num_rational::BigRational;
use num_traits::{One, Zero};

fn rat(n: i64) -> Rat {
    BigRational::from_integer(num_bigint::BigInt::from(n))
}

fn unit_z_cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    }
}

fn unit_sphere_at_origin() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    }
}

fn unit_sphere_at(cx: f64, cy: f64, cz: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(cx, cy, cz),
        radius: 1.0,
    }
}

fn z0_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

#[test]
fn rank1_diagonal_matrix_detected() {
    // diag(0,0,1,0) — only nonzero entry on the diagonal at index 2.
    // Every 2×2 minor is zero (column/row index 2 is the only one
    // with non-zero entries, and any 2-row × 2-col minor must avoid
    // duplicates so it always picks at least one zero row/col).
    let mut m = Mat4SymRat::zero();
    m.0[2][2] = rat(1);
    assert!(is_rank1_or_less(&m));
    let n = extract_plane_from_rank1(&m).expect("plane normal");
    assert_eq!(n, [Rat::zero(), Rat::zero(), rat(1), Rat::zero()]);
}

#[test]
fn rank2_diagonal_matrix_not_rank1() {
    // diag(1,1,0,0) — rank-2.
    let mut m = Mat4SymRat::zero();
    m.0[0][0] = rat(1);
    m.0[1][1] = rat(1);
    assert!(!is_rank1_or_less(&m));
}

#[test]
fn outer_product_matrix_is_rank1() {
    // n = [1, 2, 3, -1]; D = n·nᵀ.
    let n = [rat(1), rat(2), rat(3), rat(-1)];
    let mut d = Mat4SymRat::zero();
    for i in 0..4 {
        for j in 0..4 {
            d.0[i][j] = &n[i] * &n[j];
        }
    }
    assert!(is_rank1_or_less(&d));
    let recovered = extract_plane_from_rank1(&d).expect("plane normal");
    // Recovered normal is parallel to n. The pivot column is the
    // diagonal index with the largest |d_jj| = n_j² — index 2 here
    // (n_2 = 3 ⇒ d_22 = 9). column[2] = n_2·n = 3·n.
    assert_eq!(recovered, [rat(3), rat(6), rat(9), rat(-3)]);
}

#[test]
fn cylinder_with_inscribed_sphere_is_tangent_at_equator() {
    let cyl = unit_z_cylinder();
    let sphere = unit_sphere_at_origin();
    let plane = surfaces_are_tangent(&cyl, &sphere).expect("tangent plane");
    // Plane normal proportional to [0, 0, 1, 0] (the z=0 plane).
    // The fibre at λ=−1 gives D = diag(0,0,1,0) directly, so the
    // recovered column is exactly [0, 0, 1, 0].
    assert_eq!(plane[0], Rat::zero());
    assert_eq!(plane[1], Rat::zero());
    assert!(!plane[2].is_zero());
    assert_eq!(plane[3], Rat::zero());
}

#[test]
fn two_transverse_spheres_have_no_tangent_plane() {
    // Spheres at (0,0,0) and (1,0,0), both radius 1. Distance = 1 <
    // sum-of-radii = 2 — they intersect along a circle (not tangent).
    let a = unit_sphere_at_origin();
    let b = unit_sphere_at(1.0, 0.0, 0.0);
    assert!(surfaces_are_tangent(&a, &b).is_none());
}

#[test]
fn two_kissing_spheres_have_no_tangent_plane_in_levin_sense() {
    // Spheres at (0,0,0) and (2,0,0), both radius 1 — kiss at
    // (1,0,0). Intersection is a single point (0-D). Levin's
    // tangent-along-curve detection should NOT classify this as
    // tangent: rational roots {-1, +1}, both fibres rank-2.
    let a = unit_sphere_at_origin();
    let b = unit_sphere_at(2.0, 0.0, 0.0);
    assert!(surfaces_are_tangent(&a, &b).is_none());
}

#[test]
fn plane_quadric_is_rank1_by_construction() {
    let p = z0_plane();
    let q = surface_to_quadric_rat(&p).expect("plane is rational");
    assert!(is_rank1_or_less(&q));
    let n = extract_plane_from_rank1(&q).expect("plane normal");
    // The plane z=0 has 4-vector (n, -n·o) = ([0,0,1], 0). The
    // quadric is q = v·vᵀ for v = [0,0,1,0]. Diagonal: only
    // d_22 = 1 ≠ 0; column 2 = [0, 0, 1, 0].
    assert_eq!(n, [Rat::zero(), Rat::zero(), Rat::one(), Rat::zero()]);
}

#[test]
fn cylinder_alone_is_not_rank1() {
    // cylinder quadric is rank-3 (2 nonzero radial eigenvalues, plus
    // a homogeneous-coord eigenvalue from the −r² shift). is_rank1
    // must return false.
    let q = surface_to_quadric_rat(&unit_z_cylinder()).expect("cyl is rational");
    assert!(!is_rank1_or_less(&q));
}

#[test]
fn non_quadric_surface_returns_none() {
    // Torus is *not* a quadric (degree-4 surface). The rational
    // pipeline rejects it at `surface_to_quadric_rat`, so
    // `surfaces_are_tangent` must return `None` rather than mis-using
    // a partial quadric form.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let sphere = unit_sphere_at_origin();
    assert!(surfaces_are_tangent(&torus, &sphere).is_none());
}
