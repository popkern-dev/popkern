//! Tests for `classify_degenerate`: rank classification + plane
//! extraction from rank-deficient quadrics.

use super::decompose::{DegenerateShape, classify_degenerate, plane_to_surface};
use super::pencil::find_degenerate_fiber;
use super::quadric::{Mat4Sym, surface_to_quadric};
use crate::geometry::{Point3, Vec3, surface::SurfaceGeometry};

const TOL: f64 = 1e-9;

fn vec_close(a: &[f64; 4], b: &[f64; 4]) -> bool {
    (0..4).all(|i| (a[i] - b[i]).abs() < TOL)
}

fn vec_close_or_neg(a: &[f64; 4], b: &[f64; 4]) -> bool {
    vec_close(a, b) || vec_close(a, &[-b[0], -b[1], -b[2], -b[3]])
}

#[test]
fn rank_one_quadric_classifies_as_single_plane() {
    // v = (1, 0, 0, -3): plane x = 3.
    // D = v vᵀ.
    let v = [1.0_f64, 0.0, 0.0, -3.0];
    let mut d = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            d[i][j] = v[i] * v[j];
        }
    }
    let shape = classify_degenerate(&Mat4Sym(d), 1.0);
    let plane = match shape {
        DegenerateShape::SinglePlane(p) => p,
        other => panic!("expected SinglePlane, got {other:?}"),
    };
    // Normalised plane equation (n unit): (1, 0, 0) is already unit.
    assert!(vec_close_or_neg(&plane, &[1.0, 0.0, 0.0, -3.0]));
}

#[test]
fn rank_two_mixed_signs_extracts_two_planes() {
    // Construct D = a aᵀ − b bᵀ with two known unit-normal planes.
    // Plane 1: x = 0 → a_form = (1, 0, 0, 0)
    // Plane 2: y = 0 → b_form = (0, 1, 0, 0)
    // D = (1,0,0,0)(1,0,0,0)ᵀ − (0,1,0,0)(0,1,0,0)ᵀ
    //   = diag(1, -1, 0, 0)
    // From the factorisation D = ½((a+b)(a−b)ᵀ + (a−b)(a+b)ᵀ):
    //   plane 1 = a + b = (1, 1, 0, 0) — i.e. x + y = 0
    //   plane 2 = a − b = (1, -1, 0, 0) — i.e. x − y = 0
    let mut d = [[0.0_f64; 4]; 4];
    d[0][0] = 1.0;
    d[1][1] = -1.0;
    let shape = classify_degenerate(&Mat4Sym(d), 1.0);
    let (p1, p2) = match shape {
        DegenerateShape::TwoPlanes(p1, p2) => (p1, p2),
        other => panic!("expected TwoPlanes, got {other:?}"),
    };
    let inv_root2 = 1.0 / 2.0_f64.sqrt();
    let expected_a = [inv_root2, inv_root2, 0.0, 0.0];
    let expected_b = [inv_root2, -inv_root2, 0.0, 0.0];
    let matches_either = (vec_close_or_neg(&p1, &expected_a) && vec_close_or_neg(&p2, &expected_b))
        || (vec_close_or_neg(&p1, &expected_b) && vec_close_or_neg(&p2, &expected_a));
    assert!(
        matches_either,
        "got p1={p1:?}, p2={p2:?}; expected normalised (1, ±1, 0, 0)/√2"
    );
}

#[test]
fn rank_two_same_signs_classified_semidefinite() {
    // D = aaᵀ + bbᵀ — positive semidefinite, no real plane factorisation.
    let mut d = [[0.0_f64; 4]; 4];
    d[0][0] = 1.0;
    d[1][1] = 1.0;
    let shape = classify_degenerate(&Mat4Sym(d), 1.0);
    assert!(matches!(shape, DegenerateShape::Rank2SemiDefinite));
}

#[test]
fn rank_three_classified_as_cone() {
    // D = diag(1, 1, -1, 0) — rank 3, mixed-sign cone.
    let mut d = [[0.0_f64; 4]; 4];
    d[0][0] = 1.0;
    d[1][1] = 1.0;
    d[2][2] = -1.0;
    let shape = classify_degenerate(&Mat4Sym(d), 1.0);
    assert!(matches!(shape, DegenerateShape::Rank3Cone { .. }));
}

#[test]
fn zero_matrix_classified_as_zero() {
    let d = Mat4Sym([[0.0; 4]; 4]);
    let shape = classify_degenerate(&d, 1.0);
    assert!(matches!(shape, DegenerateShape::Zero));
}

#[test]
fn full_rank_classified_non_degenerate() {
    let d = Mat4Sym([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, -1.0],
    ]);
    let shape = classify_degenerate(&d, 1.0);
    assert!(matches!(shape, DegenerateShape::NonDegenerate));
}

#[test]
fn steinmetz_pencil_yields_two_planes_x_eq_pm_z() {
    // x² + y² = 1 (cylinder along z) and y² + z² = 1 (cylinder along x).
    // Their intersection lies in the planes x = ±z. Levin should
    // find the rank-2 plane-pair fiber.
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
    let fiber = find_degenerate_fiber(&a, &b).expect("Steinmetz fiber");
    let shape = classify_degenerate(&fiber.quadric, 1.0);
    let (p1, p2) = match shape {
        DegenerateShape::TwoPlanes(a, b) => (a, b),
        other => panic!("expected TwoPlanes, got {other:?}"),
    };
    // Either plane should match (1, 0, ±1, 0)/√2 (i.e., x ± z = 0).
    let inv_root2 = 1.0 / 2.0_f64.sqrt();
    let plus_z = [inv_root2, 0.0, inv_root2, 0.0];
    let minus_z = [inv_root2, 0.0, -inv_root2, 0.0];
    let matches_either = (vec_close_or_neg(&p1, &plus_z) && vec_close_or_neg(&p2, &minus_z))
        || (vec_close_or_neg(&p1, &minus_z) && vec_close_or_neg(&p2, &plus_z));
    assert!(
        matches_either,
        "Steinmetz planes p1={p1:?}, p2={p2:?}; expected x ± z = 0"
    );
}

#[test]
fn plane_to_surface_round_trip() {
    // x = 3 plane: (1, 0, 0, -3).
    let plane = plane_to_surface([1.0, 0.0, 0.0, -3.0]);
    let SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis: _,
        v_axis: _,
    } = plane
    else {
        panic!("expected Plane, got {plane:?}");
    };
    assert!((origin.x - 3.0).abs() < TOL);
    assert!(origin.y.abs() < TOL);
    assert!(origin.z.abs() < TOL);
    assert!((normal.x - 1.0).abs() < TOL);
}
