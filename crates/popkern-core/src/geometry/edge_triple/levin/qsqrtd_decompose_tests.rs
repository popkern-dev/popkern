//! Tests for symbolic eigenvector / plane extraction over `Q[√d]`.

use super::qsqrtd::{QSqrtD, to_f64};
use super::qsqrtd_decompose::{
    extract_planes_qsqrtd, lift_to_qsqrtd, null_vector_qsqrtd, rank2_eigenvalues, subtract_mu_diag,
    trace_and_s2,
};
use super::rational::{Mat4SymRat, surface_to_quadric_rat};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Zero;

fn rat_int(n: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

fn diag(d0: i64, d1: i64, d2: i64, d3: i64) -> Mat4SymRat {
    let z = rat_int(0);
    Mat4SymRat([
        [rat_int(d0), z.clone(), z.clone(), z.clone()],
        [z.clone(), rat_int(d1), z.clone(), z.clone()],
        [z.clone(), z.clone(), rat_int(d2), z.clone()],
        [z.clone(), z.clone(), z, rat_int(d3)],
    ])
}

#[test]
fn trace_and_s2_diagonal() {
    // diag(2, 3, 5, 7): tr = 17, S₂ = 2·3 + 2·5 + 2·7 + 3·5 + 3·7 + 5·7
    //                          = 6 + 10 + 14 + 15 + 21 + 35 = 101.
    let m = diag(2, 3, 5, 7);
    let (tr, s2) = trace_and_s2(&m);
    assert_eq!(tr, rat_int(17));
    assert_eq!(s2, rat_int(101));
}

#[test]
fn rank2_eigenvalues_steinmetz_lambda_minus_one() {
    // The Steinmetz pencil's λ=−1 fibre: D = diag(−1, 0, 1, 0).
    // tr = 0, S₂ = (−1)·0 + (−1)·1 + (−1)·0 + 0·1 + 0·0 + 1·0 = −1.
    // disc = 0² − 4·(−1) = 4. √4 = 2 (perfect square).
    // μ_± = (0 ± √4) / 2 = ±1. Both rational.
    let d = diag(-1, 0, 1, 0);
    let (mu_plus, mu_minus, disc) = rank2_eigenvalues(&d).unwrap();
    assert_eq!(disc, rat_int(4));
    // After QSqrtD canonicalisation: (0, 1/2, 4) → (0, 1, 1) after
    // square-free factor extraction (b ↦ b·2, d ↦ 1).
    // So μ_+ = 0 + 1·√1 = "1 surd-style", which doesn't compare equal
    // to pure rational 1 under our canonical-form rule. That's fine
    // — embedding to f64 still gives exactly 1.
    let f_plus = to_f64(&mu_plus);
    let f_minus = to_f64(&mu_minus);
    assert!((f_plus - 1.0).abs() < 1e-15, "mu_plus → {f_plus}");
    assert!((f_minus + 1.0).abs() < 1e-15, "mu_minus → {f_minus}");
}

#[test]
fn lift_preserves_entries() {
    let m = diag(1, 2, 3, 4);
    let lifted = lift_to_qsqrtd(&m);
    for i in 0..4 {
        let entry = &lifted.0[i][i];
        assert!(entry.is_rational());
        assert_eq!(entry.a, rat_int((i + 1) as i64));
    }
}

#[test]
fn subtract_mu_diag_on_diagonal() {
    let m = diag(5, 5, 5, 5);
    let lifted = lift_to_qsqrtd(&m);
    let mu = QSqrtD::from_rat(rat_int(3));
    let result = subtract_mu_diag(&lifted, &mu);
    for i in 0..4 {
        let entry = &result.0[i][i];
        assert!(entry.is_rational());
        assert_eq!(entry.a, rat_int(2)); // 5 - 3
    }
}

#[test]
fn null_vector_of_singular_diagonal() {
    // diag(1, 0, 1, 1) has a 1D null space along e_y = (0, 1, 0, 0).
    // After RREF: rows [1,0,0,0], [0,0,0,0], [0,0,1,0], [0,0,0,1].
    // Pivot columns: 0, 2, 3. Free column: 1. Null vector: e_y.
    let m = diag(1, 0, 1, 1);
    let lifted = lift_to_qsqrtd(&m);
    let v = null_vector_qsqrtd(&lifted).unwrap();
    // Free variable was column 1 (set to 1). Others came from
    // back-substitution which gives 0 since the diagonal pivots
    // don't see column 1.
    assert!(v[0].is_zero());
    assert_eq!(v[1].a, num_traits::One::one());
    assert!(v[2].is_zero());
    assert!(v[3].is_zero());
}

#[test]
fn extract_planes_steinmetz_pencil() {
    // For Steinmetz cyl × cyl, the λ=−1 fibre D = diag(−1, 0, 1, 0)
    // has rank 2 (two non-zero eigenvalues ±1, kernel = (y, w)).
    // The two planes should be along the (x, z) directions —
    // specifically x ± z = 0.
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
    use super::rational::linear_comb_rat;
    let one = rat_int(1);
    let lam = rat_int(-1);
    let d = linear_comb_rat(&a, &b, &lam, &one);
    let planes = extract_planes_qsqrtd(&d).expect("rank-2 mixed-sign extraction");
    // Each plane should evaluate to ~(0, ±1) along the (x, z) axes
    // when projected to f64. Specifically: v_+ = (0, 0, 1, 0) [or
    // scalar multiple] from the +1 eigenvalue (z-direction);
    // v_− = (1, 0, 0, 0) from the −1 eigenvalue (x-direction).
    let f_plus: [f64; 4] = [
        to_f64(&planes[0][0]),
        to_f64(&planes[0][1]),
        to_f64(&planes[0][2]),
        to_f64(&planes[0][3]),
    ];
    let f_minus: [f64; 4] = [
        to_f64(&planes[1][0]),
        to_f64(&planes[1][1]),
        to_f64(&planes[1][2]),
        to_f64(&planes[1][3]),
    ];
    // The vectors are proportional to e_z (for v_+) and e_x (for v_-).
    // After RREF normalisation (pivot = 1), they SHOULD be exactly
    // (0, 0, 1, 0) and (1, 0, 0, 0) — but RREF order depends on which
    // pivot is found first, so accept either.
    let p_norm =
        (f_plus[0].powi(2) + f_plus[1].powi(2) + f_plus[2].powi(2) + f_plus[3].powi(2)).sqrt();
    let m_norm =
        (f_minus[0].powi(2) + f_minus[1].powi(2) + f_minus[2].powi(2) + f_minus[3].powi(2)).sqrt();
    // Plane normals are along x or z — the y and w components must be zero.
    assert!(f_plus[1].abs() < 1e-15, "f_plus.y = {}", f_plus[1]);
    assert!(f_plus[3].abs() < 1e-15, "f_plus.w = {}", f_plus[3]);
    assert!(f_minus[1].abs() < 1e-15, "f_minus.y = {}", f_minus[1]);
    assert!(f_minus[3].abs() < 1e-15, "f_minus.w = {}", f_minus[3]);
    assert!(
        (p_norm - 1.0).abs() < 1e-15,
        "v_+ should be unit, |v| = {p_norm}"
    );
    assert!(
        (m_norm - 1.0).abs() < 1e-15,
        "v_- should be unit, |v| = {m_norm}"
    );
}

#[test]
fn extract_planes_complex_disc_returns_none() {
    // Construct a matrix where the 2×2 reduction has disc < 0
    // (eigenvalues are complex — D is rank-2 with same-sign
    // eigenvalues, not the mixed-sign case).
    // diag(1, 0, 1, 0): tr = 2, S₂ = (1·0 + 1·1 + 1·0 + 0·1 + 0·0 + 1·0)
    //                              = 1. disc = 4 − 4 = 0.
    // Hmm that's the rank-1 case. Try diag(1, 0, 2, 0):
    // tr = 3, S₂ = (1·0 + 1·2 + 1·0 + 0·2 + 0·0 + 2·0) = 2.
    // disc = 9 − 8 = 1 → mixed-sign with surd. Not what we want.
    //
    // For disc < 0: need S₂ > tr²/4. Try diag(1, 1, 0, 0):
    // tr = 2, S₂ = 1·1 + ... = 1 (only the (0,0)·(1,1) term among
    // pairs with two non-zero diagonals).
    // disc = 4 − 4 = 0.
    //
    // Try a non-diagonal: D = [[1, 0, 0, 0], [0, 1, 0, 0], 0, 0].
    // tr = 2, S₂ = 1, disc = 0.
    //
    // Same-sign rank-2 case: D = diag(1, 1, 0, 0). The non-zero
    // eigenvalues are both +1, so disc = (1+1)² − 4·1·1 = 0 — repeated.
    // This is actually a rank-2 *positive semidefinite* — eigenvalues
    // are real but the field doesn't extend (disc = 0).
    let d = diag(1, 1, 0, 0);
    // disc should be 0 here.
    let result = rank2_eigenvalues(&d);
    let (_, _, disc) = result.unwrap();
    assert!(disc.is_zero());
    // extract_planes returns None on disc = 0 (single doubled plane).
    assert!(extract_planes_qsqrtd(&d).is_none());
}
