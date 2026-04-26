//! Tests for the `Q[√d]` field type.

use super::qsqrtd::{QSqrtD, to_f64};
use super::rational::Rat;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

fn rat(n: i64, d: i64) -> BigRational {
    BigRational::new(BigInt::from(n), BigInt::from(d))
}

fn rat_int(n: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

fn qs(a: i64, b: i64, d: i64) -> QSqrtD {
    QSqrtD::new(rat_int(a), rat_int(b), rat_int(d))
}

#[test]
fn rational_only_is_canonical() {
    let v = QSqrtD::new(rat_int(3), rat_int(0), rat_int(7));
    assert!(v.is_rational());
    assert_eq!(v.b, Rat::zero());
    assert_eq!(v.d, Rat::zero());
}

#[test]
fn zero_radicand_collapses_to_rational() {
    let v = QSqrtD::new(rat_int(5), rat_int(7), rat_int(0));
    assert!(v.is_rational());
    assert_eq!(v.a, rat_int(5));
}

#[test]
fn perfect_square_d_collapses_to_pure_rational() {
    // 1 + 1·√4. Square-free reduction extracts √4 = 2·√1, so the
    // intermediate form is (1, 2, 1). Then d = 1 collapses (√1 = 1)
    // to pure rational 1 + 2·1 = 3.
    let v = QSqrtD::new(rat_int(1), rat_int(1), rat_int(4));
    assert!(v.is_rational());
    assert_eq!(v.a, rat_int(3));
    assert!(v.b.is_zero());
    assert!(v.d.is_zero());
}

#[test]
fn d_eight_factors_to_d_two() {
    // 1 + 1·√8 → 1 + 2·√2 (since 8 = 4·2).
    let v = QSqrtD::new(rat_int(1), rat_int(1), rat_int(8));
    assert_eq!(v.b, rat_int(2));
    assert_eq!(v.d, rat_int(2));
}

#[test]
fn d_50_factors_to_d_2() {
    // 50 = 25·2 ⇒ √50 = 5·√2.
    let v = QSqrtD::new(rat_int(0), rat_int(1), rat_int(50));
    assert_eq!(v.b, rat_int(5));
    assert_eq!(v.d, rat_int(2));
}

#[test]
fn add_with_matching_radicand() {
    // (1 + 2√3) + (3 + 4√3) = 4 + 6√3.
    let lhs = qs(1, 2, 3);
    let rhs = qs(3, 4, 3);
    let r = lhs + rhs;
    assert_eq!(r.a, rat_int(4));
    assert_eq!(r.b, rat_int(6));
    assert_eq!(r.d, rat_int(3));
}

#[test]
fn add_pure_rational_with_surd() {
    // 5 + (1 + √2) = 6 + √2.
    let lhs = QSqrtD::from_rat(rat_int(5));
    let rhs = qs(1, 1, 2);
    let r = lhs + rhs;
    assert_eq!(r.a, rat_int(6));
    assert_eq!(r.b, rat_int(1));
    assert_eq!(r.d, rat_int(2));
}

#[test]
fn mul_classical_identity() {
    // (1 + √2)·(1 − √2) = 1 − 2 = −1.
    let lhs = qs(1, 1, 2);
    let rhs = qs(1, -1, 2);
    let r = lhs * rhs;
    assert_eq!(r.a, rat_int(-1));
    assert!(r.b.is_zero());
    assert!(r.is_rational());
}

#[test]
fn mul_general() {
    // (2 + 3√5)·(4 + 7√5) = 8 + 14√5 + 12√5 + 21·5 = (8+105) + (14+12)√5
    //                    = 113 + 26√5.
    let lhs = qs(2, 3, 5);
    let rhs = qs(4, 7, 5);
    let r = lhs * rhs;
    assert_eq!(r.a, rat_int(113));
    assert_eq!(r.b, rat_int(26));
    assert_eq!(r.d, rat_int(5));
}

#[test]
fn norm_is_rational() {
    // norm(1 + √7) = 1 − 7 = −6.
    let v = qs(1, 1, 7);
    assert_eq!(v.norm(), rat_int(-6));
}

#[test]
fn div_inverse() {
    // (1 + √2)·1/(1 + √2) = 1.
    // 1/(1 + √2) = (1 − √2) / (1 − 2) = (1 − √2) / −1 = −1 + √2.
    let v = qs(1, 1, 2);
    let one = QSqrtD::from_rat(Rat::one());
    let inv = one.clone() / v.clone();
    assert_eq!(inv.a, rat_int(-1));
    assert_eq!(inv.b, rat_int(1));
    let prod = v * inv;
    assert_eq!(prod, one);
}

#[test]
fn div_general() {
    // (3 + √5) / (1 + √5) = (3 + √5)·(1 − √5) / (1 − 5)
    //                    = (3 − 3√5 + √5 − 5) / −4
    //                    = (−2 − 2√5) / −4
    //                    = (1 + √5) / 2.
    let lhs = qs(3, 1, 5);
    let rhs = qs(1, 1, 5);
    let r = lhs / rhs;
    assert_eq!(r.a, rat(1, 2));
    assert_eq!(r.b, rat(1, 2));
}

#[test]
fn equality_canonical() {
    // 1 + 1·√4 == 1 + 2·1 == 3 — perfect-square radicand collapses
    // through canonicalisation, so two values that *should* be equal
    // (one written symbolically, one written as pure rational) DO
    // compare equal.
    let v1 = QSqrtD::new(rat_int(1), rat_int(1), rat_int(4));
    let v2 = QSqrtD::from_rat(rat_int(3));
    assert_eq!(v1, v2);

    // Surds with non-collapsing d compare equal iff (a, b, d) match.
    let v3 = QSqrtD::new(rat_int(1), rat_int(2), rat_int(3));
    let v4 = QSqrtD::new(rat_int(1), rat_int(2), rat_int(3));
    assert_eq!(v3, v4);
    let v5 = QSqrtD::new(rat_int(1), rat_int(2), rat_int(5));
    assert_ne!(v3, v5);
}

#[test]
fn negate() {
    let v = qs(2, -3, 7);
    let neg = -v;
    assert_eq!(neg.a, rat_int(-2));
    assert_eq!(neg.b, rat_int(3));
    assert_eq!(neg.d, rat_int(7));
}

#[test]
fn to_f64_matches_real_value() {
    // 1 + √2 ≈ 1 + 1.414213562 = 2.414213562.
    let v = qs(1, 1, 2);
    let f = to_f64(&v);
    assert!((f - (1.0 + 2.0_f64.sqrt())).abs() < 1e-15);
}

#[test]
fn ring_axiom_associativity() {
    // (a + b)·c = a·c + b·c (distributivity check, sanity).
    let a = qs(1, 2, 3);
    let b = qs(4, -1, 3);
    let c = qs(2, 5, 3);
    let lhs = (a.clone() + b.clone()) * c.clone();
    let rhs = a * c.clone() + b * c;
    assert_eq!(lhs, rhs);
}

#[test]
fn discriminant_for_known_2x2_rank2() {
    // For Steinmetz: cyl_z and cyl_x give D at λ=−1 = diag(−1, 0, 1, 0).
    // tr = 0, det of 2×2 sub-block = (−1)·(1) = −1.
    // disc = tr² − 4·det = 0 − 4·(−1) = 4 → √4 = 2 (perfect square).
    // μ_± = (0 ± 2) / 2 = ±1. The Q[√4] adjunction collapses through
    // the d=1 reduction: Q[√4] = Q. So a value `0 + 1·√4` in the
    // raw form becomes `2 + 0·√0` in canonical form.
    let disc = QSqrtD::new(rat_int(0), rat_int(1), rat_int(4));
    assert!(disc.is_rational());
    assert_eq!(disc.a, rat_int(2));
}
