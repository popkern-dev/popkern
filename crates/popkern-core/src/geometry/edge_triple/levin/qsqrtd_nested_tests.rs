//! Stage 2b.43 tests: `Q[√d, √α]` arithmetic.

use num_rational::BigRational;
use num_traits::One;

use super::qsqrtd::QSqrtD;
use super::qsqrtd_nested::QSqrtDSqrtL;
use super::rational::Rat;

fn rat(n: i64) -> Rat {
    BigRational::from_integer(num_bigint::BigInt::from(n))
}

fn qsqrtd_rat(n: i64) -> QSqrtD {
    QSqrtD::from_rat(rat(n))
}

/// `(1 + 1·√α) + (2 + 3·√α) = 3 + 4·√α` — additive componentwise.
#[test]
fn add_componentwise() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(1), qsqrtd_rat(1), alpha.clone());
    let b = QSqrtDSqrtL::new(qsqrtd_rat(2), qsqrtd_rat(3), alpha.clone());
    let c = a + b;
    assert_eq!(c.p0, qsqrtd_rat(3));
    assert_eq!(c.p1, qsqrtd_rat(4));
    assert_eq!(c.alpha, alpha);
}

/// `(5 + 2·√α) − (3 + 4·√α) = 2 − 2·√α`.
#[test]
fn sub_componentwise() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(5), qsqrtd_rat(2), alpha.clone());
    let b = QSqrtDSqrtL::new(qsqrtd_rat(3), qsqrtd_rat(4), alpha.clone());
    let c = a - b;
    assert_eq!(c.p0, qsqrtd_rat(2));
    assert_eq!(c.p1, qsqrtd_rat(-2));
}

/// Negation flips both components.
#[test]
fn neg_flips_signs() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(5), qsqrtd_rat(-2), alpha.clone());
    let n = -a;
    assert_eq!(n.p0, qsqrtd_rat(-5));
    assert_eq!(n.p1, qsqrtd_rat(2));
}

/// `(√α)² = α`. The defining property of the second extension.
#[test]
fn sqrt_alpha_squared_is_alpha() {
    let alpha = qsqrtd_rat(7);
    let s = QSqrtDSqrtL::new(QSqrtD::zero(), qsqrtd_rat(1), alpha.clone());
    let s_sq = s.clone() * s;
    assert_eq!(s_sq.p0, alpha);
    assert_eq!(s_sq.p1, QSqrtD::zero());
}

/// `(2 + 3·√α)(2 − 3·√α) = 4 − 9α` — difference of squares pattern.
/// With `α = 7`, the result is `4 − 63 = −59`.
#[test]
fn mul_conjugate_pair() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(2), qsqrtd_rat(3), alpha.clone());
    let b = QSqrtDSqrtL::new(qsqrtd_rat(2), qsqrtd_rat(-3), alpha.clone());
    let prod = a * b;
    assert_eq!(prod.p0, qsqrtd_rat(-59));
    assert_eq!(prod.p1, QSqrtD::zero());
}

/// Division by self is multiplicative identity.
#[test]
fn div_self_is_one() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(2), qsqrtd_rat(3), alpha.clone());
    let one = a.clone() / a;
    assert_eq!(one.p0, qsqrtd_rat(1));
    assert_eq!(one.p1, QSqrtD::zero());
}

/// `1 / (a + b√α) = (a − b√α) / (a² − α·b²)` — verify directly.
/// With α = 2, a = 3, b = 1: 1/(3 + √2) = (3 − √2)/(9 − 2) = (3/7) − (1/7)√2.
#[test]
fn reciprocal_formula() {
    let alpha = qsqrtd_rat(2);
    let denom = QSqrtDSqrtL::new(qsqrtd_rat(3), qsqrtd_rat(1), alpha.clone());
    let one = QSqrtDSqrtL::from_rat(Rat::one(), alpha.clone());
    let recip = one / denom;
    assert_eq!(recip.p0, QSqrtD::from_rat(rat(3) / rat(7)));
    assert_eq!(recip.p1, QSqrtD::from_rat(-rat(1) / rat(7)));
}

/// Conjugate flips only the `√α` part.
#[test]
fn conjugate_flips_sqrt_alpha_only() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(5), qsqrtd_rat(2), alpha.clone());
    let conj = a.conjugate();
    assert_eq!(conj.p0, qsqrtd_rat(5));
    assert_eq!(conj.p1, qsqrtd_rat(-2));
}

/// Field norm: `(5 + 2√7)(5 − 2√7) = 25 − 28 = −3`. (Yes, norms can
/// be negative when `α` is positive — `Q[√d, √α]` is a field, not a
/// metric space.)
#[test]
fn norm_difference_of_squares_negative() {
    let alpha = qsqrtd_rat(7);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(5), qsqrtd_rat(2), alpha);
    assert_eq!(a.norm(), qsqrtd_rat(-3));
}

/// Nested radical with non-trivial base radicand: `α = (1 + √5)/2`,
/// the golden ratio. `(√α)² = α` must hold in `Q[√5, √α]`.
#[test]
fn golden_ratio_sqrt_squared() {
    let d = rat(5);
    let half = Rat::one() / rat(2);
    let alpha = QSqrtD::new(half.clone(), half, d);
    let s = QSqrtDSqrtL::new(QSqrtD::zero(), QSqrtD::from_rat(Rat::one()), alpha.clone());
    let s_sq = s.clone() * s;
    assert_eq!(s_sq.p0, alpha);
    assert_eq!(s_sq.p1, QSqrtD::zero());
}

/// Mixed-coefficient multiplication over `Q[√5, √α]`. With
/// `α = (1 + √5)/2` and `x = (1 + √α)`, compute `x · x = 1 + 2√α + α`
/// and verify both components.
#[test]
fn golden_ratio_one_plus_sqrt_alpha_squared() {
    let d = rat(5);
    let half = Rat::one() / rat(2);
    let alpha = QSqrtD::new(half.clone(), half, d);
    let x = QSqrtDSqrtL::new(
        QSqrtD::from_rat(Rat::one()),
        QSqrtD::from_rat(Rat::one()),
        alpha.clone(),
    );
    let x_sq = x.clone() * x;
    // Expected p0 = 1 + α = 1 + (1 + √5)/2 = (3 + √5)/2.
    let three_halves = rat(3) / rat(2);
    let expected_p0 = QSqrtD::new(three_halves, Rat::one() / rat(2), rat(5));
    assert_eq!(x_sq.p0, expected_p0);
    // Expected p1 = 2.
    assert_eq!(x_sq.p1, qsqrtd_rat(2));
}

/// Mismatched `α` panics — the type rejects mixed-extension arithmetic.
#[test]
#[should_panic(expected = "mismatched second radicand")]
fn mismatched_alpha_panics_on_add() {
    let alpha_a = qsqrtd_rat(7);
    let alpha_b = qsqrtd_rat(11);
    let a = QSqrtDSqrtL::new(qsqrtd_rat(1), qsqrtd_rat(0), alpha_a);
    let b = QSqrtDSqrtL::new(qsqrtd_rat(1), qsqrtd_rat(0), alpha_b);
    let _ = a + b;
}

/// Lifting from `Q[√d]` is identity-on-the-`p0`-axis.
#[test]
fn from_qsqrtd_keeps_p0() {
    let alpha = qsqrtd_rat(7);
    let q = QSqrtD::new(rat(2), rat(3), rat(5)); // 2 + 3√5
    let lifted = QSqrtDSqrtL::from_qsqrtd(q.clone(), alpha);
    assert_eq!(lifted.p0, q);
    assert_eq!(lifted.p1, QSqrtD::zero());
}

/// Zero detection.
#[test]
fn zero_is_zero() {
    let alpha = qsqrtd_rat(7);
    let z = QSqrtDSqrtL::zero(alpha);
    assert!(z.is_zero());
}
