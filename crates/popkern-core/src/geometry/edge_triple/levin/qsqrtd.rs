//! `Q[√d]` — rationals adjoined with one fixed algebraic square root.
//!
//! The Q foundation in `rational.rs` handles the rational pencil and
//! degenerate-fibre detection. But the **rank-2 plane decomposition**
//! step needs eigenvalues of a 2×2 sub-block, whose values are
//!
//! ```text
//! μ_± = (tr ± √(tr² − 4·det)) / 2
//! ```
//!
//! Generically `disc = tr² − 4·det` is rational but **not** a perfect
//! square — `√disc` is irrational. To represent the planes' normals
//! exactly, we extend the base field by a single algebraic surd: a
//! number `α = a + b·√d`, where `a, b ∈ Q` and `d ∈ Q` is fixed by
//! the discriminant of the case at hand.
//!
//! This module defines the `QSqrtD` type and its full field
//! arithmetic (`+`, `−`, `*`, `÷`, `==`). Operations between values
//! with **different** `d` are out of scope (would require `Q[√d₁,
//! √d₂]`); attempting them returns an error or panics, as flagged
//! per call site.
//!
//! ## Canonical form
//!
//! `(a, b, d)` is canonicalised by:
//! 1. Pulling out perfect-square integer factors of `d`'s numerator
//!    and denominator into `b` (so `d = 8` becomes `d = 2` with
//!    `b ↦ 2·b`).
//! 2. Setting `b = 0` and `d = 0` when the value is purely rational.
//! 3. After canonicalisation, two `QSqrtD` values are equal iff their
//!    `(a, b, d)` triples match exactly.

#![allow(dead_code)] // foundation — rank-2 eigenvector extraction layered on top

use super::rational::Rat;
use num_bigint::{BigInt, Sign};
use num_traits::{One, Signed, Zero};
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A value `a + b·√d` with `a, b ∈ Q` and `d ∈ Q` (canonicalised
/// square-free numerator and denominator).
#[derive(Clone, Debug)]
pub(crate) struct QSqrtD {
    /// Rational part.
    pub(crate) a: Rat,
    /// Coefficient of `√d`.
    pub(crate) b: Rat,
    /// The radicand (square-free in numerator and denominator).
    /// `b = 0` ⇒ `d` is irrelevant (we hold it at `0`).
    pub(crate) d: Rat,
}

impl QSqrtD {
    /// Canonical zero (`d = 0`, `b = 0`).
    pub(crate) fn zero() -> Self {
        Self {
            a: Rat::zero(),
            b: Rat::zero(),
            d: Rat::zero(),
        }
    }

    /// Pure rational `a` (no surd part).
    pub(crate) fn from_rat(a: Rat) -> Self {
        Self {
            a,
            b: Rat::zero(),
            d: Rat::zero(),
        }
    }

    /// `a + b·√d` after canonicalisation.
    pub(crate) fn new(a: Rat, b: Rat, d: Rat) -> Self {
        let mut v = Self { a, b, d };
        v.canonicalise();
        v
    }

    /// True iff this value is exactly the rational zero.
    pub(crate) fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    /// True iff this value is purely rational (no surd part).
    pub(crate) fn is_rational(&self) -> bool {
        self.b.is_zero()
    }

    /// Conjugate `(a + b√d)* = a − b√d`. Used by division and norm.
    pub(crate) fn conjugate(&self) -> Self {
        Self {
            a: self.a.clone(),
            b: -self.b.clone(),
            d: self.d.clone(),
        }
    }

    /// Norm `(a + b√d)·(a − b√d) = a² − b²·d`. Always rational.
    pub(crate) fn norm(&self) -> Rat {
        &self.a * &self.a - &self.b * &self.b * &self.d
    }

    /// Pull out the largest perfect-square integer factor of `d`'s
    /// numerator and denominator into `b`. After this, `d` is
    /// "square-free" up to the rational structure; two `QSqrtD`
    /// values that should agree have matching `d`.
    fn canonicalise(&mut self) {
        if self.b.is_zero() {
            self.b = Rat::zero();
            self.d = Rat::zero();
            return;
        }
        if self.d.is_zero() {
            // `b·√0 = 0`, so the surd part vanishes.
            self.b = Rat::zero();
            self.d = Rat::zero();
            return;
        }
        // Sign of d: factor out −1 if d < 0 by flipping the radicand.
        // (This module only supports real surds — d ≥ 0.)
        if self.d.is_negative() {
            // Real-only: leave the negative `d` in place; downstream
            // checks should reject. Keeping it explicit avoids silent
            // complex-number leakage.
            return;
        }
        // Square-free reduction on numerator and denominator separately.
        let num = self.d.numer().clone();
        let den = self.d.denom().clone();
        let (num_squarefree, num_extracted) = extract_square_factor(&num);
        let (den_squarefree, den_extracted) = extract_square_factor(&den);
        let new_d = Rat::new(num_squarefree, den_squarefree);
        let extracted = Rat::new(num_extracted, den_extracted);
        self.b = &self.b * &extracted;
        self.d = new_d;
    }
}

/// For a positive integer `n`, return `(s, e)` such that `n = e²·s`
/// with `s` square-free. (Equivalently, `s` has no integer factor
/// that is a perfect square > 1.)
fn extract_square_factor(n: &BigInt) -> (BigInt, BigInt) {
    let abs_n = n.abs();
    if abs_n <= BigInt::one() {
        return (n.clone(), BigInt::one());
    }
    let mut s = abs_n.clone();
    let mut e = BigInt::one();
    let mut p = BigInt::from(2);
    while &p * &p <= s {
        let p_sq = &p * &p;
        while (&s % &p_sq).is_zero() {
            s /= &p_sq;
            e *= &p;
        }
        p += BigInt::one();
    }
    // Carry the original sign on the square-free part.
    let s_signed = match n.sign() {
        Sign::Minus => -s,
        _ => s,
    };
    (s_signed, e)
}

impl PartialEq for QSqrtD {
    fn eq(&self, other: &Self) -> bool {
        // Two canonical values are equal iff (a, b, d) triples match.
        // If one or both have b = 0, only `a` needs to match.
        if self.b.is_zero() && other.b.is_zero() {
            return self.a == other.a;
        }
        if self.b.is_zero() != other.b.is_zero() {
            return false;
        }
        self.a == other.a && self.b == other.b && self.d == other.d
    }
}

impl Eq for QSqrtD {}

impl Neg for QSqrtD {
    type Output = QSqrtD;
    fn neg(self) -> Self::Output {
        QSqrtD {
            a: -self.a,
            b: -self.b,
            d: self.d,
        }
    }
}

impl Add for QSqrtD {
    type Output = QSqrtD;
    fn add(self, rhs: Self) -> Self::Output {
        // Pure rationals always combine.
        if self.b.is_zero() {
            return QSqrtD::new(self.a + rhs.a, rhs.b, rhs.d);
        }
        if rhs.b.is_zero() {
            return QSqrtD::new(self.a + rhs.a, self.b, self.d);
        }
        // Both have surds — `d` must agree.
        assert_eq!(
            self.d, rhs.d,
            "QSqrtD::add: incompatible radicands d_lhs={}, d_rhs={}",
            self.d, rhs.d
        );
        QSqrtD::new(self.a + rhs.a, self.b + rhs.b, self.d)
    }
}

impl Sub for QSqrtD {
    type Output = QSqrtD;
    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for QSqrtD {
    type Output = QSqrtD;
    fn mul(self, rhs: Self) -> Self::Output {
        // Cross product: (a + b√d)·(c + e√d) = ac + bed + (ae + bc)√d.
        if self.b.is_zero() {
            // self is pure rational `a`. (a)·(c + e√d) = ac + ae·√d.
            return QSqrtD::new(&self.a * &rhs.a, &self.a * &rhs.b, rhs.d);
        }
        if rhs.b.is_zero() {
            return QSqrtD::new(&self.a * &rhs.a, &self.b * &rhs.a, self.d);
        }
        assert_eq!(
            self.d, rhs.d,
            "QSqrtD::mul: incompatible radicands d_lhs={}, d_rhs={}",
            self.d, rhs.d
        );
        let rat = &self.a * &rhs.a + &self.b * &rhs.b * &self.d;
        let surd = &self.a * &rhs.b + &self.b * &rhs.a;
        QSqrtD::new(rat, surd, self.d)
    }
}

impl Div for QSqrtD {
    type Output = QSqrtD;
    fn div(self, rhs: Self) -> Self::Output {
        // Multiply by conjugate to rationalise the denominator.
        // (a + b√d) / (c + e√d) = (a + b√d)·(c − e√d) / (c² − e²·d).
        let n = rhs.norm();
        assert!(!n.is_zero(), "QSqrtD::div by zero");
        let conj = rhs.conjugate();
        let num = self * conj;
        let inv_n = Rat::one() / n;
        QSqrtD::new(num.a * inv_n.clone(), num.b * inv_n, num.d)
    }
}

/// Convert this `Q[√d]` value to `f64` via the natural embedding.
/// Truncation error of double-precision arithmetic — useful only for
/// final output to numerical layers.
pub(crate) fn to_f64(v: &QSqrtD) -> f64 {
    use num_traits::ToPrimitive;
    let a = v.a.to_f64().unwrap_or(f64::NAN);
    if v.b.is_zero() {
        return a;
    }
    let b = v.b.to_f64().unwrap_or(f64::NAN);
    let d = v.d.to_f64().unwrap_or(f64::NAN);
    a + b * d.sqrt()
}
