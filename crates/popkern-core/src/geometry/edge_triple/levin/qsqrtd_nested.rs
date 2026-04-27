//! Stage 2b.43: nested-radical extension `Q[√d, √α]`.
//!
//! The single-surd field `Q[√d]` (in `qsqrtd.rs`) handles the
//! discriminant of the rank-2 reduction's 2×2 sub-block. To
//! parameterise the plane normals via `a = √λ_+ · v_+`, however, we
//! sometimes need `√(λ_+)` where `λ_+ ∈ Q[√d]` is itself **not** a
//! perfect square in `Q[√d]` (e.g., the golden-ratio class
//! `λ_+ = (1 + √5)/2`). The square root then lives in a **second**
//! algebraic extension on top of the first.
//!
//! This module gives the field type for that nested extension. A
//! value is `p₀ + p₁·√α` with `p₀, p₁, α ∈ Q[√d]` for a fixed shared
//! base radicand `d`. The second radicand `α` is per-instance —
//! operations between values with mismatched `α` panic (caller should
//! lift into a common `α` first, or use `qsqrtd::QSqrtD` if both
//! operands are flat in `Q[√d]`).
//!
//! ## Field operations
//!
//! - `+`, `-`, `-(unary)` work componentwise.
//! - Multiplication uses `(p₀ + p₁√α)(q₀ + q₁√α) = (p₀·q₀ + α·p₁·q₁)
//!   + (p₀·q₁ + p₁·q₀)·√α`.
//! - Division multiplies through by the conjugate `(q₀ − q₁√α)` and
//!   divides each component by the resulting `Q[√d]` norm
//!   `q₀² − α·q₁²`.
//!
//! ## Equality
//!
//! Equality is **structural** on `(p₀, p₁, α)` — i.e., this module
//! assumes `√α` is genuinely irrational over `Q[√d]`. If the caller
//! passes a perfect-square `α`, the value collapses to `Q[√d]` and
//! should be represented as a `QSqrtD` instead; structural equality
//! would otherwise miss algebraically-equal pairs.
//!
//! Out of scope for this stage: integration with `qsqrtd_decompose`'s
//! plane extraction (will land in 2b.44 once the field foundation is
//! verified by the unit tests in `qsqrtd_nested_tests.rs`).

#![allow(dead_code)]

use std::ops::{Add, Div, Mul, Neg, Sub};

use super::qsqrtd::QSqrtD;
use super::rational::Rat;

/// A value `p₀ + p₁·√α` with `p₀, p₁, α ∈ Q[√d]`.
#[derive(Clone, Debug)]
pub(crate) struct QSqrtDSqrtL {
    /// Constant part (in `Q[√d]`).
    pub(crate) p0: QSqrtD,
    /// Coefficient of `√α` (in `Q[√d]`).
    pub(crate) p1: QSqrtD,
    /// The second radicand, `α ∈ Q[√d]`. Two `QSqrtDSqrtL` values
    /// must share the same `α` for arithmetic to apply.
    pub(crate) alpha: QSqrtD,
}

impl QSqrtDSqrtL {
    /// Construct `p₀ + p₁·√α`.
    pub(crate) fn new(p0: QSqrtD, p1: QSqrtD, alpha: QSqrtD) -> Self {
        Self { p0, p1, alpha }
    }

    /// Construct the additive zero in `Q[√d, √α]`.
    pub(crate) fn zero(alpha: QSqrtD) -> Self {
        Self {
            p0: QSqrtD::zero(),
            p1: QSqrtD::zero(),
            alpha,
        }
    }

    /// Lift a `Q[√d]` value into `Q[√d, √α]` (with `p₁ = 0`).
    pub(crate) fn from_qsqrtd(q: QSqrtD, alpha: QSqrtD) -> Self {
        Self {
            p0: q,
            p1: QSqrtD::zero(),
            alpha,
        }
    }

    /// Lift a rational into `Q[√d, √α]` (with `p₁ = 0` and pure-rat
    /// `p₀`). Convenient when the caller already has a base radicand
    /// `α` and wants integer / fractional coefficients.
    pub(crate) fn from_rat(r: Rat, alpha: QSqrtD) -> Self {
        Self {
            p0: QSqrtD::from_rat(r),
            p1: QSqrtD::zero(),
            alpha,
        }
    }

    /// True iff the value is exactly zero (both components vanish).
    pub(crate) fn is_zero(&self) -> bool {
        self.p0.is_zero() && self.p1.is_zero()
    }

    /// `(p₀ − p₁·√α)` — the conjugate over the second extension.
    pub(crate) fn conjugate(&self) -> Self {
        Self {
            p0: self.p0.clone(),
            p1: -self.p1.clone(),
            alpha: self.alpha.clone(),
        }
    }

    /// `p₀² − α·p₁²` — the field norm down to `Q[√d]`. Equals
    /// `self * self.conjugate()`'s `p0` component (with `p1 = 0`),
    /// and is the denominator that `Div` factors out via
    /// rationalisation.
    pub(crate) fn norm(&self) -> QSqrtD {
        self.p0.clone() * self.p0.clone()
            - self.alpha.clone() * self.p1.clone() * self.p1.clone()
    }
}

impl PartialEq for QSqrtDSqrtL {
    fn eq(&self, other: &Self) -> bool {
        // Structural — assumes `√α` is irrational over `Q[√d]` (see
        // module-level docs). Mismatched `α` is unequal even if both
        // values happen to reduce to `Q[√d]`.
        self.alpha == other.alpha && self.p0 == other.p0 && self.p1 == other.p1
    }
}

impl Eq for QSqrtDSqrtL {}

impl Neg for QSqrtDSqrtL {
    type Output = QSqrtDSqrtL;
    fn neg(self) -> Self::Output {
        Self {
            p0: -self.p0,
            p1: -self.p1,
            alpha: self.alpha,
        }
    }
}

impl Add for QSqrtDSqrtL {
    type Output = QSqrtDSqrtL;
    fn add(self, other: Self) -> Self::Output {
        assert_eq!(
            self.alpha, other.alpha,
            "QSqrtDSqrtL::add — mismatched second radicand α"
        );
        Self {
            p0: self.p0 + other.p0,
            p1: self.p1 + other.p1,
            alpha: self.alpha,
        }
    }
}

impl Sub for QSqrtDSqrtL {
    type Output = QSqrtDSqrtL;
    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(
            self.alpha, other.alpha,
            "QSqrtDSqrtL::sub — mismatched second radicand α"
        );
        Self {
            p0: self.p0 - other.p0,
            p1: self.p1 - other.p1,
            alpha: self.alpha,
        }
    }
}

impl Mul for QSqrtDSqrtL {
    type Output = QSqrtDSqrtL;
    fn mul(self, other: Self) -> Self::Output {
        assert_eq!(
            self.alpha, other.alpha,
            "QSqrtDSqrtL::mul — mismatched second radicand α"
        );
        // (p₀ + p₁·√α)(q₀ + q₁·√α)
        //   = p₀·q₀ + α·p₁·q₁ + (p₀·q₁ + p₁·q₀)·√α
        let p0q0 = self.p0.clone() * other.p0.clone();
        let p1q1 = self.p1.clone() * other.p1.clone();
        let alpha_p1q1 = self.alpha.clone() * p1q1;
        let new_p0 = p0q0 + alpha_p1q1;
        let new_p1 = self.p0 * other.p1 + self.p1 * other.p0;
        Self {
            p0: new_p0,
            p1: new_p1,
            alpha: self.alpha,
        }
    }
}

impl Div for QSqrtDSqrtL {
    type Output = QSqrtDSqrtL;
    fn div(self, other: Self) -> Self::Output {
        assert_eq!(
            self.alpha, other.alpha,
            "QSqrtDSqrtL::div — mismatched second radicand α"
        );
        // num · conj(denom) / norm(denom)
        let conj = other.conjugate();
        let denom_norm = other.norm();
        assert!(
            !denom_norm.is_zero(),
            "QSqrtDSqrtL::div — divisor has zero norm (degenerate or zero)"
        );
        let prod = self * conj;
        Self {
            p0: prod.p0 / denom_norm.clone(),
            p1: prod.p1 / denom_norm,
            alpha: prod.alpha,
        }
    }
}
