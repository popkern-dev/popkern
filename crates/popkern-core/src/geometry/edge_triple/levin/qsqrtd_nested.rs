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
//! ## Helpers (stage 2b.44)
//!
//! - `is_non_negative(λ)` — sign of `a + b·√d` with `√d` taken on the
//!   positive branch (the standard convention).
//! - `sqrt_qsqrtd(λ) -> Option<QSqrtDSqrtL>` — symbolic square root.
//!   Always succeeds for non-negative `λ`, embedding `√λ` as the
//!   canonical `0 + 1·√λ` in `Q[√d, √λ]`. `(√λ)² = λ` by the
//!   defining identity. Returns `None` for negative inputs.
//!
//! Plumbing into plane extraction (`extract_planes_qsqrtd_nested`)
//! lands in a follow-up stage so the helpers above can be exercised
//! standalone first.

#![allow(dead_code)]

use std::ops::{Add, Div, Mul, Neg, Sub};

use num_traits::{One, Signed, Zero};

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
        self.p0.clone() * self.p0.clone() - self.alpha.clone() * self.p1.clone() * self.p1.clone()
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

/// True iff `λ ≥ 0` as a real number, taking `√d` to be the positive
/// branch (the standard convention; `Q[√d]` carries no orientation
/// otherwise). Used by `sqrt_qsqrtd` to reject negative inputs before
/// embedding them as a real-valued nested radical.
///
/// Decision logic for `λ = a + b·√d` with `d ≥ 0`:
/// - `b = 0`: sign matches `a`'s sign.
/// - `b > 0`: non-negative iff `a ≥ 0`, or (when `a < 0`) `a² ≤ b²·d`.
/// - `b < 0`: non-negative iff `a ≥ 0` AND `a² ≥ b²·d`.
pub(crate) fn is_non_negative(lambda: &QSqrtD) -> bool {
    if lambda.b.is_zero() {
        return !lambda.a.is_negative();
    }
    let a_sq = &lambda.a * &lambda.a;
    let b_sq_d = &lambda.b * &lambda.b * &lambda.d;
    if lambda.b.is_positive() {
        !lambda.a.is_negative() || a_sq <= b_sq_d
    } else {
        !lambda.a.is_negative() && a_sq >= b_sq_d
    }
}

/// Compute `√λ` symbolically as a `QSqrtDSqrtL`.
///
/// Always succeeds for non-negative `λ`. Returns `None` only when
/// `λ` is negative (real square root undefined).
///
/// Stage 2b.43 keeps the implementation minimal: the result embeds
/// `λ` as the second radicand and writes `√λ = 0 + 1·√λ`. By the
/// nested field's defining identity `(√α)² = α`, squaring the result
/// recovers `λ` exactly. Denesting (returning a value in `Q[√d]`
/// when one exists) is a separate optimisation deferred to the stage
/// that wires this into plane extraction.
pub(crate) fn sqrt_qsqrtd(lambda: &QSqrtD) -> Option<QSqrtDSqrtL> {
    if !is_non_negative(lambda) {
        return None;
    }
    if lambda.is_zero() {
        // √0 = 0; the result has no algebraic dependency on a
        // particular `α`, so we use the trivial radicand.
        return Some(QSqrtDSqrtL::zero(QSqrtD::zero()));
    }
    Some(QSqrtDSqrtL::new(
        QSqrtD::zero(),
        QSqrtD::from_rat(Rat::one()),
        lambda.clone(),
    ))
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
