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
//! ## Plane extraction (stage 2b.45)
//!
//! - `try_sqrt_rat_in_qsqrtd(r, d) -> Option<QSqrtD>` — try to
//!   express `√r` as a value in `Q[√d]` (success cases: `r` a
//!   rational square, or `r/d` a rational square ⇒ `√r = c·√d`).
//! - `extract_planes_qsqrtd_nested(D) -> Option<[[QSqrtDSqrtL; 4]; 2]>`
//!   — companion to `qsqrtd_decompose::extract_planes_qsqrtd` that
//!   uses the nested field. Closes the golden-ratio class where the
//!   single-surd path returns `None` because the eigenvalues aren't
//!   perfect squares in `Q[√d]`. Limitation: requires `√(−S₂) ∈ Q[√d]`
//!   so the second eigenvalue's sqrt embeds in the same nested field
//!   as `√μ_+`; otherwise a triple tower would be needed (out of
//!   scope).

#![allow(dead_code)]

use std::ops::{Add, Div, Mul, Neg, Sub};

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use super::qsqrtd::QSqrtD;
use super::qsqrtd_decompose::{
    Mat4SymQSqrtD, lift_to_qsqrtd, null_vector_qsqrtd, rank2_eigenvalues, subtract_mu_diag,
};
use super::rational::{Mat4SymRat, Rat};

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

// ─── Stage 2b.45: nested-radical plane extraction ─────────────────

/// Exact integer square root: returns `Some(s)` iff `n = s²` for some
/// non-negative integer `s`; `None` otherwise. Newton's method.
fn isqrt_exact(n: &BigInt) -> Option<BigInt> {
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    if n.is_negative() {
        return None;
    }
    let mut x = n.clone();
    let mut y = (&x + BigInt::one()) >> 1;
    while y < x {
        x = y;
        y = (&x + n / &x) >> 1;
    }
    if &x * &x == *n { Some(x) } else { None }
}

/// Exact rational square root: returns `Some(s)` iff `r = s²` for
/// some non-negative rational `s`; `None` otherwise. Requires both
/// numerator and denominator of `r` (in lowest terms) to be perfect
/// squares.
fn isqrt_rat(r: &Rat) -> Option<Rat> {
    if r.is_negative() {
        return None;
    }
    if r.is_zero() {
        return Some(Rat::zero());
    }
    let n_sqrt = isqrt_exact(r.numer())?;
    let d_sqrt = isqrt_exact(r.denom())?;
    Some(Rat::new(n_sqrt, d_sqrt))
}

/// Try to express `√r` as a value in `Q[√d]` for rational `r ≥ 0`
/// and rational `d > 0`. Two cases succeed:
/// 1. `r` is a perfect rational square ⇒ `√r ∈ Q ⊂ Q[√d]`.
/// 2. `r/d` is a perfect rational square ⇒ `√r = c·√d` for `c ∈ Q`.
///
/// Returns `None` when `√r` requires a fresh radical extension
/// distinct from `√d` (e.g., `r = 2`, `d = 3` — `√2` and `√3` are
/// independent surds, both irrational and unrelated).
pub(crate) fn try_sqrt_rat_in_qsqrtd(r: &Rat, d: &Rat) -> Option<QSqrtD> {
    if r.is_negative() {
        return None;
    }
    // Case 1: √r ∈ Q.
    if let Some(s) = isqrt_rat(r) {
        return Some(QSqrtD::from_rat(s));
    }
    // Case 2: √r = c·√d.
    if d.is_zero() {
        return None;
    }
    let r_div_d = r / d;
    if let Some(c) = isqrt_rat(&r_div_d) {
        return Some(QSqrtD::new(Rat::zero(), c, d.clone()));
    }
    None
}

/// Symbolically extract the two **plane normals** from a rank-2
/// mixed-sign quadric matrix `D`, using the nested-radical extension
/// `Q[√d, √μ_+]` to handle eigenvalues that are *not* perfect squares
/// in `Q[√d]`. Companion to `qsqrtd_decompose::extract_planes_qsqrtd`,
/// which only succeeds when `√μ_+` and `√(−μ_-)` already live in
/// `Q[√d]`.
///
/// ## When this succeeds
///
/// The returned planes live in `Q[√d, √μ_+]` (with `α = μ_+` the
/// second radicand). Both square roots embed in this single field
/// iff `√(−S₂) ∈ Q[√d]`, where `S₂ = μ_+ · μ_- ∈ Q` is the product
/// of the non-zero eigenvalues. Then:
/// - `√μ_+ = 0 + 1·√α`,
/// - `√(−μ_-) = √(−S₂) / √μ_+ = (√(−S₂) / α)·√α` (one component, in `Q[√d]`).
///
/// ## When this returns `None`
///
/// - `disc < 0`: complex eigenvalues (no real plane decomposition).
/// - `disc = 0`: rank-1 (doubled plane — use `tangent::extract_plane_from_rank1`).
/// - `μ_+ = 0`: degenerate.
/// - Null vectors not extractable (matrix not rank-2 mixed-sign).
/// - `√(−S₂)` is genuinely irrational over `Q[√d]` (would need a
///   triple tower `Q[√d, √μ_+, √(−μ_-)]` — out of scope).
///
/// The "golden-ratio" class — eigenvalues `(1 ± √5)/2` from an
/// integer-trace, integer-determinant 2×2 reduction — has
/// `S₂ = −1` and so `−S₂ = 1` is trivially a square in `Q ⊂ Q[√d]`,
/// so this function succeeds where the single-surd path fails.
pub(crate) fn extract_planes_qsqrtd_nested(d: &Mat4SymRat) -> Option<[[QSqrtDSqrtL; 4]; 2]> {
    let (mu_plus, mu_minus, disc) = rank2_eigenvalues(d)?;
    if disc.is_zero() {
        return None;
    }
    if mu_plus.is_zero() {
        return None;
    }

    let alpha = mu_plus.clone();
    let lifted: Mat4SymQSqrtD = lift_to_qsqrtd(d);
    let m_plus = subtract_mu_diag(&lifted, &mu_plus);
    let v_plus = null_vector_qsqrtd(&m_plus)?;
    let m_minus = subtract_mu_diag(&lifted, &mu_minus);
    let v_minus = null_vector_qsqrtd(&m_minus)?;

    // S₂ = μ_+·μ_-; the cross-terms in (a + b√d)(a' + b'√d) with
    // mu_minus = (tr/2) − (1/2)√disc cancel, so the product is pure-
    // rational by Vieta. Verify before extracting.
    let s2_qsqrtd = mu_plus.clone() * mu_minus.clone();
    if !s2_qsqrtd.b.is_zero() {
        // Shouldn't happen for the well-formed rank-2 case, but guard
        // anyway against numeric surprises in degenerate inputs.
        return None;
    }
    let neg_s2: Rat = -s2_qsqrtd.a.clone();

    // sqrt_neg_s2 must live in Q[√d] — i.e., either √(−S₂) ∈ Q
    // (case 1) or √(−S₂) = c·√d (case 2). Reject otherwise.
    let sqrt_neg_s2 = try_sqrt_rat_in_qsqrtd(&neg_s2, &alpha.d)?;

    // √μ_+ = 0 + 1·√α (canonical embedding of the second radicand).
    let sqrt_lp = QSqrtDSqrtL::new(QSqrtD::zero(), QSqrtD::from_rat(Rat::one()), alpha.clone());
    // √(−μ_-) = √(−S₂)/√α = (√(−S₂)/α)·√α.
    let sqrt_lm = QSqrtDSqrtL::new(QSqrtD::zero(), sqrt_neg_s2 / alpha.clone(), alpha.clone());

    let scale_v = |s: &QSqrtDSqrtL, v: &[QSqrtD; 4]| -> [QSqrtDSqrtL; 4] {
        [
            s.clone() * QSqrtDSqrtL::from_qsqrtd(v[0].clone(), alpha.clone()),
            s.clone() * QSqrtDSqrtL::from_qsqrtd(v[1].clone(), alpha.clone()),
            s.clone() * QSqrtDSqrtL::from_qsqrtd(v[2].clone(), alpha.clone()),
            s.clone() * QSqrtDSqrtL::from_qsqrtd(v[3].clone(), alpha.clone()),
        ]
    };
    let a = scale_v(&sqrt_lp, &v_plus);
    let b = scale_v(&sqrt_lm, &v_minus);
    let p1: [QSqrtDSqrtL; 4] = [
        a[0].clone() + b[0].clone(),
        a[1].clone() + b[1].clone(),
        a[2].clone() + b[2].clone(),
        a[3].clone() + b[3].clone(),
    ];
    let p2: [QSqrtDSqrtL; 4] = [
        a[0].clone() - b[0].clone(),
        a[1].clone() - b[1].clone(),
        a[2].clone() - b[2].clone(),
        a[3].clone() - b[3].clone(),
    ];
    Some([p1, p2])
}
