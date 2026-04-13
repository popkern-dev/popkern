//! Adaptive precision for geometric predicates.
//!
//! When a geometric predicate evaluates to a value within floating-point
//! uncertainty of zero, the result is unreliable. This module provides
//! adaptive evaluation strategies:
//!
//! 1. **Fast path**: evaluate with standard f64. If the result is clearly
//!    nonzero (|result| > ε_robust), accept it.
//! 2. **Fallback**: if near zero, use an alternative strategy (perturbed
//!    evaluation, winding number, or multiple rays).
//!
//! This provides the correctness of exact arithmetic with the performance
//! of floating-point for the common (non-degenerate) case.
//!
//! ## Predicate primitives
//!
//! Low-level:
//! - [`adaptive_sign`] — sign of a scalar with an ambiguity zone.
//! - [`DoublePrecision`] — double-double arithmetic (≈106 bits of
//!   significand). Constructors: [`DoublePrecision::from_f64`],
//!   [`DoublePrecision::from_sub`], [`DoublePrecision::from_mul`].
//!   Operations: [`dd_add`](DoublePrecision::dd_add),
//!   [`dd_sub`](DoublePrecision::dd_sub),
//!   [`dd_mul`](DoublePrecision::dd_mul),
//!   [`dd_neg`](DoublePrecision::dd_neg).
//!
//! High-level (geometric):
//! - [`adaptive_dist_sq_lt`] — `|p1 - p2|² < threshold²` with an
//!   automatic escalation to double-precision when the f64 answer
//!   lands inside the ambiguity margin.
//! - [`adaptive_cross_ratio_sq_lt`] — `|n1 × n2|² / (|n1|²·|n2|²) <
//!   threshold²` (the sine-of-angle tangency test) with the same
//!   escalation pattern. Used for near-parallel-normal detection
//!   where catastrophic cancellation in the cross-product
//!   components would otherwise flip the sign.
//! - [`adaptive_abs_plane_distance_lt`] — `|(point - origin) ·
//!   normal| < threshold` (signed distance from a point to a
//!   plane) with escalation. Cancellation in the per-axis
//!   subtraction is captured via [`DoublePrecision::from_sub`],
//!   and the dot product is re-accumulated in DP so the sign
//!   flips that f64 can produce at the threshold boundary become
//!   impossible.

/// Evaluate a sign predicate with adaptive precision.
///
/// Computes the sign of `value`. If `|value| < epsilon`,
/// the result is ambiguous and `None` is returned.
///
/// This is the building block for higher-level geometric predicates
/// like orientation tests and containment checks.
pub fn adaptive_sign(value: f64, epsilon: f64) -> Option<std::cmp::Ordering> {
    if value > epsilon {
        Some(std::cmp::Ordering::Greater)
    } else if value < -epsilon {
        Some(std::cmp::Ordering::Less)
    } else {
        None // ambiguous
    }
}

/// Double-double arithmetic: represent a high-precision value as two f64s.
///
/// `value = hi + lo` where `|lo| <= 0.5 * ulp(hi)`.
/// This gives roughly 106 bits of significand (vs 53 for f64).
#[derive(Debug, Clone, Copy)]
pub struct DoublePrecision {
    /// High-order component.
    pub hi: f64,
    /// Low-order component (correction term).
    pub lo: f64,
}

impl DoublePrecision {
    /// Create from a single f64.
    pub fn from_f64(x: f64) -> Self {
        Self { hi: x, lo: 0.0 }
    }

    /// Compute `a - b` exactly as a double-precision value.
    ///
    /// Uses the error-free subtraction `(a - b) + err` pattern so
    /// catastrophic cancellation between nearby inputs is captured
    /// in the `lo` term. For typical inputs the `lo` term is zero
    /// or a few ulps of `hi`, but for subtractions like `1.0 - (1.0
    /// - 1e-16)` the difference between what f64 reports and what
    /// the mathematical result actually is lands in `lo`.
    pub fn from_sub(a: f64, b: f64) -> Self {
        let diff = a - b;
        // Error term: the bit-level discrepancy between the f64
        // result and the exact mathematical difference. When
        // |a| ≈ |b|, `diff` loses significant bits; `err` recovers
        // them.
        let err = (a - diff) - b;
        Self { hi: diff, lo: err }
    }

    /// Compute `a * b` exactly as a double-precision value via
    /// Dekker's TwoProduct with FMA.
    pub fn from_mul(a: f64, b: f64) -> Self {
        let p = a * b;
        let err = a.mul_add(b, -p);
        Self { hi: p, lo: err }
    }

    /// Negate.
    pub fn dd_neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }

    /// Add two double-precision values.
    pub fn dd_add(self, other: Self) -> Self {
        // Knuth's two-sum algorithm
        let s = self.hi + other.hi;
        let v = s - self.hi;
        let t = (self.hi - (s - v)) + (other.hi - v);
        let lo = t + self.lo + other.lo;
        // Normalize
        let hi = s + lo;
        let lo2 = lo - (hi - s);
        Self { hi, lo: lo2 }
    }

    /// Subtract `other` from `self`.
    pub fn dd_sub(self, other: Self) -> Self {
        self.dd_add(other.dd_neg())
    }

    /// Multiply two double-precision values.
    pub fn dd_mul(self, other: Self) -> Self {
        // Dekker's multiplication
        let p = self.hi * other.hi;
        let err = self.hi.mul_add(other.hi, -p); // FMA-based error
        let lo = err + self.hi * other.lo + self.lo * other.hi;
        let hi = p + lo;
        let lo2 = lo - (hi - p);
        Self { hi, lo: lo2 }
    }

    /// Get the approximate f64 value.
    pub fn to_f64(self) -> f64 {
        self.hi + self.lo
    }

    /// Get the sign, or None if ambiguous (|value| < epsilon).
    pub fn sign(self, epsilon: f64) -> Option<std::cmp::Ordering> {
        adaptive_sign(self.hi + self.lo, epsilon)
    }
}

/// Adaptive `|p1 - p2|² < threshold_sq` predicate.
///
/// Two-stage evaluation:
///
/// 1. **Fast path**: compute the squared distance and the threshold's
///    ambiguity margin in f64. If the result is outside the margin,
///    return the f64 answer directly.
///
/// 2. **Slow path**: if the f64 answer lies inside the margin,
///    recompute the squared distance using [`DoublePrecision`] for
///    the three subtractions and three products, then compare the
///    hi-precision result against the threshold. The DP recomputation
///    eliminates catastrophic cancellation in the `p1.i - p2.i`
///    subtractions, which is the dominant error source when `p1`
///    and `p2` are numerically close.
///
/// The f64 margin is `4 * f64::EPSILON * max(dist_sq, threshold_sq)`,
/// a standard conservative bound for the worst-case error of three
/// multiplies + two adds + one subtract in f64.
#[inline]
pub fn adaptive_dist_sq_lt(p1: (f64, f64, f64), p2: (f64, f64, f64), threshold_sq: f64) -> bool {
    let dx = p1.0 - p2.0;
    let dy = p1.1 - p2.1;
    let dz = p1.2 - p2.2;
    let dist_sq = dx * dx + dy * dy + dz * dz;
    let margin = 4.0 * f64::EPSILON * dist_sq.max(threshold_sq);

    if dist_sq + margin < threshold_sq {
        return true;
    }
    if dist_sq > threshold_sq + margin {
        return false;
    }

    // Slow path goes through a non-inline helper so the
    // inlined fast path stays compact.
    adaptive_dist_sq_lt_slow(p1, p2, threshold_sq)
}

#[inline(never)]
fn adaptive_dist_sq_lt_slow(p1: (f64, f64, f64), p2: (f64, f64, f64), threshold_sq: f64) -> bool {
    let dp_dx = DoublePrecision::from_sub(p1.0, p2.0);
    let dp_dy = DoublePrecision::from_sub(p1.1, p2.1);
    let dp_dz = DoublePrecision::from_sub(p1.2, p2.2);
    let dp_dist_sq = dp_dx
        .dd_mul(dp_dx)
        .dd_add(dp_dy.dd_mul(dp_dy))
        .dd_add(dp_dz.dd_mul(dp_dz));
    dp_dist_sq
        .dd_sub(DoublePrecision::from_f64(threshold_sq))
        .to_f64()
        < 0.0
}

/// Adaptive `|n1 × n2|² / (|n1|²·|n2|²) < threshold²` predicate.
///
/// This is the sine-of-angle tangency test: it returns `true`
/// when the two vectors are nearly parallel (within angular
/// tolerance `threshold`, measured as `sin(θ) < threshold`).
///
/// Rewrites the comparison to avoid the division:
///
/// ```text
/// |n1 × n2| / (|n1| * |n2|)      < threshold
/// ⇔ |n1 × n2|²                    < threshold² · |n1|² · |n2|²
/// ```
///
/// Two-stage evaluation:
///
/// 1. **Fast path**: compute the cross product components,
///    `|cross|²`, and `threshold² · |n1|² · |n2|²` in f64. If the
///    comparison is outside the margin (`12 * ε * max`), return
///    the f64 answer.
///
/// 2. **Slow path**: recompute the cross-product components using
///    [`DoublePrecision::from_mul`] and [`DoublePrecision::dd_sub`],
///    then sum their squares. Compare against
///    `threshold² · |n1|² · |n2|²` likewise computed in DP. Only
///    pays the slow-path cost when the f64 answer lands inside
///    the ambiguity margin.
///
/// Returns `false` if either input has effectively zero magnitude
/// (`|n1|² · |n2|²` below `f64::EPSILON`) — mirrors the existing
/// intersect.rs guard against dividing by a tiny number.
#[inline]
pub fn adaptive_cross_ratio_sq_lt(
    n1: (f64, f64, f64),
    n2: (f64, f64, f64),
    threshold: f64,
) -> bool {
    // Fast path: compute cross product and magnitudes in f64.
    let cx = n1.1 * n2.2 - n1.2 * n2.1;
    let cy = n1.2 * n2.0 - n1.0 * n2.2;
    let cz = n1.0 * n2.1 - n1.1 * n2.0;
    let cross_len_sq = cx * cx + cy * cy + cz * cz;

    let n1_len_sq = n1.0 * n1.0 + n1.1 * n1.1 + n1.2 * n1.2;
    let n2_len_sq = n2.0 * n2.0 + n2.1 * n2.1 + n2.2 * n2.2;
    let norms_sq = n1_len_sq * n2_len_sq;

    // Degenerate: one of the vectors has zero magnitude. Return
    // false so callers don't classify as tangent.
    if norms_sq <= f64::EPSILON {
        return false;
    }

    let threshold_sq = threshold * threshold;
    let rhs = threshold_sq * norms_sq;

    // Margin: 3 multiplies + 2 adds for cross² (≈8 ulps) + 1
    // multiply on the threshold (≈1 ulp). Conservative bound:
    // 12 * ε * max(cross_len_sq, rhs).
    let margin = 12.0 * f64::EPSILON * cross_len_sq.max(rhs);

    if cross_len_sq + margin < rhs {
        return true;
    }
    if cross_len_sq > rhs + margin {
        return false;
    }

    // Slow path: recompute cross product in double precision.
    let dp_n1_0 = DoublePrecision::from_f64(n1.0);
    let dp_n1_1 = DoublePrecision::from_f64(n1.1);
    let dp_n1_2 = DoublePrecision::from_f64(n1.2);
    let dp_n2_0 = DoublePrecision::from_f64(n2.0);
    let dp_n2_1 = DoublePrecision::from_f64(n2.1);
    let dp_n2_2 = DoublePrecision::from_f64(n2.2);

    // cx = n1.1 * n2.2 - n1.2 * n2.1 — use from_mul to get exact
    // products, then dd_sub.
    let dp_cx = DoublePrecision::from_mul(n1.1, n2.2).dd_sub(DoublePrecision::from_mul(n1.2, n2.1));
    let dp_cy = DoublePrecision::from_mul(n1.2, n2.0).dd_sub(DoublePrecision::from_mul(n1.0, n2.2));
    let dp_cz = DoublePrecision::from_mul(n1.0, n2.1).dd_sub(DoublePrecision::from_mul(n1.1, n2.0));

    let dp_cross_len_sq = dp_cx
        .dd_mul(dp_cx)
        .dd_add(dp_cy.dd_mul(dp_cy))
        .dd_add(dp_cz.dd_mul(dp_cz));

    let dp_n1_len_sq = dp_n1_0
        .dd_mul(dp_n1_0)
        .dd_add(dp_n1_1.dd_mul(dp_n1_1))
        .dd_add(dp_n1_2.dd_mul(dp_n1_2));
    let dp_n2_len_sq = dp_n2_0
        .dd_mul(dp_n2_0)
        .dd_add(dp_n2_1.dd_mul(dp_n2_1))
        .dd_add(dp_n2_2.dd_mul(dp_n2_2));
    let dp_norms_sq = dp_n1_len_sq.dd_mul(dp_n2_len_sq);
    let dp_rhs = dp_norms_sq.dd_mul(DoublePrecision::from_f64(threshold_sq));

    dp_cross_len_sq.dd_sub(dp_rhs).to_f64() < 0.0
}

/// Adaptive `|(point - origin) · normal| < threshold` predicate.
///
/// Computes the signed distance from `point` to the plane defined
/// by `origin` and `normal` (assumed to be a unit vector), then
/// compares the magnitude against `threshold`.
///
/// Two-stage evaluation:
/// 1. **Fast path**: compute the f64 signed distance and check
///    against `threshold` with a conservative 6·ε margin.
/// 2. **Slow path**: recompute each `point.i - origin.i`
///    subtraction with [`DoublePrecision::from_sub`] so any
///    cancellation bits are recovered, then re-accumulate the
///    dot product in DP and compare against `threshold`.
#[inline]
pub fn adaptive_abs_plane_distance_lt(
    point: (f64, f64, f64),
    origin: (f64, f64, f64),
    normal: (f64, f64, f64),
    threshold: f64,
) -> bool {
    // Fast path.
    let dx = point.0 - origin.0;
    let dy = point.1 - origin.1;
    let dz = point.2 - origin.2;
    let f = dx * normal.0 + dy * normal.1 + dz * normal.2;
    let f_abs = f.abs();
    let margin = 6.0 * f64::EPSILON * f_abs.max(threshold);

    if f_abs + margin < threshold {
        return true;
    }
    if f_abs > threshold + margin {
        return false;
    }

    // Slow path: recompute in double precision.
    let dp_dx = DoublePrecision::from_sub(point.0, origin.0);
    let dp_dy = DoublePrecision::from_sub(point.1, origin.1);
    let dp_dz = DoublePrecision::from_sub(point.2, origin.2);
    let dp_f = dp_dx
        .dd_mul(DoublePrecision::from_f64(normal.0))
        .dd_add(dp_dy.dd_mul(DoublePrecision::from_f64(normal.1)))
        .dd_add(dp_dz.dd_mul(DoublePrecision::from_f64(normal.2)));
    dp_f.to_f64().abs() < threshold
}

#[cfg(test)]
#[path = "adaptive_precision_tests.rs"]
mod adaptive_precision_tests;
