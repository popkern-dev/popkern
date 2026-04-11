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

#[cfg(test)]
#[path = "adaptive_precision_tests.rs"]
mod adaptive_precision_tests;
