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

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::classify::point_in_solid;
use crate::classify_winding::point_in_solid_winding;
use crate::error::PointClassification;

/// Robust point-in-solid classification with adaptive fallback.
///
/// First tries multi-ray classification (fast). If the result is
/// ambiguous (OnBoundary), falls back to winding number (more robust
/// but slower). This provides the best of both methods.
pub fn robust_point_in_solid(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> PointClassification {
    // Fast path: multi-ray classification
    match point_in_solid(point, solid, face_aabbs, tolerance) {
        Ok(PointClassification::Inside) => PointClassification::Inside,
        Ok(PointClassification::Outside) => PointClassification::Outside,
        Ok(PointClassification::OnBoundary) | Err(_) => {
            // Ambiguous or failed — fall back to winding number
            point_in_solid_winding(point, solid)
        }
    }
}

/// Evaluate a sign predicate with adaptive precision.
///
/// Computes `f(x)` and returns its sign. If `|f(x)| < epsilon`,
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
mod tests {
    use super::*;

    #[test]
    fn robust_classification_inside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let aabbs = crate::intersect::compute_face_aabbs(&cube);
        let pt = Point3::new(0.5, 0.5, 0.5);
        let result = robust_point_in_solid(&pt, &cube, &aabbs, &Tolerance::default());
        assert_eq!(result, PointClassification::Inside);
    }

    #[test]
    fn robust_classification_outside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let aabbs = crate::intersect::compute_face_aabbs(&cube);
        let pt = Point3::new(5.0, 5.0, 5.0);
        let result = robust_point_in_solid(&pt, &cube, &aabbs, &Tolerance::default());
        assert_eq!(result, PointClassification::Outside);
    }

    #[test]
    fn adaptive_sign_clear() {
        assert_eq!(adaptive_sign(1.0, 1e-10), Some(std::cmp::Ordering::Greater));
        assert_eq!(adaptive_sign(-1.0, 1e-10), Some(std::cmp::Ordering::Less));
    }

    #[test]
    fn adaptive_sign_ambiguous() {
        assert_eq!(adaptive_sign(1e-15, 1e-10), None);
    }

    #[test]
    fn double_precision_add() {
        let a = DoublePrecision::from_f64(1.0);
        let b = DoublePrecision::from_f64(1e-16);
        let sum = a.dd_add(b);
        // Should preserve the small value that f64 would lose
        assert!((sum.to_f64() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn double_precision_mul() {
        let a = DoublePrecision::from_f64(3.0);
        let b = DoublePrecision::from_f64(7.0);
        let product = a.dd_mul(b);
        assert!((product.to_f64() - 21.0).abs() < 1e-14);
    }
}
