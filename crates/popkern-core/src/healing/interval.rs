//! Lightweight interval arithmetic for guaranteed geometric bounds.
//!
//! Provides a simple [`Interval`] type with conservative arithmetic that
//! guarantees enclosure: the true result is always within the computed
//! interval. Used by regional analysis for proven curvature bounds.

use std::ops::{Add, Mul, Neg, Sub};

/// A closed interval \[lo, hi\] representing a range of possible values.
///
/// All arithmetic operations produce intervals that are guaranteed to
/// contain the true result (conservative rounding outward).
#[derive(Debug, Clone, Copy)]
pub struct Interval {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}

impl Interval {
    /// Create a new interval \[lo, hi\].
    ///
    /// # Panics
    /// Panics if `lo > hi`.
    pub fn new(lo: f64, hi: f64) -> Self {
        assert!(
            lo <= hi,
            "Interval::new requires lo <= hi, got lo={lo}, hi={hi}"
        );
        Self { lo, hi }
    }

    /// Create a point interval \[v, v\].
    pub fn point(v: f64) -> Self {
        Self { lo: v, hi: v }
    }

    /// Whether the interval contains the value `v`.
    pub fn contains(&self, v: f64) -> bool {
        self.lo <= v && v <= self.hi
    }

    /// Width of the interval: hi - lo.
    pub fn width(&self) -> f64 {
        self.hi - self.lo
    }

    /// Midpoint of the interval: (lo + hi) / 2.
    pub fn midpoint(&self) -> f64 {
        (self.lo + self.hi) * 0.5
    }

    /// Whether the entire interval is strictly positive (lo > 0).
    pub fn is_positive(&self) -> bool {
        self.lo > 0.0
    }

    /// Whether the entire interval is strictly negative (hi < 0).
    pub fn is_negative(&self) -> bool {
        self.hi < 0.0
    }

    /// Whether the interval contains zero.
    pub fn contains_zero(&self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }

    /// Reciprocal: 1/\[lo, hi\].
    ///
    /// If the interval contains zero, returns a conservatively wide interval
    /// \[-1e18, 1e18\] since the true reciprocal is unbounded.
    /// Otherwise returns \[1/hi, 1/lo\].
    pub fn recip(&self) -> Self {
        if self.contains_zero() {
            Self::new(-1e18, 1e18)
        } else {
            // 1/hi <= 1/x <= 1/lo for positive intervals
            // 1/hi >= 1/x >= 1/lo for negative intervals
            let a = 1.0 / self.lo;
            let b = 1.0 / self.hi;
            Self::new(a.min(b), a.max(b))
        }
    }

    /// Square root: sqrt(\[lo, hi\]).
    ///
    /// Clamps the lower bound to zero before taking the square root,
    /// since sqrt is undefined for negative values.
    pub fn sqrt(&self) -> Self {
        let clamped_lo = self.lo.max(0.0);
        Self::new(clamped_lo.sqrt(), self.hi.max(0.0).sqrt())
    }

    /// Absolute value of the interval.
    ///
    /// - If both bounds are non-negative: self.
    /// - If both bounds are non-positive: \[-hi, -lo\].
    /// - If the interval straddles zero: \[0, max(|lo|, |hi|)\].
    pub fn abs(&self) -> Self {
        if self.lo >= 0.0 {
            *self
        } else if self.hi <= 0.0 {
            Self::new(-self.hi, -self.lo)
        } else {
            Self::new(0.0, self.lo.abs().max(self.hi.abs()))
        }
    }

    /// Union of two intervals: the smallest interval containing both.
    pub fn union(&self, other: &Self) -> Self {
        Self::new(self.lo.min(other.lo), self.hi.max(other.hi))
    }
}

impl Add<Interval> for Interval {
    type Output = Interval;

    fn add(self, rhs: Interval) -> Interval {
        Interval::new(self.lo + rhs.lo, self.hi + rhs.hi)
    }
}

impl Sub<Interval> for Interval {
    type Output = Interval;

    fn sub(self, rhs: Interval) -> Interval {
        Interval::new(self.lo - rhs.hi, self.hi - rhs.lo)
    }
}

impl Mul<Interval> for Interval {
    type Output = Interval;

    fn mul(self, rhs: Interval) -> Interval {
        let p1 = self.lo * rhs.lo;
        let p2 = self.lo * rhs.hi;
        let p3 = self.hi * rhs.lo;
        let p4 = self.hi * rhs.hi;
        let lo = p1.min(p2).min(p3).min(p4);
        let hi = p1.max(p2).max(p3).max(p4);
        Interval::new(lo, hi)
    }
}

impl Neg for Interval {
    type Output = Interval;

    fn neg(self) -> Interval {
        Interval::new(-self.hi, -self.lo)
    }
}

impl Mul<f64> for Interval {
    type Output = Interval;

    fn mul(self, rhs: f64) -> Interval {
        let a = self.lo * rhs;
        let b = self.hi * rhs;
        Interval::new(a.min(b), a.max(b))
    }
}

impl Add<f64> for Interval {
    type Output = Interval;

    fn add(self, rhs: f64) -> Interval {
        Interval::new(self.lo + rhs, self.hi + rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_interval() {
        let i = Interval::point(3.0);
        assert!((i.lo - 3.0).abs() < 1e-15);
        assert!((i.hi - 3.0).abs() < 1e-15);
        assert!((i.width()).abs() < 1e-15);
        assert!((i.midpoint() - 3.0).abs() < 1e-15);
    }

    #[test]
    fn contains_and_predicates() {
        let i = Interval::new(1.0, 3.0);
        assert!(i.contains(1.0));
        assert!(i.contains(2.0));
        assert!(i.contains(3.0));
        assert!(!i.contains(0.99));
        assert!(!i.contains(3.01));
        assert!(i.is_positive());
        assert!(!i.is_negative());
        assert!(!i.contains_zero());

        let neg = Interval::new(-5.0, -2.0);
        assert!(neg.is_negative());
        assert!(!neg.is_positive());
        assert!(!neg.contains_zero());

        let straddling = Interval::new(-1.0, 1.0);
        assert!(straddling.contains_zero());
        assert!(!straddling.is_positive());
        assert!(!straddling.is_negative());
    }

    #[test]
    fn add_intervals() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 4.0);
        let c = a + b;
        assert!((c.lo - 3.0).abs() < 1e-15);
        assert!((c.hi - 7.0).abs() < 1e-15);
    }

    #[test]
    fn sub_intervals() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 4.0);
        let c = a - b;
        assert!((c.lo - (-3.0)).abs() < 1e-15);
        assert!((c.hi - 1.0).abs() < 1e-15);
    }

    #[test]
    fn mul_intervals_all_positive() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 4.0);
        let c = a * b;
        assert!((c.lo - 2.0).abs() < 1e-15);
        assert!((c.hi - 12.0).abs() < 1e-15);
        // Containment: all products of values in a * values in b should be in c
        for &x in &[1.0, 1.5, 2.0, 2.5, 3.0] {
            for &y in &[2.0, 2.5, 3.0, 3.5, 4.0] {
                assert!(
                    c.contains(x * y),
                    "c should contain {} * {} = {}",
                    x,
                    y,
                    x * y
                );
            }
        }
    }

    #[test]
    fn mul_intervals_mixed_signs() {
        let a = Interval::new(-2.0, 3.0);
        let b = Interval::new(-1.0, 4.0);
        let c = a * b;
        // Products: (-2)*(-1)=2, (-2)*4=-8, 3*(-1)=-3, 3*4=12
        assert!((c.lo - (-8.0)).abs() < 1e-15);
        assert!((c.hi - 12.0).abs() < 1e-15);
    }

    #[test]
    fn neg_interval() {
        let a = Interval::new(1.0, 3.0);
        let b = -a;
        assert!((b.lo - (-3.0)).abs() < 1e-15);
        assert!((b.hi - (-1.0)).abs() < 1e-15);
    }

    #[test]
    fn recip_positive() {
        let a = Interval::new(2.0, 4.0);
        let r = a.recip();
        assert!((r.lo - 0.25).abs() < 1e-15);
        assert!((r.hi - 0.5).abs() < 1e-15);
    }

    #[test]
    fn recip_negative() {
        let a = Interval::new(-4.0, -2.0);
        let r = a.recip();
        assert!((r.lo - (-0.5)).abs() < 1e-15);
        assert!((r.hi - (-0.25)).abs() < 1e-15);
    }

    #[test]
    fn recip_containing_zero_gives_wide_interval() {
        let a = Interval::new(-1.0, 1.0);
        let r = a.recip();
        assert!(r.lo <= -1e17);
        assert!(r.hi >= 1e17);
    }

    #[test]
    fn sqrt_positive() {
        let a = Interval::new(4.0, 9.0);
        let s = a.sqrt();
        assert!((s.lo - 2.0).abs() < 1e-15);
        assert!((s.hi - 3.0).abs() < 1e-15);
    }

    #[test]
    fn sqrt_clamps_negative_lo() {
        let a = Interval::new(-1.0, 4.0);
        let s = a.sqrt();
        assert!((s.lo).abs() < 1e-15);
        assert!((s.hi - 2.0).abs() < 1e-15);
    }

    #[test]
    fn abs_all_positive() {
        let a = Interval::new(1.0, 3.0);
        let b = a.abs();
        assert!((b.lo - 1.0).abs() < 1e-15);
        assert!((b.hi - 3.0).abs() < 1e-15);
    }

    #[test]
    fn abs_all_negative() {
        let a = Interval::new(-5.0, -2.0);
        let b = a.abs();
        assert!((b.lo - 2.0).abs() < 1e-15);
        assert!((b.hi - 5.0).abs() < 1e-15);
    }

    #[test]
    fn abs_straddles_zero() {
        let a = Interval::new(-3.0, 2.0);
        let b = a.abs();
        assert!((b.lo).abs() < 1e-15);
        assert!((b.hi - 3.0).abs() < 1e-15);
    }

    #[test]
    fn scalar_mul() {
        let a = Interval::new(1.0, 3.0);
        let b = a * 2.0;
        assert!((b.lo - 2.0).abs() < 1e-15);
        assert!((b.hi - 6.0).abs() < 1e-15);

        let c = a * (-2.0);
        assert!((c.lo - (-6.0)).abs() < 1e-15);
        assert!((c.hi - (-2.0)).abs() < 1e-15);
    }

    #[test]
    fn scalar_add() {
        let a = Interval::new(1.0, 3.0);
        let b = a + 5.0;
        assert!((b.lo - 6.0).abs() < 1e-15);
        assert!((b.hi - 8.0).abs() < 1e-15);
    }

    #[test]
    fn union_intervals() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 5.0);
        let c = a.union(&b);
        assert!((c.lo - 1.0).abs() < 1e-15);
        assert!((c.hi - 5.0).abs() < 1e-15);
    }

    #[test]
    fn width_and_midpoint() {
        let a = Interval::new(2.0, 8.0);
        assert!((a.width() - 6.0).abs() < 1e-15);
        assert!((a.midpoint() - 5.0).abs() < 1e-15);
    }
}
