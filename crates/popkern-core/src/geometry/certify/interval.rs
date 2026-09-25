//! Outward-rounded interval arithmetic.
//!
//! Every IEEE-754 basic operation (`+ − × ÷ √`) is correctly rounded
//! to nearest, so the exact result lies within half an ulp of the
//! computed one; stepping the lower end one ulp down and the upper
//! end one ulp up therefore encloses it. Transcendentals come from
//! libm, which is not correctly rounded: their enclosures are widened
//! by [`TRIG_ULPS`] ulps plus [`TRIG_ABS_SLACK`] (the documented
//! libm assumption of the certified-seam design).
//!
//! Any operation whose result is not representable (an overflow, a
//! NaN, a division by an interval containing zero) yields
//! [`Iv::ENTIRE`], which consumers treat as inconclusive.

use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Ulps of widening applied to every libm transcendental.
pub const TRIG_ULPS: u32 = 4;

/// Absolute widening applied to every libm transcendental, covering
/// results that underflow to subnormals.
pub const TRIG_ABS_SLACK: f64 = 1e-300;

/// A closed interval `[lo, hi]` of reals, enclosing an unknown value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Iv {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}

fn down(x: f64) -> f64 {
    x.next_down()
}

fn up(x: f64) -> f64 {
    x.next_up()
}

impl Iv {
    /// The whole real line: the result of any unrepresentable
    /// operation.
    pub const ENTIRE: Iv = Iv {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    /// The degenerate interval `[x, x]` (exact: `x` is a float).
    pub fn point(x: f64) -> Iv {
        if x.is_nan() {
            Iv::ENTIRE
        } else {
            Iv { lo: x, hi: x }
        }
    }

    /// `[lo, hi]`, normalised so that a reversed or NaN pair becomes
    /// the hull or [`Iv::ENTIRE`].
    pub fn new(lo: f64, hi: f64) -> Iv {
        if lo.is_nan() || hi.is_nan() {
            Iv::ENTIRE
        } else {
            Iv {
                lo: lo.min(hi),
                hi: lo.max(hi),
            }
        }
    }

    /// Outward-rounded enclosure of a correctly rounded result pair.
    fn outward(lo: f64, hi: f64) -> Iv {
        if lo.is_nan() || hi.is_nan() {
            Iv::ENTIRE
        } else {
            Iv {
                lo: down(lo),
                hi: up(hi),
            }
        }
    }

    /// True when both ends are finite.
    pub fn is_finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }

    /// Width `hi − lo` (rounded up).
    pub fn width(self) -> f64 {
        up(self.hi - self.lo)
    }

    /// Midpoint (a float inside the interval when it is finite).
    pub fn mid(self) -> f64 {
        let m = 0.5 * self.lo + 0.5 * self.hi;
        m.clamp(self.lo, self.hi)
    }

    /// Largest absolute value in the interval.
    pub fn mag(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }

    /// Smallest absolute value in the interval (0 when it contains 0).
    pub fn mig(self) -> f64 {
        if self.contains_zero() {
            0.0
        } else {
            self.lo.abs().min(self.hi.abs())
        }
    }

    /// True when `x` lies in the interval.
    pub fn contains(self, x: f64) -> bool {
        self.lo <= x && x <= self.hi
    }

    /// True when 0 lies in the interval.
    pub fn contains_zero(self) -> bool {
        self.contains(0.0)
    }

    /// Smallest interval containing both.
    pub fn hull(self, other: Iv) -> Iv {
        Iv {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    /// `x²`, tighter than `x · x` when the interval contains 0.
    pub fn sqr(self) -> Iv {
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        if self.contains_zero() {
            Iv {
                lo: 0.0,
                hi: up(a.max(b)),
            }
        } else {
            let out = Iv::outward(a.min(b), a.max(b));
            Iv {
                lo: out.lo.max(0.0),
                hi: out.hi,
            }
        }
    }

    /// `√x` over the non-negative part. A caller passes an interval
    /// whose true value is known to be non-negative (a sum of
    /// squares); a wholly negative interval yields [`Iv::ENTIRE`].
    pub fn sqrt(self) -> Iv {
        if self.hi < 0.0 || self.lo.is_nan() || self.hi.is_nan() {
            return Iv::ENTIRE;
        }
        let lo = if self.lo <= 0.0 {
            0.0
        } else {
            down(self.lo.sqrt()).max(0.0)
        };
        Iv {
            lo,
            hi: up(self.hi.sqrt()),
        }
    }

    /// `|x|`.
    pub fn abs(self) -> Iv {
        if self.contains_zero() {
            Iv {
                lo: 0.0,
                hi: self.mag(),
            }
        } else {
            Iv {
                lo: self.mig(),
                hi: self.mag(),
            }
        }
    }

    /// `1 / x`; [`Iv::ENTIRE`] when the interval contains 0.
    pub fn recip(self) -> Iv {
        if self.contains_zero() {
            return Iv::ENTIRE;
        }
        let (a, b) = (1.0 / self.hi, 1.0 / self.lo);
        Iv::outward(a.min(b), a.max(b))
    }

    /// `k · x` for an exact float `k`.
    pub fn scale(self, k: f64) -> Iv {
        self * Iv::point(k)
    }

    fn from_candidates(values: &[f64]) -> Iv {
        if values.iter().any(|x| x.is_nan()) {
            return Iv::ENTIRE;
        }
        let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Iv::outward(lo, hi)
    }

    /// Widen a libm result by the documented transcendental slack.
    fn widen_libm(lo: f64, hi: f64) -> Iv {
        let (mut lo, mut hi) = (lo, hi);
        for _ in 0..TRIG_ULPS {
            lo = down(lo);
            hi = up(hi);
        }
        Iv::new(lo - TRIG_ABS_SLACK, hi + TRIG_ABS_SLACK)
    }

    /// `sin x`, enclosing every extremum the interval may contain.
    pub fn sin(self) -> Iv {
        // sin x = cos(x − π/2); the shift is exact enough only as an
        // interval, so do it outward-rounded.
        (self - Iv::point(FRAC_PI_2).hull(Iv::point(FRAC_PI_2.next_up()))).cos()
    }

    /// `cos x`, enclosing every extremum the interval may contain.
    ///
    /// The maxima of `cos` lie at `2kπ` and the minima at `(2k+1)π`.
    /// An extremum is included whenever it could lie in the interval
    /// given the rounding of `k·π` (a slack of `1e-12·(1 + |x|)`),
    /// which may widen but never narrows the enclosure.
    pub fn cos(self) -> Iv {
        if !self.is_finite() || self.width() >= TAU {
            return Iv::new(-1.0, 1.0);
        }
        let (a, b) = (self.lo.cos(), self.hi.cos());
        let mut enc = Iv::widen_libm(a.min(b), a.max(b));
        let slack = 1e-12 * (1.0 + self.mag());
        let first = ((self.lo - slack) / PI).ceil() as i64;
        let last = ((self.hi + slack) / PI).floor() as i64;
        for k in first..=last {
            if k.rem_euclid(2) == 0 {
                enc.hi = 1.0;
            } else {
                enc.lo = -1.0;
            }
        }
        Iv::new(enc.lo.max(-1.0), enc.hi.min(1.0))
    }

    /// `atan x` (monotone increasing).
    pub fn atan(self) -> Iv {
        if self.lo.is_nan() || self.hi.is_nan() {
            return Iv::ENTIRE;
        }
        let enc = Iv::widen_libm(self.lo.atan(), self.hi.atan());
        Iv::new(
            enc.lo.max(-FRAC_PI_2.next_up()),
            enc.hi.min(FRAC_PI_2.next_up()),
        )
    }

    /// `tan x` for a point argument (used on constant half-angles).
    pub fn tan_point(x: f64) -> Iv {
        let t = x.tan();
        Iv::widen_libm(t, t)
    }
}

impl Add for Iv {
    type Output = Iv;
    fn add(self, rhs: Iv) -> Iv {
        Iv::outward(self.lo + rhs.lo, self.hi + rhs.hi)
    }
}

impl Sub for Iv {
    type Output = Iv;
    fn sub(self, rhs: Iv) -> Iv {
        Iv::outward(self.lo - rhs.hi, self.hi - rhs.lo)
    }
}

impl Neg for Iv {
    type Output = Iv;
    fn neg(self) -> Iv {
        Iv {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}

impl Mul for Iv {
    type Output = Iv;
    fn mul(self, rhs: Iv) -> Iv {
        // 0 · ∞ is NaN; an infinite end only arises from ENTIRE, so
        // the product is ENTIRE unless one factor is exactly [0, 0].
        if self == Iv::point(0.0) || rhs == Iv::point(0.0) {
            return Iv::point(0.0);
        }
        Iv::from_candidates(&[
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        ])
    }
}

impl Div for Iv {
    type Output = Iv;
    /// `self / other`; [`Iv::ENTIRE`] when `other` contains 0.
    fn div(self, other: Iv) -> Iv {
        if other.contains_zero() {
            return Iv::ENTIRE;
        }
        let quotients = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ];
        Iv::from_candidates(&quotients)
    }
}
