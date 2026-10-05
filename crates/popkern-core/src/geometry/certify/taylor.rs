//! Truncated Taylor arithmetic in the seam parameter, with interval
//! coefficients.
//!
//! A [`Taylor<M>`] holds normalized coefficients `c_k = q⁽ᵏ⁾(s)/k!`,
//! `k < M`, each an enclosure over the same s-set. Evaluated at the
//! midpoint `m` of an s-interval `T` and over `T`, it supplies the
//! Taylor forms
//! `q(T) ⊆ Σ_{k<j} c_k(m)(T − m)^k + c_j(T)(T − m)^j`, `1 ≤ j < M`,
//! which the certifier intersects. On an exact seam the midpoint
//! coefficients are roundoff and the top form tightens as
//! `O(h^{M})` (PK-499 exploration E1/E2: the second-order form's
//! `O(h³)` constant made a full-circle rail cost ~4k intervals).
//!
//! Every recurrence below is inclusion-monotone, so coefficients
//! computed from enclosures over `T` enclose the true coefficients at
//! every `s ∈ T`.
#![allow(clippy::needless_range_loop)] // Coefficient recurrences index several series at once (`c_j · q_{k−j}`); index loops state them as written.

use super::interval::Iv;
use super::interval_vec::IvVec3;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Taylor order of the certifier's seam forms: coefficients `0..=4`.
pub const SEAM_TAYLOR_ORDER: usize = 4;

/// Coefficient count of a seam quantity series.
pub const SEAM_TERMS: usize = SEAM_TAYLOR_ORDER + 1;

/// A truncated Taylor series with `M` interval coefficients.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Taylor<const M: usize> {
    /// Normalized coefficients `q⁽ᵏ⁾/k!`.
    pub c: [Iv; M],
}

impl<const M: usize> Taylor<M> {
    /// The series with the given coefficients.
    pub fn new(c: [Iv; M]) -> Self {
        Self { c }
    }

    /// A constant.
    pub fn constant(v: Iv) -> Self {
        let mut c = [Iv::point(0.0); M];
        c[0] = v;
        Self { c }
    }

    /// The value (coefficient 0).
    pub fn value(&self) -> Iv {
        self.c[0]
    }

    /// Scale by an interval.
    pub fn scale(self, k: Iv) -> Self {
        Self {
            c: self.c.map(|x| x * k),
        }
    }

    /// `q²` (coefficient 0 via the tight square).
    pub fn sqr(self) -> Self {
        let mut out = self * self;
        out.c[0] = self.c[0].sqr();
        out
    }

    /// `√q`; every coefficient beyond 0 is [`Iv::ENTIRE`] when the
    /// value can reach 0.
    pub fn sqrt(self) -> Self {
        let mut r = [Iv::point(0.0); M];
        r[0] = self.c[0].sqrt();
        let twice = r[0].scale(2.0);
        for k in 1..M {
            let mut acc = self.c[k];
            for j in 1..k {
                acc = acc - r[j] * r[k - j];
            }
            r[k] = acc / twice;
        }
        Self { c: r }
    }

    /// `atan q`: `θ′ = q′/(1 + q²)` integrated term by term.
    pub fn atan(self) -> Self {
        let derivative = self.derivative_series();
        let one_plus = Self::constant(Iv::point(1.0)) + self.sqr();
        let rate = derivative / one_plus;
        let mut c = [Iv::point(0.0); M];
        c[0] = Iv::atan(self.c[0]);
        for k in 1..M {
            c[k] = rate.c[k - 1] / Iv::point(k as f64);
        }
        Self { c }
    }

    /// `dq/ds` as a series (its top coefficient is unknown and set to
    /// [`Iv::ENTIRE`]; callers integrate or truncate it away).
    pub fn derivative_series(self) -> Self {
        let mut c = [Iv::ENTIRE; M];
        for k in 0..M - 1 {
            c[k] = self.c[k + 1].scale((k + 1) as f64);
        }
        Self { c }
    }

    /// Integer power (`n ≥ 0`), truncated.
    pub fn powi(self, n: usize) -> Self {
        let mut out = Self::constant(Iv::point(1.0));
        for _ in 0..n {
            out = out * self;
        }
        out
    }

    /// The first `R` coefficients (`R ≤ M`).
    pub fn truncate<const R: usize>(self) -> Taylor<R> {
        Taylor {
            c: std::array::from_fn(|k| self.c[k]),
        }
    }

    /// Hull with another series.
    pub fn hull(self, o: Self) -> Self {
        Self {
            c: std::array::from_fn(|k| self.c[k].hull(o.c[k])),
        }
    }

    /// Enclosure of the function over `x₀ + offset` from its series at
    /// the point `x₀` (`self`) and over the whole set (`over`): the
    /// intersection of the natural enclosure `over.c[0]` and every
    /// Taylor form `Σ_{k<j} c_k(x₀) offset^k + over.c[j] offset^j`
    /// (non-finite forms are skipped).
    pub fn enclose_over(&self, over: &Self, offset: Iv) -> Iv {
        // offset^k, even powers through the tight square.
        let mut powers = [Iv::point(1.0); M];
        for k in 1..M {
            powers[k] = if k % 2 == 0 {
                powers[k / 2].sqr()
            } else {
                powers[k - 1] * offset
            };
        }
        let natural = over.c[0];
        let (mut lo, mut hi) = (natural.lo, natural.hi);
        let mut head = self.c[0];
        for j in 1..M {
            let form = head + over.c[j] * powers[j];
            if form.is_finite() {
                lo = lo.max(form.lo);
                hi = hi.min(form.hi);
            }
            head = head + self.c[j] * powers[j];
        }
        if lo <= hi { Iv::new(lo, hi) } else { natural }
    }

    /// Compose `self(σ)` with `σ(s) = σ₀ + Δ(s)`: `Σ_j c_j Δ^j`, where
    /// `delta` has zero constant term.
    pub fn compose(self, delta: Self) -> Self {
        let mut out = Self::constant(self.c[0]);
        let mut power = Self::constant(Iv::point(1.0));
        for j in 1..M {
            power = power * delta;
            out = out + power.scale(self.c[j]);
        }
        out
    }
}

impl<const M: usize> Add for Taylor<M> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            c: std::array::from_fn(|k| self.c[k] + o.c[k]),
        }
    }
}

impl<const M: usize> Sub for Taylor<M> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            c: std::array::from_fn(|k| self.c[k] - o.c[k]),
        }
    }
}

impl<const M: usize> Neg for Taylor<M> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            c: self.c.map(|x| -x),
        }
    }
}

impl<const M: usize> Mul for Taylor<M> {
    type Output = Self;
    /// Truncated Cauchy product.
    fn mul(self, o: Self) -> Self {
        let mut c = [Iv::point(0.0); M];
        for (k, slot) in c.iter_mut().enumerate() {
            let mut acc = self.c[0] * o.c[k];
            for j in 1..=k {
                acc = acc + self.c[j] * o.c[k - j];
            }
            *slot = acc;
        }
        Self { c }
    }
}

impl<const M: usize> Div for Taylor<M> {
    type Output = Self;
    /// `q_k = (a_k − Σ_{j≥1} b_j q_{k−j}) / b_0`; every coefficient is
    /// [`Iv::ENTIRE`] when `b_0` contains 0.
    fn div(self, o: Self) -> Self {
        let mut q = [Iv::point(0.0); M];
        for k in 0..M {
            let mut acc = self.c[k];
            for j in 1..=k {
                acc = acc - o.c[j] * q[k - j];
            }
            q[k] = acc / o.c[0];
        }
        Self { c: q }
    }
}

/// A 3-vector of [`Taylor`] series.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TaylorVec3<const M: usize> {
    /// x component.
    pub x: Taylor<M>,
    /// y component.
    pub y: Taylor<M>,
    /// z component.
    pub z: Taylor<M>,
}

impl<const M: usize> TaylorVec3<M> {
    /// From per-coefficient vector boxes.
    pub fn from_coefficients(c: [IvVec3; M]) -> Self {
        Self {
            x: Taylor::new(c.map(|v| v.x)),
            y: Taylor::new(c.map(|v| v.y)),
            z: Taylor::new(c.map(|v| v.z)),
        }
    }

    /// A constant vector.
    pub fn constant(v: IvVec3) -> Self {
        Self {
            x: Taylor::constant(v.x),
            y: Taylor::constant(v.y),
            z: Taylor::constant(v.z),
        }
    }

    /// Coefficient `k` as a vector box.
    pub fn coefficient(&self, k: usize) -> IvVec3 {
        IvVec3::new(self.x.c[k], self.y.c[k], self.z.c[k])
    }

    /// The value (coefficient 0).
    pub fn value(&self) -> IvVec3 {
        self.coefficient(0)
    }

    fn map(self, f: impl Fn(Taylor<M>) -> Taylor<M>) -> Self {
        Self {
            x: f(self.x),
            y: f(self.y),
            z: f(self.z),
        }
    }

    /// Scale by a scalar series.
    pub fn scale(self, k: Taylor<M>) -> Self {
        self.map(|c| c * k)
    }

    /// Dot product.
    pub fn dot(self, o: Self) -> Taylor<M> {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product.
    pub fn cross(self, o: Self) -> Self {
        Self {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    /// Length.
    pub fn norm(self) -> Taylor<M> {
        (self.x.sqr() + self.y.sqr() + self.z.sqr()).sqrt()
    }

    /// Derivative series (top coefficient unknown).
    pub fn derivative_series(self) -> Self {
        self.map(Taylor::derivative_series)
    }

    /// First `R` coefficients.
    pub fn truncate<const R: usize>(self) -> TaylorVec3<R> {
        TaylorVec3 {
            x: self.x.truncate(),
            y: self.y.truncate(),
            z: self.z.truncate(),
        }
    }

    /// Hull with another series.
    pub fn hull(self, o: Self) -> Self {
        Self {
            x: self.x.hull(o.x),
            y: self.y.hull(o.y),
            z: self.z.hull(o.z),
        }
    }

    /// Compose with `σ(s) = σ₀ + Δ(s)` (see [`Taylor::compose`]).
    pub fn compose(self, delta: Taylor<M>) -> Self {
        self.map(|c| c.compose(delta))
    }
}

impl<const M: usize> Add for TaylorVec3<M> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            x: self.x + o.x,
            y: self.y + o.y,
            z: self.z + o.z,
        }
    }
}

impl<const M: usize> Sub for TaylorVec3<M> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            x: self.x - o.x,
            y: self.y - o.y,
            z: self.z - o.z,
        }
    }
}
