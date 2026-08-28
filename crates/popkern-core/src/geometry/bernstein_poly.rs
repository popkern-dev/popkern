//! Polynomials in the Bernstein basis on `[0, 1]`.
//!
//! A polynomial of degree `n` is stored as its `n + 1` Bernstein
//! coefficients `(c_0, c_1, …, c_n)` with
//!
//! ```text
//!   p(t) = Σ_{i = 0..n} c_i · B_i^n(t)
//!   B_i^n(t) = C(n, i) · t^i · (1 − t)^{n − i}
//! ```
//!
//! This basis is numerically stable on `[0, 1]` and admits closed-
//! form coefficient-level formulas for the operations we need
//! (evaluation, derivative, degree elevation, polynomial product).
//!
//! The polynomial-identity G1/G2 solver (the polynomial-identity solver) uses this
//! module to represent cross-boundary tangent polynomials
//! `xtan(t)`, along-boundary tangents `T(t)`, and the scalar
//! blending polynomial `β(t)` that appears in the K-P continuity
//! identity `xtan_A(t) + α · xtan_B(t) + β(t) · T(t) ≡ 0`. Zeroing
//! every Bernstein coefficient of the identity's LHS gives **exact
//! G1 everywhere along the shared boundary** — not merely at
//! discrete samples, which is the structural limit of a sampled
//! solver.

use crate::geometry::Vec3;

/// Coefficient type admitted as the value of a Bernstein
/// polynomial. The solver uses `f64` for scalar blending
/// polynomials and [`Vec3`] for vector-valued tangents; the trait
/// is open so other coefficient types (matrices, affine maps, …)
/// can plug into the same operations without duplicating code.
pub trait PolyCoeff: Copy {
    /// The additive identity.
    fn zero() -> Self;
    /// Scalar multiplication `self · s`.
    fn scale(self, s: f64) -> Self;
    /// Addition `self + other`.
    fn add(self, other: Self) -> Self;
    /// Subtraction `self − other`.
    fn sub(self, other: Self) -> Self;
}

impl PolyCoeff for f64 {
    fn zero() -> Self {
        0.0
    }
    fn scale(self, s: f64) -> Self {
        self * s
    }
    fn add(self, other: Self) -> Self {
        self + other
    }
    fn sub(self, other: Self) -> Self {
        self - other
    }
}

impl PolyCoeff for Vec3 {
    fn zero() -> Self {
        Vec3::zero()
    }
    fn scale(self, s: f64) -> Self {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
    fn add(self, other: Self) -> Self {
        Vec3::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
    fn sub(self, other: Self) -> Self {
        Vec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

/// A polynomial of degree `len(coeffs) − 1` expressed in the
/// Bernstein basis on `[0, 1]`.
#[derive(Clone, Debug)]
pub struct BernsteinPoly<T: PolyCoeff> {
    coeffs: Vec<T>,
}

impl<T: PolyCoeff> BernsteinPoly<T> {
    /// Construct from Bernstein coefficients. Must be non-empty;
    /// a zero polynomial is degree 0 with the single coefficient
    /// [`PolyCoeff::zero`].
    pub fn new(coeffs: Vec<T>) -> Self {
        assert!(
            !coeffs.is_empty(),
            "BernsteinPoly requires at least one coefficient"
        );
        Self { coeffs }
    }

    /// Polynomial degree `n` (`= len − 1`).
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// The all-zero polynomial of a chosen degree.
    pub fn zero_of_degree(degree: usize) -> Self {
        Self {
            coeffs: vec![T::zero(); degree + 1],
        }
    }

    /// Access the `i`-th Bernstein coefficient.
    pub fn coeff(&self, i: usize) -> T {
        self.coeffs[i]
    }

    /// Access all Bernstein coefficients.
    pub fn coeffs(&self) -> &[T] {
        &self.coeffs
    }

    /// Evaluate `p(t)` via direct Bernstein-basis sum. Cost is
    /// `O(n)` after the basis table is formed; use [`de_casteljau`]
    /// (below) for marginally better numerical conditioning near
    /// the endpoints when `n` is large.
    pub fn evaluate(&self, t: f64) -> T {
        let n = self.degree();
        if n == 0 {
            return self.coeffs[0];
        }
        let mut acc = T::zero();
        for (i, &c) in self.coeffs.iter().enumerate() {
            acc = acc.add(c.scale(bernstein_basis(n, i, t)));
        }
        acc
    }

    /// Evaluate via the de Casteljau pyramid. Preferred for
    /// high-degree polynomials near endpoints where
    /// `t^i · (1 − t)^{n − i}` underflows the direct basis sum.
    pub fn de_casteljau(&self, t: f64) -> T {
        let mut work = self.coeffs.clone();
        let n = self.degree();
        let omt = 1.0 - t;
        for r in 1..=n {
            for i in 0..=(n - r) {
                work[i] = work[i].scale(omt).add(work[i + 1].scale(t));
            }
        }
        work[0]
    }

    /// Derivative `p'(t)` as a Bernstein polynomial of degree
    /// `n − 1`. Standard Bezier-curve formula
    /// `p'(t) = n · Σ_{i = 0..n−1} (c_{i+1} − c_i) · B_i^{n−1}(t)`.
    /// A degree-0 polynomial differentiates to the zero polynomial.
    pub fn differentiate(&self) -> Self {
        let n = self.degree();
        if n == 0 {
            return Self::zero_of_degree(0);
        }
        let n_f = n as f64;
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            out.push(self.coeffs[i + 1].sub(self.coeffs[i]).scale(n_f));
        }
        Self { coeffs: out }
    }

    /// Degree-elevate by `by` — return a polynomial of degree
    /// `n + by` whose point values `p(t)` at every `t ∈ [0, 1]`
    /// match the original. Needed when a Bernstein-product
    /// constraint mixes polynomials of different degrees and
    /// the solver wants them on a common degree.
    pub fn degree_elevate(self, by: usize) -> Self {
        let mut p = self;
        for _ in 0..by {
            p = p.elevate_once();
        }
        p
    }

    fn elevate_once(self) -> Self {
        let n = self.degree();
        let np1 = (n + 1) as f64;
        let mut out = Vec::with_capacity(n + 2);
        out.push(self.coeffs[0]);
        for i in 1..=n {
            let frac = (i as f64) / np1;
            let c = self.coeffs[i - 1]
                .scale(frac)
                .add(self.coeffs[i].scale(1.0 - frac));
            out.push(c);
        }
        out.push(self.coeffs[n]);
        Self { coeffs: out }
    }

    /// Convert to the **monomial (power) basis** on `[0, 1]`:
    /// returns coefficients `(a_0, a_1, …, a_n)` such that
    /// `p(t) = Σ a_k · t^k`.
    ///
    /// Formula (standard Bezier→power-basis conversion):
    ///
    /// ```text
    ///   a_k = Σ_{i = 0..=k} (−1)^{k − i} · C(n, i) · C(n − i, k − i) · c_i.
    /// ```
    ///
    /// Bernstein form is preferred for numerical work on `[0, 1]`,
    /// but the K-P bi-6 G¹/G² identity solver of §3.6 step 4 (see
    /// `docs/decisions/2026-06-22-kp15-bi6-deep-interior-without-maple.md`)
    /// needs the monomial coefficients of the identity polynomial
    /// in `u`. Zeroing each `a_k` (or equivalently each `c_k`) is
    /// the linear constraint the solver collects to determine the
    /// 16 deep-interior BB-coefficients per sector that the paper
    /// otherwise computes via a Maple script we cannot access.
    pub fn to_monomial_coefficients(&self) -> Vec<T> {
        let n = self.degree();
        let mut out = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let mut acc = T::zero();
            for (i, &c_i) in self.coeffs.iter().enumerate().take(k + 1) {
                let coeff = binomial(n, i) * binomial(n - i, k - i);
                let sign = if (k - i).is_multiple_of(2) { 1.0 } else { -1.0 };
                acc = acc.add(c_i.scale(sign * coeff));
            }
            out.push(acc);
        }
        out
    }

    /// Coefficient-wise add. If degrees differ the lower-degree
    /// polynomial is elevated to match before summing, so point
    /// values combine correctly.
    pub fn add_elevated(&self, other: &Self) -> Self {
        let (a, b) = match self.degree().cmp(&other.degree()) {
            std::cmp::Ordering::Equal => (self.clone(), other.clone()),
            std::cmp::Ordering::Less => (
                self.clone().degree_elevate(other.degree() - self.degree()),
                other.clone(),
            ),
            std::cmp::Ordering::Greater => (
                self.clone(),
                other.clone().degree_elevate(self.degree() - other.degree()),
            ),
        };
        let coeffs = a
            .coeffs
            .iter()
            .zip(b.coeffs.iter())
            .map(|(&x, &y)| x.add(y))
            .collect();
        Self { coeffs }
    }
}

/// `B_i^n(t) = C(n, i) · t^i · (1 − t)^{n − i}`.
pub(crate) fn bernstein_basis(n: usize, i: usize, t: f64) -> f64 {
    binomial(n, i) * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32)
}

/// Binomial coefficient `C(n, k)` as an `f64`. Stable multiplicative
/// form (no factorial overflow) — exact up to n ≈ 53 on IEEE 754.
pub(crate) fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut out = 1.0_f64;
    for i in 0..k {
        out *= (n - i) as f64 / (i + 1) as f64;
    }
    out
}

fn bernstein_product_generic<A, B, C>(
    a: &BernsteinPoly<A>,
    b: &BernsteinPoly<B>,
    mul: impl Fn(A, B) -> C,
) -> BernsteinPoly<C>
where
    A: PolyCoeff,
    B: PolyCoeff,
    C: PolyCoeff,
{
    let m = a.degree();
    let n = b.degree();
    let out_deg = m + n;
    let mut out = vec![C::zero(); out_deg + 1];
    let c_mn_k: Vec<f64> = (0..=out_deg).map(|k| binomial(m + n, k)).collect();
    for i in 0..=m {
        let ci_m = binomial(m, i);
        for j in 0..=n {
            let cj_n = binomial(n, j);
            let k = i + j;
            let weight = ci_m * cj_n / c_mn_k[k];
            let prod = mul(a.coeff(i), b.coeff(j));
            out[k] = out[k].add(prod.scale(weight));
        }
    }
    BernsteinPoly { coeffs: out }
}

/// Product of two scalar Bernstein polynomials via the Bezier
/// product formula: if `a(t) = Σ a_i B_i^m`, `b(t) = Σ b_j B_j^n`,
/// then `(a · b)(t) = Σ c_k B_k^{m+n}` with
///
/// ```text
///   c_k = Σ_{i + j = k} (C(m, i) · C(n, j) / C(m + n, k)) · a_i · b_j.
/// ```
pub fn bernstein_product_scalar(
    a: &BernsteinPoly<f64>,
    b: &BernsteinPoly<f64>,
) -> BernsteinPoly<f64> {
    bernstein_product_generic(a, b, |x, y| x * y)
}

/// Scalar · vector Bernstein product. Used to form `β(t) · T(t)`
/// where `β` is the K-P blending scalar polynomial and `T` is
/// the boundary tangent vector polynomial.
pub fn bernstein_product_scalar_vec3(
    a: &BernsteinPoly<f64>,
    b: &BernsteinPoly<Vec3>,
) -> BernsteinPoly<Vec3> {
    bernstein_product_generic(a, b, |x, y| y.scale(x))
}

#[cfg(test)]
#[path = "bernstein_poly_tests.rs"]
mod tests;
