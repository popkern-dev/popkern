//! Portable transcendental functions: the same bits on every target.
//!
//! `f64::sin` and friends call the platform's C math library: glibc on
//! x86_64 Linux, the `libm` port inside the Rust runtime on
//! `wasm32-unknown-unknown`. The two agree to an ulp or so, not bit
//! for bit — and a kernel decision fed by them (a refinement
//! threshold, a near-zero coordinate, a classification sign) can land
//! on opposite sides on two platforms. popkern live builds natively
//! for a self-hosted page and in wasm for the public page, so the same
//! script built different parts (PK-674: STEP coordinates differing in
//! the last digit, `-3.1e-33` against `0.`, and the propeller's blade
//! faces tessellating to different triangle counts).
//!
//! [`PortableMath`] routes every transcendental through the pure-Rust
//! [`libm`] crate (musl's algorithms, at most 1 ulp), compiled from the
//! same source on every target, so the result depends only on the
//! input bits. Arithmetic (`+ - * /`, `sqrt`, `mul_add`, `powi`,
//! rounding) is already IEEE-exact or target-independent and stays on
//! the inherent methods.
//!
//! Every method is the inherent method's name with a `p` (portable)
//! prefix: `theta.psin()`, `y.patan2(x)`. Library code must not call
//! the inherent transcendentals; `structural_portable_math` enforces it.
//!
//! ```
//! use popkern_math::PortableMath;
//! let (s, c) = 0.5f64.psin_cos();
//! assert_eq!(s, 0.5f64.psin());
//! assert!((s * s + c * c - 1.0).abs() < 1e-15);
//! ```

#![forbid(unsafe_code)]

/// Transcendental functions with the same bits on every target.
pub trait PortableMath: Sized + Copy {
    /// Portable sine.
    fn psin(self) -> Self;
    /// Portable cosine.
    fn pcos(self) -> Self;
    /// Portable tangent.
    fn ptan(self) -> Self;
    /// Portable arcsine.
    fn pasin(self) -> Self;
    /// Portable arccosine.
    fn pacos(self) -> Self;
    /// Portable arctangent.
    fn patan(self) -> Self;
    /// Portable hyperbolic sine.
    fn psinh(self) -> Self;
    /// Portable hyperbolic cosine.
    fn pcosh(self) -> Self;
    /// Portable hyperbolic tangent.
    fn ptanh(self) -> Self;
    /// Portable inverse hyperbolic sine.
    fn pasinh(self) -> Self;
    /// Portable inverse hyperbolic cosine.
    fn pacosh(self) -> Self;
    /// Portable inverse hyperbolic tangent.
    fn patanh(self) -> Self;
    /// Portable `e^self`.
    fn pexp(self) -> Self;
    /// Portable `e^self - 1`, accurate near zero.
    fn pexp_m1(self) -> Self;
    /// Portable natural logarithm.
    fn pln(self) -> Self;
    /// Portable `ln(1 + self)`, accurate near zero.
    fn pln_1p(self) -> Self;
    /// Portable base-2 logarithm.
    fn plog2(self) -> Self;
    /// Portable base-10 logarithm.
    fn plog10(self) -> Self;
    /// Portable cube root.
    fn pcbrt(self) -> Self;
    /// Portable `(sin, cos)` of `self`.
    fn psin_cos(self) -> (Self, Self);
    /// Portable four-quadrant arctangent of `self` (y) and `other` (x).
    fn patan2(self, other: Self) -> Self;
    /// Portable `sqrt(self² + other²)` without undue overflow.
    fn phypot(self, other: Self) -> Self;
    /// Portable `self` raised to the power `other`.
    fn ppowf(self, other: Self) -> Self;
}

impl PortableMath for f64 {
    #[inline]
    fn psin(self) -> Self {
        libm::sin(self)
    }
    #[inline]
    fn pcos(self) -> Self {
        libm::cos(self)
    }
    #[inline]
    fn ptan(self) -> Self {
        libm::tan(self)
    }
    #[inline]
    fn pasin(self) -> Self {
        libm::asin(self)
    }
    #[inline]
    fn pacos(self) -> Self {
        libm::acos(self)
    }
    #[inline]
    fn patan(self) -> Self {
        libm::atan(self)
    }
    #[inline]
    fn psinh(self) -> Self {
        libm::sinh(self)
    }
    #[inline]
    fn pcosh(self) -> Self {
        libm::cosh(self)
    }
    #[inline]
    fn ptanh(self) -> Self {
        libm::tanh(self)
    }
    #[inline]
    fn pasinh(self) -> Self {
        libm::asinh(self)
    }
    #[inline]
    fn pacosh(self) -> Self {
        libm::acosh(self)
    }
    #[inline]
    fn patanh(self) -> Self {
        libm::atanh(self)
    }
    #[inline]
    fn pexp(self) -> Self {
        libm::exp(self)
    }
    #[inline]
    fn pexp_m1(self) -> Self {
        libm::expm1(self)
    }
    #[inline]
    fn pln(self) -> Self {
        libm::log(self)
    }
    #[inline]
    fn pln_1p(self) -> Self {
        libm::log1p(self)
    }
    #[inline]
    fn plog2(self) -> Self {
        libm::log2(self)
    }
    #[inline]
    fn plog10(self) -> Self {
        libm::log10(self)
    }
    #[inline]
    fn pcbrt(self) -> Self {
        libm::cbrt(self)
    }
    #[inline]
    fn psin_cos(self) -> (Self, Self) {
        libm::sincos(self)
    }
    #[inline]
    fn patan2(self, other: Self) -> Self {
        libm::atan2(self, other)
    }
    #[inline]
    fn phypot(self, other: Self) -> Self {
        libm::hypot(self, other)
    }
    #[inline]
    fn ppowf(self, other: Self) -> Self {
        libm::pow(self, other)
    }
}

impl PortableMath for f32 {
    #[inline]
    fn psin(self) -> Self {
        libm::sinf(self)
    }
    #[inline]
    fn pcos(self) -> Self {
        libm::cosf(self)
    }
    #[inline]
    fn ptan(self) -> Self {
        libm::tanf(self)
    }
    #[inline]
    fn pasin(self) -> Self {
        libm::asinf(self)
    }
    #[inline]
    fn pacos(self) -> Self {
        libm::acosf(self)
    }
    #[inline]
    fn patan(self) -> Self {
        libm::atanf(self)
    }
    #[inline]
    fn psinh(self) -> Self {
        libm::sinhf(self)
    }
    #[inline]
    fn pcosh(self) -> Self {
        libm::coshf(self)
    }
    #[inline]
    fn ptanh(self) -> Self {
        libm::tanhf(self)
    }
    #[inline]
    fn pasinh(self) -> Self {
        libm::asinhf(self)
    }
    #[inline]
    fn pacosh(self) -> Self {
        libm::acoshf(self)
    }
    #[inline]
    fn patanh(self) -> Self {
        libm::atanhf(self)
    }
    #[inline]
    fn pexp(self) -> Self {
        libm::expf(self)
    }
    #[inline]
    fn pexp_m1(self) -> Self {
        libm::expm1f(self)
    }
    #[inline]
    fn pln(self) -> Self {
        libm::logf(self)
    }
    #[inline]
    fn pln_1p(self) -> Self {
        libm::log1pf(self)
    }
    #[inline]
    fn plog2(self) -> Self {
        libm::log2f(self)
    }
    #[inline]
    fn plog10(self) -> Self {
        libm::log10f(self)
    }
    #[inline]
    fn pcbrt(self) -> Self {
        libm::cbrtf(self)
    }
    #[inline]
    fn psin_cos(self) -> (Self, Self) {
        libm::sincosf(self)
    }
    #[inline]
    fn patan2(self, other: Self) -> Self {
        libm::atan2f(self, other)
    }
    #[inline]
    fn phypot(self, other: Self) -> Self {
        libm::hypotf(self, other)
    }
    #[inline]
    fn ppowf(self, other: Self) -> Self {
        libm::powf(self, other)
    }
}

#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
