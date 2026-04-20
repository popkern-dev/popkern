//! Degree-9 Bernstein basis + derivatives — hot-path evaluation
//! for bidegree-9 patches (Phase H).
//!
//! Tight, specialised implementation alongside the generic
//! [`BernsteinPoly`](crate::geometry::bernstein_poly::BernsteinPoly)
//! in `bernstein_poly`. The generic module covers polynomial
//! arithmetic (products, degree elevation, etc.); this module
//! handles raw scalar evaluation of `B_i^9(t)`, `B_i^9'(t)`,
//! `B_i^9''(t)` — the innermost loop of
//! [`RectangularGregoryBiNonic::point_at`](super::rectangular_binonic::RectangularGregoryBiNonic::point_at)
//! and its partial-derivative siblings (to be landed in H.2).
//!
//! # Source
//!
//! Standard Bernstein-basis identities (Farin 1993 Ch. 5):
//!
//! ```text
//!   B_i^n(t)    = C(n, i) · t^i · (1 − t)^{n − i}
//!   B_i^n'(t)   = n · (B_{i-1}^{n-1}(t) − B_i^{n-1}(t))
//!   B_i^n''(t)  = n(n − 1) · (B_{i-2}^{n-2}(t) − 2·B_{i-1}^{n-2}(t) + B_i^{n-2}(t))
//! ```
//!
//! At degree 9: `n = 9`, `n(n-1) = 72`. Out-of-range basis
//! values are treated as 0 (the usual convention that extends
//! the recursion naturally).
//!
//! # Performance
//!
//! Binomial tables for `C(9, ·)`, `C(8, ·)`, `C(7, ·)` are
//! `const` to avoid any allocation or branch on every call.
//! Evaluation uses `f64::powi` which is a single C-library
//! call (typically vectorisable); the whole `bernstein_binonic`
//! invocation is on the order of 4 multiplications + 1
//! indexed load, well within the per-corner budget.

/// Binomial coefficients `C(9, i)` for `i ∈ 0..=9`.
const C9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];

/// Binomial coefficients `C(8, i)` for `i ∈ 0..=8`.
/// Used in the first derivative.
const C8: [f64; 9] = [1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0];

/// Binomial coefficients `C(7, i)` for `i ∈ 0..=7`.
/// Used in the second derivative.
const C7: [f64; 8] = [1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0];

/// `B_i^9(t) = C(9, i) · t^i · (1 − t)^{9 − i}` for `i ∈ 0..=9`.
/// Returns `0.0` for `i` outside the range (the natural
/// extension of the Bernstein basis).
#[inline]
pub fn bernstein_binonic(i: usize, t: f64) -> f64 {
    if i > 9 {
        return 0.0;
    }
    let omt = 1.0 - t;
    C9[i] * t.powi(i as i32) * omt.powi((9 - i) as i32)
}

/// Degree-8 basis, `i ∈ 0..=8`. Helper for the first derivative.
/// `i` is `i32` here to allow negative inputs (from the
/// `i - 1` in the derivative recurrence).
#[inline]
fn bernstein_8(i: i32, t: f64) -> f64 {
    if !(0..=8).contains(&i) {
        return 0.0;
    }
    let idx = i as usize;
    let omt = 1.0 - t;
    C8[idx] * t.powi(i) * omt.powi(8 - i)
}

/// Degree-7 basis, `i ∈ 0..=7`. Helper for the second derivative.
#[inline]
fn bernstein_7(i: i32, t: f64) -> f64 {
    if !(0..=7).contains(&i) {
        return 0.0;
    }
    let idx = i as usize;
    let omt = 1.0 - t;
    C7[idx] * t.powi(i) * omt.powi(7 - i)
}

/// `B_i^9'(t) = 9 · (B_{i-1}^8(t) − B_i^8(t))`. Returns `0.0`
/// for `i > 9`.
#[inline]
pub fn bernstein_binonic_deriv(i: usize, t: f64) -> f64 {
    if i > 9 {
        return 0.0;
    }
    let ii = i as i32;
    9.0 * (bernstein_8(ii - 1, t) - bernstein_8(ii, t))
}

/// `B_i^9''(t) = 72 · (B_{i-2}^7(t) − 2·B_{i-1}^7(t) + B_i^7(t))`.
/// Returns `0.0` for `i > 9`.
#[inline]
pub fn bernstein_binonic_second_deriv(i: usize, t: f64) -> f64 {
    if i > 9 {
        return 0.0;
    }
    let ii = i as i32;
    72.0 * (bernstein_7(ii - 2, t) - 2.0 * bernstein_7(ii - 1, t) + bernstein_7(ii, t))
}

/// Convenience: evaluate all 10 basis functions `B_0^9(t)..B_9^9(t)`
/// at once into a fixed-size array. For hot paths that need
/// every basis value at a single `t` (the typical surface-
/// evaluation inner loop), this is faster than 10 separate
/// `bernstein_binonic` calls because `t.powi(·)` and
/// `(1-t).powi(·)` are reused.
#[inline]
pub fn bernstein_binonic_all(t: f64) -> [f64; 10] {
    let omt = 1.0 - t;
    // Build two tables: t_pow[i] = t^i, omt_pow[i] = (1-t)^i.
    let mut t_pow = [1.0_f64; 10];
    let mut omt_pow = [1.0_f64; 10];
    for i in 1..10 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    let mut out = [0.0_f64; 10];
    for i in 0..10 {
        out[i] = C9[i] * t_pow[i] * omt_pow[9 - i];
    }
    out
}

/// Convenience: evaluate all 10 derivative values
/// `B_0^9'(t)..B_9^9'(t)` at once.
#[inline]
pub fn bernstein_binonic_deriv_all(t: f64) -> [f64; 10] {
    let omt = 1.0 - t;
    let mut t_pow = [1.0_f64; 9];
    let mut omt_pow = [1.0_f64; 9];
    for i in 1..9 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    // b8[i] = B_i^8(t) for i = 0..=8.
    let mut b8 = [0.0_f64; 9];
    for i in 0..9 {
        b8[i] = C8[i] * t_pow[i] * omt_pow[8 - i];
    }
    // B_i^9'(t) = 9 · (B_{i-1}^8 − B_i^8), with out-of-range = 0.
    let mut out = [0.0_f64; 10];
    out[0] = 9.0 * (0.0 - b8[0]);
    for i in 1..9 {
        out[i] = 9.0 * (b8[i - 1] - b8[i]);
    }
    out[9] = 9.0 * (b8[8] - 0.0);
    out
}

/// Convenience: all 10 second-derivative values at once.
#[inline]
pub fn bernstein_binonic_second_deriv_all(t: f64) -> [f64; 10] {
    let omt = 1.0 - t;
    let mut t_pow = [1.0_f64; 8];
    let mut omt_pow = [1.0_f64; 8];
    for i in 1..8 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    let mut b7 = [0.0_f64; 8];
    for i in 0..8 {
        b7[i] = C7[i] * t_pow[i] * omt_pow[7 - i];
    }
    // B_i^9''(t) = 72 · (B_{i-2}^7 − 2·B_{i-1}^7 + B_i^7).
    // Out-of-range b7 terms are 0.
    let mut out = [0.0_f64; 10];
    let bm = |i: i32| -> f64 {
        if (0..=7).contains(&i) {
            b7[i as usize]
        } else {
            0.0
        }
    };
    for i in 0..10_i32 {
        let idx = i as usize;
        out[idx] = 72.0 * (bm(i - 2) - 2.0 * bm(i - 1) + bm(i));
    }
    out
}

#[cfg(test)]
#[path = "bernstein_binonic_tests.rs"]
mod tests;
