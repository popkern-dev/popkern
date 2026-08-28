//! Degree-3 Bernstein basis + derivatives through 3rd
//! order. The tensor-product of two instances produces a bi-cubic
//! (bi-3) patch — the output primitive of the Polyhedral-Net Splines
//! (PnS) n-gon construction (`docs/papers/14bi234.pdf` §4.2,
//! Appendix C) and used for 2×2-per-sector splits at higher valence.
//!
//! Follows the pattern of [`bernstein_bi5`](super::bernstein_bi5) /
//! [`bernstein_bi6`](super::bernstein_bi6) /
//! [`bernstein_binonic`](super::bernstein_binonic): specialised
//! scalar tables, no heap allocation, bit-identical outputs across
//! repeated invocations.
//!
//! # Identities
//!
//! ```text
//!   B_i^3(t)    = C(3, i) · t^i · (1 − t)^{3 − i}
//!   B_i^3'(t)   = 3 · (B_{i-1}^2(t) − B_i^2(t))
//!   B_i^3''(t)  = 6 · (B_{i-2}^1(t) − 2·B_{i-1}^1(t) + B_i^1(t))
//!   B_i^3'''(t) = 6 · (B_{i-3}^0(t) − 3·B_{i-2}^0(t)
//!                        + 3·B_{i-1}^0(t) − B_i^0(t))
//! ```
//!
//! where `3! = 6`, `3 · 2 = 6`. Basis values outside `0..=n` are
//! taken as 0 (standard convention).

/// Binomial coefficients `C(3, i)` for `i ∈ 0..=3`.
const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];

/// Binomial coefficients `C(2, i)` for `i ∈ 0..=2`.
const C2: [f64; 3] = [1.0, 2.0, 1.0];

/// `B_i^3(t)` for `i ∈ 0..=3`. Returns `0.0` for `i > 3`.
#[inline]
pub fn bernstein_bi3(i: usize, t: f64) -> f64 {
    if i > 3 {
        return 0.0;
    }
    let omt = 1.0 - t;
    C3[i] * t.powi(i as i32) * omt.powi((3 - i) as i32)
}

/// `B_i^3'(t) = 3 · (B_{i-1}^2(t) − B_i^2(t))`.
#[inline]
pub fn bernstein_bi3_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    3.0 * (b2(ii - 1, t) - b2(ii, t))
}

/// `B_i^3''(t) = 6 · (B_{i-2}^1(t) − 2·B_{i-1}^1(t) + B_i^1(t))`.
#[inline]
pub fn bernstein_bi3_second_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    6.0 * (b1(ii - 2, t) - 2.0 * b1(ii - 1, t) + b1(ii, t))
}

/// `B_i^3'''(t) = 6 · (B_{i-3}^0(t) − 3·B_{i-2}^0(t) + 3·B_{i-1}^0(t) − B_i^0(t))`.
/// At degree 3 the third derivative is a piecewise-constant on each
/// basis function (degree-0 step functions): `B_j^0(t) = 1` for
/// `j == 0`, else `0`. Useful for F_3 quality-functional gradients.
#[inline]
pub fn bernstein_bi3_third_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    6.0 * (b0(ii - 3, t) - 3.0 * b0(ii - 2, t) + 3.0 * b0(ii - 1, t) - b0(ii, t))
}

/// All 4 basis values `B_0^3(t)..B_3^3(t)` at once — faster than 4
/// separate calls because `t^i` and `(1-t)^i` are reused.
#[inline]
pub fn bernstein_bi3_all(t: f64) -> [f64; 4] {
    let omt = 1.0 - t;
    let mut t_pow = [1.0_f64; 4];
    let mut omt_pow = [1.0_f64; 4];
    for i in 1..4 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    let mut out = [0.0_f64; 4];
    for i in 0..4 {
        out[i] = C3[i] * t_pow[i] * omt_pow[3 - i];
    }
    out
}

/// All 4 first-derivative values `B_0^3'(t)..B_3^3'(t)` at once.
#[inline]
pub fn bernstein_bi3_deriv_all(t: f64) -> [f64; 4] {
    let b2_arr = bernstein_bi2_all(t);
    [
        3.0 * (0.0 - b2_arr[0]),
        3.0 * (b2_arr[0] - b2_arr[1]),
        3.0 * (b2_arr[1] - b2_arr[2]),
        3.0 * (b2_arr[2] - 0.0),
    ]
}

/// All 4 second-derivative values `B_0^3''(t)..B_3^3''(t)` at once.
#[inline]
pub fn bernstein_bi3_second_deriv_all(t: f64) -> [f64; 4] {
    let b1_arr = bernstein_bi1_all(t);
    [
        6.0 * b1_arr[0],
        6.0 * (-2.0 * b1_arr[0] + b1_arr[1]),
        6.0 * (b1_arr[0] - 2.0 * b1_arr[1]),
        6.0 * b1_arr[1],
    ]
}

// ── Helpers: lower-degree basis values for the recursion tables ──

/// `B_i^2(t)` for `i ∈ 0..=2`, else `0`.
#[inline]
fn b2(i: i32, t: f64) -> f64 {
    if !(0..=2).contains(&i) {
        return 0.0;
    }
    let omt = 1.0 - t;
    C2[i as usize] * t.powi(i) * omt.powi(2 - i)
}

/// `B_i^1(t)` for `i ∈ {0, 1}`, else `0`.
#[inline]
fn b1(i: i32, t: f64) -> f64 {
    match i {
        0 => 1.0 - t,
        1 => t,
        _ => 0.0,
    }
}

/// `B_i^0(t) = 1` if `i == 0`, else `0` (the constant basis).
#[inline]
fn b0(i: i32, _t: f64) -> f64 {
    if i == 0 { 1.0 } else { 0.0 }
}

/// All 3 degree-2 basis values.
#[inline]
fn bernstein_bi2_all(t: f64) -> [f64; 3] {
    let omt = 1.0 - t;
    [omt * omt, 2.0 * t * omt, t * t]
}

/// All 2 degree-1 basis values.
#[inline]
fn bernstein_bi1_all(t: f64) -> [f64; 2] {
    [1.0 - t, t]
}

#[cfg(test)]
#[path = "bernstein_bi3_tests.rs"]
mod tests;
