//! Phase I — Degree-5 Bernstein basis + derivatives through 3rd
//! order.
//!
//! Shared primitive for the K-P 2013 biquintic polar construction
//! (`docs/journal/polar-kp-math.md` §3) and the pre-existing
//! biquintic guide patch in
//! [`kp_guide_patch_biquintic`](super::kp_guide_patch_biquintic).
//!
//! Follows the pattern established by
//! [`bernstein_binonic`](super::bernstein_binonic) — specialised
//! scalar tables, no heap allocation, bit-identical outputs across
//! repeated invocations.
//!
//! # Identities
//!
//! ```text
//!   B_i^5(t)    = C(5, i) · t^i · (1 − t)^{5 − i}
//!   B_i^5'(t)   = 5 · (B_{i-1}^4(t) − B_i^4(t))
//!   B_i^5''(t)  = 20 · (B_{i-2}^3(t) − 2·B_{i-1}^3(t) + B_i^3(t))
//!   B_i^5'''(t) = 60 · (B_{i-3}^2(t) − 3·B_{i-2}^2(t)
//!                        + 3·B_{i-1}^2(t) − B_i^2(t))
//! ```
//!
//! where `n · (n − 1) = 20`, `n · (n − 1) · (n − 2) = 60`. Basis
//! values outside `0..=n` are taken as 0 (standard convention).
//! 3rd derivatives are needed for F_3 quality-functional gradients
//! in the K-P 2013 σ step.

/// Binomial coefficients `C(5, i)` for `i ∈ 0..=5`.
const C5: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];

/// Binomial coefficients `C(4, i)` for `i ∈ 0..=4`.
const C4: [f64; 5] = [1.0, 4.0, 6.0, 4.0, 1.0];

/// Binomial coefficients `C(3, i)` for `i ∈ 0..=3`.
const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];

/// Binomial coefficients `C(2, i)` for `i ∈ 0..=2`.
const C2: [f64; 3] = [1.0, 2.0, 1.0];

/// `B_i^5(t)` for `i ∈ 0..=5`. Returns `0.0` for `i > 5`.
#[inline]
pub fn bernstein_bi5(i: usize, t: f64) -> f64 {
    if i > 5 {
        return 0.0;
    }
    let omt = 1.0 - t;
    C5[i] * t.powi(i as i32) * omt.powi((5 - i) as i32)
}

/// `B_i^5'(t) = 5 · (B_{i-1}^4(t) − B_i^4(t))`.
#[inline]
pub fn bernstein_bi5_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    5.0 * (b4(ii - 1, t) - b4(ii, t))
}

/// `B_i^5''(t) = 20 · (B_{i-2}^3(t) − 2·B_{i-1}^3(t) + B_i^3(t))`.
#[inline]
pub fn bernstein_bi5_second_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    20.0 * (b3(ii - 2, t) - 2.0 * b3(ii - 1, t) + b3(ii, t))
}

/// `B_i^5'''(t) = 60 · (B_{i-3}^2(t) − 3·B_{i-2}^2(t) + 3·B_{i-1}^2(t) − B_i^2(t))`.
#[inline]
pub fn bernstein_bi5_third_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    60.0 * (b2(ii - 3, t) - 3.0 * b2(ii - 2, t) + 3.0 * b2(ii - 1, t) - b2(ii, t))
}

/// All 6 basis values `B_0^5(t)..B_5^5(t)` at once — faster than 6
/// separate calls because `t^i` and `(1-t)^i` are reused.
#[inline]
pub fn bernstein_bi5_all(t: f64) -> [f64; 6] {
    let omt = 1.0 - t;
    let mut t_pow = [1.0_f64; 6];
    let mut omt_pow = [1.0_f64; 6];
    for i in 1..6 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    let mut out = [0.0_f64; 6];
    for i in 0..6 {
        out[i] = C5[i] * t_pow[i] * omt_pow[5 - i];
    }
    out
}

/// All 6 first-derivative values at once.
#[inline]
pub fn bernstein_bi5_deriv_all(t: f64) -> [f64; 6] {
    let b = bernstein_bi4_all(t);
    let mut out = [0.0_f64; 6];
    out[0] = 5.0 * (0.0 - b[0]);
    for i in 1..5 {
        out[i] = 5.0 * (b[i - 1] - b[i]);
    }
    out[5] = 5.0 * (b[4] - 0.0);
    out
}

/// All 6 second-derivative values at once.
#[inline]
pub fn bernstein_bi5_second_deriv_all(t: f64) -> [f64; 6] {
    let b = super::bernstein_bi3::bernstein_bi3_all(t);
    let bm = |i: i32| -> f64 {
        if (0..=3).contains(&i) {
            b[i as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 6];
    for i in 0..6_i32 {
        out[i as usize] = 20.0 * (bm(i - 2) - 2.0 * bm(i - 1) + bm(i));
    }
    out
}

/// All 6 third-derivative values at once.
#[inline]
pub fn bernstein_bi5_third_deriv_all(t: f64) -> [f64; 6] {
    let b = bernstein_bi2_all(t);
    let bm = |i: i32| -> f64 {
        if (0..=2).contains(&i) {
            b[i as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 6];
    for i in 0..6_i32 {
        out[i as usize] = 60.0 * (bm(i - 3) - 3.0 * bm(i - 2) + 3.0 * bm(i - 1) - bm(i));
    }
    out
}

/// All 5 basis values `B_0^4(t)..B_4^4(t)` at once.
#[inline]
pub fn bernstein_bi4_all(t: f64) -> [f64; 5] {
    let omt = 1.0 - t;
    let mut t_pow = [1.0_f64; 5];
    let mut omt_pow = [1.0_f64; 5];
    for i in 1..5 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    let mut out = [0.0_f64; 5];
    for i in 0..5 {
        out[i] = C4[i] * t_pow[i] * omt_pow[4 - i];
    }
    out
}

// `bernstein_bi3_all` now lives in the dedicated `bernstein_bi3`
// module (Phase IV PnS support); imported below as needed.

/// All 3 basis values `B_0^2(t)..B_2^2(t)` at once.
#[inline]
pub fn bernstein_bi2_all(t: f64) -> [f64; 3] {
    let omt = 1.0 - t;
    [C2[0] * omt * omt, C2[1] * t * omt, C2[2] * t * t]
}

/// Scalar helper `B_i^4(t)` for `i : i32` (negative / out-of-range
/// → 0.0). Used by the derivative-recurrence formulas; not
/// typically called directly.
#[inline]
fn b4(i: i32, t: f64) -> f64 {
    if !(0..=4).contains(&i) {
        return 0.0;
    }
    let idx = i as usize;
    let omt = 1.0 - t;
    C4[idx] * t.powi(i) * omt.powi(4 - i)
}

#[inline]
fn b3(i: i32, t: f64) -> f64 {
    if !(0..=3).contains(&i) {
        return 0.0;
    }
    let idx = i as usize;
    let omt = 1.0 - t;
    C3[idx] * t.powi(i) * omt.powi(3 - i)
}

#[inline]
fn b2(i: i32, t: f64) -> f64 {
    if !(0..=2).contains(&i) {
        return 0.0;
    }
    let idx = i as usize;
    let omt = 1.0 - t;
    C2[idx] * t.powi(i) * omt.powi(2 - i)
}

#[cfg(test)]
#[path = "bernstein_bi5_tests.rs"]
mod tests;
