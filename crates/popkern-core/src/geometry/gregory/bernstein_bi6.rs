//! Phase I — Degree-6 Bernstein basis + derivatives through 3rd
//! order.
//!
//! Shared primitive for the K-P 2015 minimal-bi-6 polar
//! construction (`docs/journal/polar-kp-math.md` §2). Follows the
//! same pattern as
//! [`bernstein_binonic`](super::bernstein_binonic) and
//! [`bernstein_bi5`](super::bernstein_bi5) — specialised scalar
//! tables, no heap allocation.
//!
//! # Identities
//!
//! ```text
//!   B_i^6(t)    = C(6, i) · t^i · (1 − t)^{6 − i}
//!   B_i^6'(t)   =  6 · (B_{i-1}^5(t) − B_i^5(t))
//!   B_i^6''(t)  = 30 · (B_{i-2}^4(t) − 2·B_{i-1}^4(t) + B_i^4(t))
//!   B_i^6'''(t) = 120 · (B_{i-3}^3(t) − 3·B_{i-2}^3(t)
//!                         + 3·B_{i-1}^3(t) − B_i^3(t))
//! ```
//!
//! where `n · (n − 1) = 30`, `n · (n − 1) · (n − 2) = 120`. Basis
//! values outside `0..=n` are 0 (standard convention). 3rd
//! derivatives support F_5 quality-functional gradients in the K-P
//! 2015 step 4.
//!
//! # Usage footprint
//!
//! 7 basis values per `t` (indices 0..=6). A bi-6 tensor-product
//! patch has 7×7 = 49 control points, so a single
//! `point_at(u, v)` involves 49 CP × 1 scalar weight each and
//! partials involve a second 7-value table. Evaluation cost is
//! comparable to bi-5 (36 CPs) and cheaper than bi-9 (100 CPs).

/// Binomial coefficients `C(6, i)` for `i ∈ 0..=6`.
const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];

/// Binomial coefficients `C(5, i)` for `i ∈ 0..=5`.
const C5: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];

/// Binomial coefficients `C(4, i)` for `i ∈ 0..=4`.
const C4: [f64; 5] = [1.0, 4.0, 6.0, 4.0, 1.0];

/// Binomial coefficients `C(3, i)` for `i ∈ 0..=3`.
const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];

/// `B_i^6(t)` for `i ∈ 0..=6`. Returns `0.0` for `i > 6`.
#[inline]
pub fn bernstein_bi6(i: usize, t: f64) -> f64 {
    if i > 6 {
        return 0.0;
    }
    let omt = 1.0 - t;
    C6[i] * t.powi(i as i32) * omt.powi((6 - i) as i32)
}

/// `B_i^6'(t) = 6 · (B_{i-1}^5(t) − B_i^5(t))`.
#[inline]
pub fn bernstein_bi6_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    6.0 * (b5(ii - 1, t) - b5(ii, t))
}

/// `B_i^6''(t) = 30 · (B_{i-2}^4(t) − 2·B_{i-1}^4(t) + B_i^4(t))`.
#[inline]
pub fn bernstein_bi6_second_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    30.0 * (b4(ii - 2, t) - 2.0 * b4(ii - 1, t) + b4(ii, t))
}

/// `B_i^6'''(t) = 120 · (B_{i-3}^3(t) − 3·B_{i-2}^3(t) + 3·B_{i-1}^3(t) − B_i^3(t))`.
#[inline]
pub fn bernstein_bi6_third_deriv(i: usize, t: f64) -> f64 {
    let ii = i as i32;
    120.0 * (b3(ii - 3, t) - 3.0 * b3(ii - 2, t) + 3.0 * b3(ii - 1, t) - b3(ii, t))
}

/// All 7 basis values `B_0^6(t)..B_6^6(t)` at once.
#[inline]
pub fn bernstein_bi6_all(t: f64) -> [f64; 7] {
    let omt = 1.0 - t;
    let mut t_pow = [1.0_f64; 7];
    let mut omt_pow = [1.0_f64; 7];
    for i in 1..7 {
        t_pow[i] = t_pow[i - 1] * t;
        omt_pow[i] = omt_pow[i - 1] * omt;
    }
    let mut out = [0.0_f64; 7];
    for i in 0..7 {
        out[i] = C6[i] * t_pow[i] * omt_pow[6 - i];
    }
    out
}

/// All 7 first-derivative values at once.
#[inline]
pub fn bernstein_bi6_deriv_all(t: f64) -> [f64; 7] {
    let b = bernstein_b5_all(t);
    let mut out = [0.0_f64; 7];
    out[0] = 6.0 * (0.0 - b[0]);
    for i in 1..6 {
        out[i] = 6.0 * (b[i - 1] - b[i]);
    }
    out[6] = 6.0 * (b[5] - 0.0);
    out
}

/// All 7 second-derivative values at once.
#[inline]
pub fn bernstein_bi6_second_deriv_all(t: f64) -> [f64; 7] {
    let b = bernstein_b4_all(t);
    let bm = |i: i32| -> f64 {
        if (0..=4).contains(&i) {
            b[i as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 7];
    for i in 0..7_i32 {
        out[i as usize] = 30.0 * (bm(i - 2) - 2.0 * bm(i - 1) + bm(i));
    }
    out
}

/// All 7 third-derivative values at once.
#[inline]
pub fn bernstein_bi6_third_deriv_all(t: f64) -> [f64; 7] {
    let b = bernstein_b3_all(t);
    let bm = |i: i32| -> f64 {
        if (0..=3).contains(&i) {
            b[i as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 7];
    for i in 0..7_i32 {
        out[i as usize] = 120.0 * (bm(i - 3) - 3.0 * bm(i - 2) + 3.0 * bm(i - 1) - bm(i));
    }
    out
}

#[inline]
fn bernstein_b5_all(t: f64) -> [f64; 6] {
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

#[inline]
fn bernstein_b4_all(t: f64) -> [f64; 5] {
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

#[inline]
fn bernstein_b3_all(t: f64) -> [f64; 4] {
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

#[inline]
fn b5(i: i32, t: f64) -> f64 {
    if !(0..=5).contains(&i) {
        return 0.0;
    }
    let idx = i as usize;
    let omt = 1.0 - t;
    C5[idx] * t.powi(i) * omt.powi(5 - i)
}

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

#[cfg(test)]
#[path = "bernstein_bi6_tests.rs"]
mod tests;
