//! KP20 support — triangular total-degree-6 Bernstein-Bézier basis
//! and partials through third order.
//!
//! KP20 uses a piecewise total-degree-6 guide surface `g` on each
//! sector wedge. This module mirrors [`bernstein_tri5`](super::bernstein_tri5)
//! but keeps the implementation generic enough to derive first,
//! second, and third partials from the standard barycentric
//! differentiation identities.

/// Number of CPs at total degree 6: `(6 + 1)(6 + 2)/2 = 28`.
pub const TRI6_CP_COUNT: usize = 28;

/// Flat index for `(i, j)` at degree `n`, where the third
/// barycentric index is `k = n - i - j`.
#[inline]
pub const fn flat_index(n: usize, i: usize, j: usize) -> usize {
    assert!(i + j <= n);
    i * (2 * n + 3 - i) / 2 + j
}

#[inline]
fn tri_cp_count(n: usize) -> usize {
    (n + 1) * (n + 2) / 2
}

/// Evaluate all degree-6 triangular Bernstein basis functions at
/// `(u, v, w = 1 - u - v)`.
pub fn bernstein_tri6_all(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(eval_basis(6, u, v))
}

/// Evaluate all degree-6 basis functions without a simplex-domain
/// assumption. Same polynomial as [`bernstein_tri6_all`].
pub fn bernstein_tri6_all_extended(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    bernstein_tri6_all(u, v)
}

/// `∂/∂u` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_partial_u(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 1, 0, u, v))
}

/// `∂/∂v` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_partial_v(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 0, 1, u, v))
}

/// `∂²/∂u²` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_second_partial_uu(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 2, 0, u, v))
}

/// `∂²/∂u∂v` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_second_partial_uv(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 1, 1, u, v))
}

/// `∂²/∂v²` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_second_partial_vv(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 0, 2, u, v))
}

/// `∂³/∂u³` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_third_partial_uuu(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 3, 0, u, v))
}

/// `∂³/∂u²∂v` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_third_partial_uuv(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 2, 1, u, v))
}

/// `∂³/∂u∂v²` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_third_partial_uvv(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 1, 2, u, v))
}

/// `∂³/∂v³` of all degree-6 triangular basis functions.
pub fn bernstein_tri6_third_partial_vvv(u: f64, v: f64) -> [f64; TRI6_CP_COUNT] {
    to_28(derivative_basis(6, 0, 3, u, v))
}

fn derivative_basis(n: usize, du: usize, dv: usize, u: f64, v: f64) -> Vec<f64> {
    if du == 0 && dv == 0 {
        return eval_basis(n, u, v);
    }
    assert!(du + dv <= n);
    if du > 0 {
        let lower = derivative_basis(n - 1, du - 1, dv, u, v);
        lift_partial_u(n, &lower)
    } else {
        let lower = derivative_basis(n - 1, du, dv - 1, u, v);
        lift_partial_v(n, &lower)
    }
}

fn lift_partial_u(n: usize, lower: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0_f64; tri_cp_count(n)];
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            let term_i = if i >= 1 {
                lower[flat_index(n - 1, i - 1, j)]
            } else {
                0.0
            };
            let term_k = if k >= 1 {
                lower[flat_index(n - 1, i, j)]
            } else {
                0.0
            };
            out[flat_index(n, i, j)] = n as f64 * (term_i - term_k);
        }
    }
    out
}

fn lift_partial_v(n: usize, lower: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0_f64; tri_cp_count(n)];
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            let term_j = if j >= 1 {
                lower[flat_index(n - 1, i, j - 1)]
            } else {
                0.0
            };
            let term_k = if k >= 1 {
                lower[flat_index(n - 1, i, j)]
            } else {
                0.0
            };
            out[flat_index(n, i, j)] = n as f64 * (term_j - term_k);
        }
    }
    out
}

fn eval_basis(n: usize, u: f64, v: f64) -> Vec<f64> {
    let w = 1.0 - u - v;
    let mut out = vec![0.0_f64; tri_cp_count(n)];
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            out[flat_index(n, i, j)] =
                multinomial(n, i, j, k) * u.powi(i as i32) * v.powi(j as i32) * w.powi(k as i32);
        }
    }
    out
}

#[inline]
fn multinomial(n: usize, i: usize, j: usize, k: usize) -> f64 {
    debug_assert_eq!(i + j + k, n);
    factorial(n) / (factorial(i) * factorial(j) * factorial(k))
}

#[inline]
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, v| acc * v as f64)
}

fn to_28(v: Vec<f64>) -> [f64; TRI6_CP_COUNT] {
    assert_eq!(v.len(), TRI6_CP_COUNT);
    let mut out = [0.0_f64; TRI6_CP_COUNT];
    out.copy_from_slice(&v);
    out
}

#[cfg(test)]
#[path = "bernstein_tri6_tests.rs"]
mod tests;
