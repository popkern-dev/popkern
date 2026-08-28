//! Triangular total-degree-5 Bernstein-Bézier basis +
//! partials.
//!
//! Used by the K-P 2013 biquintic polar construction to represent
//! the per-sector **auxiliary triangular guide `g△`** over the
//! barycentric simplex
//!
//! ```text
//!   T = { (u, v, w) ∈ ℝ³ : u + v + w = 1, u, v, w ≥ 0 }.
//! ```
//!
//! In the Peters convention we parameterise the simplex by `(u, v)`
//! with `w := 1 − u − v`. The K-P 2013 construction uses `g△` with
//! total degree 5, so the basis has
//! `(5 + 2) choose 2 = 21` basis functions
//!
//! ```text
//!   B^5_{i, j, k}(u, v, w) = (5! / (i! j! k!)) · u^i · v^j · w^k
//!                             for i + j + k = 5.
//! ```
//!
//! # Directional partials
//!
//! Treating `(u, v)` as independent (with `w = 1 − u − v`):
//!
//! ```text
//!   ∂B^n_{i,j,k}/∂u = n · (B^{n-1}_{i-1,j,k}(u,v,w) − B^{n-1}_{i,j,k-1}(u,v,w))
//!   ∂B^n_{i,j,k}/∂v = n · (B^{n-1}_{i,j-1,k}(u,v,w) − B^{n-1}_{i,j,k-1}(u,v,w))
//! ```
//!
//! (out-of-range indices treated as 0). Second partials are
//! analogous and use the degree-(n − 2) basis.
//!
//! # Flat indexing
//!
//! Triples `(i, j, k)` with `i + j + k = n` are flattened into an
//! array of length `T(n) := (n + 1)(n + 2) / 2` in **row-major
//! order** with `i` as the outer index and `j` inner:
//!
//! ```text
//!   flat_index^n(i, j) = i · (2n + 3 − i) / 2 + j
//! ```
//!
//! For `n = 5`, `T(5) = 21`; for `n = 4`, `T(4) = 15`; for `n = 3`,
//! `T(3) = 10`. These are all that the K-P 2013 g△ construction
//! consumes (point + 1st + 2nd partial).

/// Number of triangular-Bezier control points at total degree `n`.
#[inline]
pub const fn tri_cp_count(n: usize) -> usize {
    (n + 1) * (n + 2) / 2
}

/// Number of CPs at degree 5: 21.
pub const TRI5_CP_COUNT: usize = 21;

/// Number of CPs at degree 4: 15.
pub const TRI4_CP_COUNT: usize = 15;

/// Number of CPs at degree 3: 10.
pub const TRI3_CP_COUNT: usize = 10;

/// Multinomial coefficient `n! / (i! j! k!)` in the row-major
/// ordering described in the module docs. For `n = 5`, 21 entries.
const M5: [f64; TRI5_CP_COUNT] = [
    // i = 0: (0, 0, 5), (0, 1, 4), …, (0, 5, 0)
    1.0, 5.0, 10.0, 10.0, 5.0, 1.0, // i = 1: (1, 0, 4), (1, 1, 3), …, (1, 4, 0)
    5.0, 20.0, 30.0, 20.0, 5.0, // i = 2: (2, 0, 3), (2, 1, 2), (2, 2, 1), (2, 3, 0)
    10.0, 30.0, 30.0, 10.0, // i = 3: (3, 0, 2), (3, 1, 1), (3, 2, 0)
    10.0, 20.0, 10.0, // i = 4: (4, 0, 1), (4, 1, 0)
    5.0, 5.0, // i = 5: (5, 0, 0)
    1.0,
];

/// Multinomial table at degree 4 (15 entries).
const M4: [f64; TRI4_CP_COUNT] = [
    // i = 0: (0, 0, 4), (0, 1, 3), (0, 2, 2), (0, 3, 1), (0, 4, 0)
    1.0, 4.0, 6.0, 4.0, 1.0, // i = 1: (1, 0, 3), (1, 1, 2), (1, 2, 1), (1, 3, 0)
    4.0, 12.0, 12.0, 4.0, // i = 2: (2, 0, 2), (2, 1, 1), (2, 2, 0)
    6.0, 12.0, 6.0, // i = 3: (3, 0, 1), (3, 1, 0)
    4.0, 4.0, // i = 4: (4, 0, 0)
    1.0,
];

/// Multinomial table at degree 3 (10 entries).
const M3: [f64; TRI3_CP_COUNT] = [
    // i = 0: (0, 0, 3), (0, 1, 2), (0, 2, 1), (0, 3, 0)
    1.0, 3.0, 3.0, 1.0, // i = 1: (1, 0, 2), (1, 1, 1), (1, 2, 0)
    3.0, 6.0, 3.0, // i = 2: (2, 0, 1), (2, 1, 0)
    3.0, 3.0, // i = 3: (3, 0, 0)
    1.0,
];

/// Flat index for `(i, j)` at degree `n` (triple `(i, j, n - i - j)`).
/// Panics if `i + j > n`.
#[inline]
pub const fn flat_index(n: usize, i: usize, j: usize) -> usize {
    assert!(i + j <= n);
    i * (2 * n + 3 - i) / 2 + j
}

/// Small slack above `1.0` for the `u + v ≤ 1` simplex guard —
/// accommodates FD-difference samples that step slightly past the
/// edge (`u + h + v` with h ≈ 1e-6).
const SIMPLEX_EPSILON: f64 = 1e-6;

/// Asserts `(u, v)` lies in (or within `SIMPLEX_EPSILON` of) the
/// barycentric simplex `u ≥ 0, v ≥ 0, u + v ≤ 1`. Debug-only;
/// `w = 1 − u − v` silently goes negative outside the simplex
/// and `w^k` produces meaningless values that callers will
/// rarely want.
#[inline]
fn debug_assert_simplex(u: f64, v: f64) {
    debug_assert!(
        u >= -SIMPLEX_EPSILON && v >= -SIMPLEX_EPSILON && u + v <= 1.0 + SIMPLEX_EPSILON,
        "triangular basis eval requires u ≥ 0, v ≥ 0, u + v ≤ 1; got (u, v) = ({u}, {v})"
    );
}

/// Evaluate all 21 degree-5 triangular basis functions at
/// `(u, v, w = 1 − u − v)`. Output ordering matches [`flat_index`].
///
/// # Panics (debug)
/// Panics in debug builds if `(u, v)` lies outside the barycentric
/// simplex (with a small epsilon for FD samples). Release builds
/// silently compute `w^k` with negative `w`, which is meaningless
/// for most users — see [`bernstein_tri5_all_extended`] for the
/// variant that skips the check.
pub fn bernstein_tri5_all(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    debug_assert_simplex(u, v);
    let w = 1.0 - u - v;
    eval_monomial_scaled::<TRI5_CP_COUNT>(5, u, v, w, &M5)
}

/// Same as [`bernstein_tri5_all`] but **without** the simplex
/// check — suitable for evaluating the tri-5 polynomial at points
/// *outside* the natural wedge (KP13 §3 "v_k domain extension").
/// The polynomial identity `B^5_{i,j,k} = (5!/(i!j!k!)) u^i v^j w^k`
/// is well-defined for any real `(u, v, w)`.
pub fn bernstein_tri5_all_extended(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    let w = 1.0 - u - v;
    eval_monomial_scaled::<TRI5_CP_COUNT>(5, u, v, w, &M5)
}

/// Same as [`bernstein_tri4_all`] but without the simplex check.
pub fn bernstein_tri4_all_extended(u: f64, v: f64) -> [f64; TRI4_CP_COUNT] {
    let w = 1.0 - u - v;
    eval_monomial_scaled::<TRI4_CP_COUNT>(4, u, v, w, &M4)
}

/// Same as [`bernstein_tri3_all`] but without the simplex check.
pub fn bernstein_tri3_all_extended(u: f64, v: f64) -> [f64; TRI3_CP_COUNT] {
    let w = 1.0 - u - v;
    eval_monomial_scaled::<TRI3_CP_COUNT>(3, u, v, w, &M3)
}

/// Evaluate all 15 degree-4 basis functions. See
/// [`bernstein_tri5_all`] for simplex-domain caveat.
pub fn bernstein_tri4_all(u: f64, v: f64) -> [f64; TRI4_CP_COUNT] {
    debug_assert_simplex(u, v);
    let w = 1.0 - u - v;
    eval_monomial_scaled::<TRI4_CP_COUNT>(4, u, v, w, &M4)
}

/// Evaluate all 10 degree-3 basis functions. See
/// [`bernstein_tri5_all`] for simplex-domain caveat.
pub fn bernstein_tri3_all(u: f64, v: f64) -> [f64; TRI3_CP_COUNT] {
    debug_assert_simplex(u, v);
    let w = 1.0 - u - v;
    eval_monomial_scaled::<TRI3_CP_COUNT>(3, u, v, w, &M3)
}

/// Internal: ∂/∂u lift from a degree-4 basis sample.
fn partial_u_from_b4(b4: [f64; TRI4_CP_COUNT]) -> [f64; TRI5_CP_COUNT] {
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let term_i = if i >= 1 {
                b4[flat_index(4, i - 1, j)]
            } else {
                0.0
            };
            let term_k = if k >= 1 { b4[flat_index(4, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 5.0 * (term_i - term_k);
        }
    }
    out
}

/// Internal: ∂/∂v lift from a degree-4 basis sample.
fn partial_v_from_b4(b4: [f64; TRI4_CP_COUNT]) -> [f64; TRI5_CP_COUNT] {
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let term_j = if j >= 1 {
                b4[flat_index(4, i, j - 1)]
            } else {
                0.0
            };
            let term_k = if k >= 1 { b4[flat_index(4, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 5.0 * (term_j - term_k);
        }
    }
    out
}

/// `∂/∂u` extended variant — no simplex check.
pub fn bernstein_tri5_partial_u_extended(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    partial_u_from_b4(bernstein_tri4_all_extended(u, v))
}

/// `∂/∂v` extended variant — no simplex check.
pub fn bernstein_tri5_partial_v_extended(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    partial_v_from_b4(bernstein_tri4_all_extended(u, v))
}

/// `∂/∂u` of all 21 degree-5 basis functions, per the formula
/// `∂B^5_{i,j,k}/∂u = 5 · (B^4_{i-1,j,k} − B^4_{i,j,k-1})`. Output
/// ordering matches [`flat_index`] at degree 5.
pub fn bernstein_tri5_partial_u(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    let b4 = bernstein_tri4_all(u, v);
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let term_i = if i >= 1 {
                b4[flat_index(4, i - 1, j)]
            } else {
                0.0
            };
            let term_k = if k >= 1 { b4[flat_index(4, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 5.0 * (term_i - term_k);
        }
    }
    out
}

/// `∂/∂v` of all 21 degree-5 basis functions.
pub fn bernstein_tri5_partial_v(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    let b4 = bernstein_tri4_all(u, v);
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let term_j = if j >= 1 {
                b4[flat_index(4, i, j - 1)]
            } else {
                0.0
            };
            let term_k = if k >= 1 { b4[flat_index(4, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 5.0 * (term_j - term_k);
        }
    }
    out
}

/// Internal: ∂²/∂u² lift from a degree-3 basis sample.
fn second_partial_uu_from_b3(b3: [f64; TRI3_CP_COUNT]) -> [f64; TRI5_CP_COUNT] {
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let t_ii = if i >= 2 {
                b3[flat_index(3, i - 2, j)]
            } else {
                0.0
            };
            let t_ik = if i >= 1 && k >= 1 {
                b3[flat_index(3, i - 1, j)]
            } else {
                0.0
            };
            let t_kk = if k >= 2 { b3[flat_index(3, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 20.0 * (t_ii - 2.0 * t_ik + t_kk);
        }
    }
    out
}

/// Internal: ∂²/∂v² lift from a degree-3 basis sample.
fn second_partial_vv_from_b3(b3: [f64; TRI3_CP_COUNT]) -> [f64; TRI5_CP_COUNT] {
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let t_jj = if j >= 2 {
                b3[flat_index(3, i, j - 2)]
            } else {
                0.0
            };
            let t_jk = if j >= 1 && k >= 1 {
                b3[flat_index(3, i, j - 1)]
            } else {
                0.0
            };
            let t_kk = if k >= 2 { b3[flat_index(3, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 20.0 * (t_jj - 2.0 * t_jk + t_kk);
        }
    }
    out
}

/// Internal: ∂²/∂u∂v lift from a degree-3 basis sample.
fn second_partial_uv_from_b3(b3: [f64; TRI3_CP_COUNT]) -> [f64; TRI5_CP_COUNT] {
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let t_ij = if i >= 1 && j >= 1 {
                b3[flat_index(3, i - 1, j - 1)]
            } else {
                0.0
            };
            let t_ik = if i >= 1 && k >= 1 {
                b3[flat_index(3, i - 1, j)]
            } else {
                0.0
            };
            let t_jk = if j >= 1 && k >= 1 {
                b3[flat_index(3, i, j - 1)]
            } else {
                0.0
            };
            let t_kk = if k >= 2 { b3[flat_index(3, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 20.0 * (t_ij - t_ik - t_jk + t_kk);
        }
    }
    out
}

/// `∂²/∂u²` extended variant — no simplex check.
pub fn bernstein_tri5_second_partial_uu_extended(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    second_partial_uu_from_b3(bernstein_tri3_all_extended(u, v))
}

/// `∂²/∂v²` extended variant — no simplex check.
pub fn bernstein_tri5_second_partial_vv_extended(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    second_partial_vv_from_b3(bernstein_tri3_all_extended(u, v))
}

/// `∂²/∂u∂v` extended variant — no simplex check.
pub fn bernstein_tri5_second_partial_uv_extended(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    second_partial_uv_from_b3(bernstein_tri3_all_extended(u, v))
}

/// `∂²/∂u²` of all 21 degree-5 basis functions, via
/// `∂²B^5/∂u² = 20 · (B^3_{i-2,j,k} − 2·B^3_{i-1,j,k-1} + B^3_{i,j,k-2})`.
pub fn bernstein_tri5_second_partial_uu(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    let b3 = bernstein_tri3_all(u, v);
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let t_ii = if i >= 2 {
                b3[flat_index(3, i - 2, j)]
            } else {
                0.0
            };
            let t_ik = if i >= 1 && k >= 1 {
                b3[flat_index(3, i - 1, j)]
            } else {
                0.0
            };
            let t_kk = if k >= 2 { b3[flat_index(3, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 20.0 * (t_ii - 2.0 * t_ik + t_kk);
        }
    }
    out
}

/// `∂²/∂v²` of all 21 degree-5 basis functions.
pub fn bernstein_tri5_second_partial_vv(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    let b3 = bernstein_tri3_all(u, v);
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let t_jj = if j >= 2 {
                b3[flat_index(3, i, j - 2)]
            } else {
                0.0
            };
            let t_jk = if j >= 1 && k >= 1 {
                b3[flat_index(3, i, j - 1)]
            } else {
                0.0
            };
            let t_kk = if k >= 2 { b3[flat_index(3, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 20.0 * (t_jj - 2.0 * t_jk + t_kk);
        }
    }
    out
}

/// `∂²/∂u∂v` of all 21 degree-5 basis functions, via
/// `∂²B^5/∂u∂v = 20 · (B^3_{i-1,j-1,k} − B^3_{i-1,j,k-1} − B^3_{i,j-1,k-1} + B^3_{i,j,k-2})`.
pub fn bernstein_tri5_second_partial_uv(u: f64, v: f64) -> [f64; TRI5_CP_COUNT] {
    let b3 = bernstein_tri3_all(u, v);
    let mut out = [0.0_f64; TRI5_CP_COUNT];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let k = 5 - i - j;
            let t_ij = if i >= 1 && j >= 1 {
                b3[flat_index(3, i - 1, j - 1)]
            } else {
                0.0
            };
            let t_ik = if i >= 1 && k >= 1 {
                b3[flat_index(3, i - 1, j)]
            } else {
                0.0
            };
            let t_jk = if j >= 1 && k >= 1 {
                b3[flat_index(3, i, j - 1)]
            } else {
                0.0
            };
            let t_kk = if k >= 2 { b3[flat_index(3, i, j)] } else { 0.0 };
            out[flat_index(5, i, j)] = 20.0 * (t_ij - t_ik - t_jk + t_kk);
        }
    }
    out
}

/// Private helper — builds the `n!/(i!j!k!) · u^i · v^j · w^k` table
/// from a precomputed multinomial table `m` of the right length.
fn eval_monomial_scaled<const LEN: usize>(
    n: usize,
    u: f64,
    v: f64,
    w: f64,
    m: &[f64; LEN],
) -> [f64; LEN] {
    // Pre-compute u^i, v^j, w^k tables.
    let max = n + 1;
    let mut u_pow = [1.0_f64; 16];
    let mut v_pow = [1.0_f64; 16];
    let mut w_pow = [1.0_f64; 16];
    debug_assert!(max <= 16);
    for p in 1..max {
        u_pow[p] = u_pow[p - 1] * u;
        v_pow[p] = v_pow[p - 1] * v;
        w_pow[p] = w_pow[p - 1] * w;
    }
    let mut out = [0.0_f64; LEN];
    let mut idx = 0;
    for (i, &ui) in u_pow.iter().enumerate().take(max) {
        for (j, &vj) in v_pow.iter().enumerate().take(max - i) {
            let k = n - i - j;
            out[idx] = m[idx] * ui * vj * w_pow[k];
            idx += 1;
        }
    }
    out
}

#[cfg(test)]
#[path = "bernstein_tri5_tests.rs"]
mod tests;
