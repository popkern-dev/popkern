//! Phase IV PnS — F_κ quality functional evaluator for bi-4 patches.
//!
//! `F_κ(p) := ∫∫ Σ_{i+j=κ, i,j ≥ 0} (κ! / (i!·j!)) · (∂^i_s ∂^j_t p)² ds dt`
//!
//! over `(s, t) ∈ [0, 1]²`. The functional penalises `κ`-th-order
//! variation of the patch — higher `κ` suppresses higher-order
//! oscillation, producing smoother fills. For the PnS bi-4 n-gon
//! cap (`docs/papers/14bi234.pdf` §4), the paper selects `F_3` for
//! `n ∈ {5, 6, 7}` and `F_4` for `n ∈ {8, 9, 10}` as the functional
//! minimised over the `n + 5` remaining free
//! BB-coefs after all G¹ / b∘β / central-quadratic constraints have
//! pinned the rest.
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` §4 "Construction: Multi-sided caps".
//! The `F_κ` functional is defined in a general form in §2 of the
//! same paper and in fairness-functional surveys referenced there
//! ([24]–[29] of the paper's bibliography).
//!
//! # Implementation
//!
//! We evaluate `F_κ` on a bi-4 patch by Gauss-Legendre quadrature
//! in 2D. A bi-4 patch has degree 4 in each of `s, t`. Its `κ`-th
//! partial derivatives have degree `4 − i` in `s` and `4 − j` in
//! `t`. The squared integrand therefore has degree at most `8` in
//! each variable (when `i = 0` or `j = 0`).
//!
//! Gauss-Legendre with 5 nodes per dimension is exact for
//! polynomials up to degree 9 and therefore exact for the full
//! `F_κ` integrand (κ ≥ 1) on a bi-4 patch. We use the 5-point
//! rule.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use crate::geometry::Point3;

/// Evaluate `F_κ(p)` on a bi-4 patch `p` with `5 × 5` BB-coefs.
/// Returns a scalar energy (non-negative); zero iff all κ-th
/// partial derivatives vanish identically.
///
/// `kappa` must be in `1..=4` (the functional is well-defined for
/// any `κ ≥ 0`, but `F_0` is just the squared L² norm which the
/// paper doesn't use for cap fairness).
pub fn f_kappa_bi4(grid: &[[Point3; 5]; 5], kappa: usize) -> f64 {
    assert!(
        (1..=4).contains(&kappa),
        "f_kappa_bi4: kappa must be in 1..=4, got {kappa}"
    );
    // 5-point Gauss-Legendre nodes + weights on [0, 1].
    let (nodes, weights) = gauss_legendre_5_unit_interval();

    // Accumulate the double integral.
    let mut total = 0.0_f64;
    for (is, &s) in nodes.iter().enumerate() {
        for (it, &t) in nodes.iter().enumerate() {
            // At (s, t), sum over all (i, j) with i+j=κ:
            //   (κ!/(i!·j!)) · (∂^i_s ∂^j_t p)²
            let mut pointwise = 0.0_f64;
            for i in 0..=kappa {
                let j = kappa - i;
                let coef = factorial(kappa) as f64 / (factorial(i) as f64 * factorial(j) as f64);
                let d = partial_bi4(grid, i, j, s, t);
                let d_sq = d.x * d.x + d.y * d.y + d.z * d.z;
                pointwise += coef * d_sq;
            }
            total += weights[is] * weights[it] * pointwise;
        }
    }
    total
}

/// Evaluate the mixed partial `∂^i_s ∂^j_t p` of a bi-4 patch at
/// `(s, t)`. Works for any `i, j ≥ 0` with `i ≤ 4, j ≤ 4`.
fn partial_bi4(grid: &[[Point3; 5]; 5], i: usize, j: usize, s: f64, t: f64) -> Point3 {
    let bs = bernstein_bi4_deriv_order(i, s);
    let bt = bernstein_bi4_deriv_order(j, t);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for a in 0..5 {
        for b in 0..5 {
            let w = bs[a] * bt[b];
            out.x += w * grid[a][b].x;
            out.y += w * grid[a][b].y;
            out.z += w * grid[a][b].z;
        }
    }
    out
}

/// All 5 basis values `[B_0^4^{(k)}(u), .., B_4^4^{(k)}(u)]` for
/// the `k`-th derivative (`k ∈ 0..=4`).
fn bernstein_bi4_deriv_order(k: usize, u: f64) -> [f64; 5] {
    match k {
        0 => bi4_all(u),
        1 => bi4_deriv_all(u),
        2 => bi4_second_deriv_all(u),
        3 => bi4_third_deriv_all(u),
        4 => bi4_fourth_deriv_all(u),
        _ => [0.0; 5],
    }
}

/// `B_i^4(u)` for `i ∈ 0..=4`.
fn bi4_all(u: f64) -> [f64; 5] {
    let omu = 1.0 - u;
    let c: [f64; 5] = [1.0, 4.0, 6.0, 4.0, 1.0];
    let mut out = [0.0_f64; 5];
    let mut tp = 1.0_f64;
    let mut op_p = [0.0_f64; 5];
    op_p[0] = 1.0;
    for i in 1..5 {
        op_p[i] = op_p[i - 1] * omu;
    }
    for i in 0..5 {
        // op_p holds [1, omu, omu^2, omu^3, omu^4] indexed by exponent
        let om_pow = powi_f64(omu, (4 - i) as i32);
        out[i] = c[i] * tp * om_pow;
        tp *= u;
    }
    out
}

/// `B_i^4'(u) = 4·(B_{i-1}^3(u) − B_i^3(u))`.
fn bi4_deriv_all(u: f64) -> [f64; 5] {
    let b3 = bi3_all(u);
    [
        4.0 * (0.0 - b3[0]),
        4.0 * (b3[0] - b3[1]),
        4.0 * (b3[1] - b3[2]),
        4.0 * (b3[2] - b3[3]),
        4.0 * (b3[3] - 0.0),
    ]
}

/// `B_i^4''(u) = 12·(B_{i-2}^2(u) − 2·B_{i-1}^2(u) + B_i^2(u))`.
fn bi4_second_deriv_all(u: f64) -> [f64; 5] {
    let b2 = bi2_all(u);
    let bm = |k: i32| -> f64 {
        if (0..=2).contains(&k) {
            b2[k as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 5];
    for i in 0..5 {
        let ii = i as i32;
        out[i] = 12.0 * (bm(ii - 2) - 2.0 * bm(ii - 1) + bm(ii));
    }
    out
}

/// `B_i^4'''(u) = 24·(B_{i-3}^1 − 3·B_{i-2}^1 + 3·B_{i-1}^1 − B_i^1)`.
fn bi4_third_deriv_all(u: f64) -> [f64; 5] {
    let b1 = [1.0 - u, u];
    let bm = |k: i32| -> f64 {
        if (0..=1).contains(&k) {
            b1[k as usize]
        } else {
            0.0
        }
    };
    let mut out = [0.0_f64; 5];
    for i in 0..5 {
        let ii = i as i32;
        out[i] = 24.0 * (bm(ii - 3) - 3.0 * bm(ii - 2) + 3.0 * bm(ii - 1) - bm(ii));
    }
    out
}

/// `B_i^4''''(u) = 24·(Δ^4 B_i^0)` — piecewise-constant
/// `[24, −96, 144, −96, 24]` on `(0, 1)` (the "4th difference" of
/// the delta function basis at degree 0).
fn bi4_fourth_deriv_all(_u: f64) -> [f64; 5] {
    [24.0, -96.0, 144.0, -96.0, 24.0]
}

/// `B_i^3(u)` for `i ∈ 0..=3`.
fn bi3_all(u: f64) -> [f64; 4] {
    let omu = 1.0 - u;
    [
        omu * omu * omu,
        3.0 * u * omu * omu,
        3.0 * u * u * omu,
        u * u * u,
    ]
}

/// `B_i^2(u)` for `i ∈ 0..=2`.
fn bi2_all(u: f64) -> [f64; 3] {
    let omu = 1.0 - u;
    [omu * omu, 2.0 * u * omu, u * u]
}

/// Gauss-Legendre 5-point rule on `[0, 1]` — exact for polynomials
/// up to degree 9.
fn gauss_legendre_5_unit_interval() -> ([f64; 5], [f64; 5]) {
    // Standard Gauss-Legendre on [-1, 1] with 5 nodes, then
    // transform to [0, 1] via x_unit = (x + 1) / 2,
    // w_unit = w / 2.
    // Nodes on [-1, 1]:
    //   ±1/3·√(5 − 2√(10/7)),  ±1/3·√(5 + 2√(10/7)),  0
    // Weights:
    //   (322 + 13·√70)/900,  (322 − 13·√70)/900,  128/225
    let a1 = (5.0 - 2.0 * (10.0_f64 / 7.0).sqrt()).sqrt() / 3.0;
    let a2 = (5.0 + 2.0 * (10.0_f64 / 7.0).sqrt()).sqrt() / 3.0;
    let w0 = 128.0 / 225.0;
    let w1 = (322.0 + 13.0 * 70.0_f64.sqrt()) / 900.0;
    let w2 = (322.0 - 13.0 * 70.0_f64.sqrt()) / 900.0;
    let nodes = [
        (-a2 + 1.0) * 0.5,
        (-a1 + 1.0) * 0.5,
        0.5,
        (a1 + 1.0) * 0.5,
        (a2 + 1.0) * 0.5,
    ];
    let weights = [w2 * 0.5, w1 * 0.5, w0 * 0.5, w1 * 0.5, w2 * 0.5];
    (nodes, weights)
}

fn powi_f64(x: f64, n: i32) -> f64 {
    if n == 0 { 1.0 } else { x.powi(n) }
}

/// Small-integer factorial (up to `κ = 4` is the max needed here).
fn factorial(n: usize) -> usize {
    (1..=n).product::<usize>().max(1)
}

#[cfg(test)]
#[path = "pns_f_kappa_tests.rs"]
mod tests;
