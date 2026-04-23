//! Phase IV PnS — Tensor-border input + `b ∘ β` boundary evaluator.
//!
//! The PnS n-gon cap consumes a ring of **bi-2 tensor-border**
//! strips — one per sector — representing the position and first
//! transversal derivative of the surrounding bi-2 spline complex
//! along the cap's outer boundary. This module defines the ring
//! type and evaluates `b ∘ β` to produce the cap's **outer two
//! layers** (`j = 0, 1` rows in the cap's bi-4 or bi-3 BB form).
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` §3 + Fig. 6b for the tensor-border
//! structure; §4 eq (3) for β; §4.2 eq (3³) for the bi-3 β variant.
//!
//! # Per-sector strip
//!
//! Each sector `k ∈ 0..n` has a bi-2 strip consisting of 3 × 2 = 6
//! Point3 BB-coefficients:
//! ```text
//!   b^k_{i,0}   for i ∈ {0, 1, 2}   — position row  (b(u, 0))
//!   b^k_{i,1}   for i ∈ {0, 1, 2}   — first-derivative row (∂_v b(u, 0))
//! ```
//!
//! The "position row" gives `p(u, 0) = Σ_i b^k_{i,0} B_i^2(u)` and
//! the "derivative row" gives `∂_v p(u, 0) = Σ_i b^k_{i,1} B_i^2(u)`.
//! Both are degree-2 polynomials in `u`.
//!
//! # `b ∘ β` via Taylor expansion
//!
//! With `β(u, v) = (u + b(u)·v, a(u)·v)` (eq 3), the cap's outer
//! layers correspond to:
//! ```text
//!   (b ∘ β)(u, 0)     = b(u, 0) = position row                  (degree 2 in u)
//!   ∂_v (b ∘ β)(u, 0) = b(u) · ∂_u b(u, 0) + a(u) · ∂_v b(u, 0)
//! ```
//!
//! The first expression is just the original position row,
//! degree-raised to match the cap's output degree.
//!
//! The second is a polynomial sum of two products:
//! - `b(u) · ∂_u b(u, 0)`: deg 3 · deg 1 = deg 4
//! - `a(u) · ∂_v b(u, 0)`: deg 2 · deg 2 = deg 4
//!
//! Both products are expressed in degree-4 BB form via the standard
//! BB polynomial-multiplication formula
//! `c_k = Σ_{i+j=k} (C(m,i)·C(n,j)/C(m+n,k))·a_i·b_j`.
//!
//! For the **bi-3 n ∈ {3, 5} variant**, β has simpler closed forms
//! (a deg 1, b deg 2), so the j = 1 row comes out as degree 3,
//! fitting the bi-3 cap directly without further elevation.

#![allow(dead_code)]

use super::pns_beta::BetaReparam;
use crate::geometry::Point3;

/// A bi-2 tensor-border strip for one sector.
///
/// `pos[i]` is the position BB-coef at `u = i/2` parametric
/// position (degree-2 Bernstein). `deriv[i]` is the transversal
/// derivative at the same positions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bi2BorderStrip {
    /// Position row: `b^k_{i, 0}` for `i ∈ {0, 1, 2}`.
    pub pos: [Point3; 3],
    /// Derivative row: `b^k_{i, 1}` for `i ∈ {0, 1, 2}`.
    pub deriv: [Point3; 3],
}

/// Full ring of tensor-border strips, one per sector.
#[derive(Debug, Clone)]
pub struct Bi2TensorBorderRing {
    /// `sectors[k]` is sector `k`'s strip. Length `n`.
    pub sectors: Vec<Bi2BorderStrip>,
}

impl Bi2TensorBorderRing {
    /// Valence (number of sectors).
    pub fn n(&self) -> usize {
        self.sectors.len()
    }
}

/// Apply `b ∘ β` and degree-raise to produce the cap's outer two
/// layers (rows `j = 0, 1`) at **bi-4**, for one sector.
///
/// Returns `(row0, row1)` where each is a 5-Point3 array at degree
/// 4 in `u`. Matches the cap-patch's outer boundary-row layout
/// (`p^k_{i, 0}, p^k_{i, 1}` for `i ∈ 0..=4`).
pub fn compose_bi4_outer_layers(
    strip: &Bi2BorderStrip,
    beta: &BetaReparam,
) -> ([Point3; 5], [Point3; 5]) {
    // j = 0 row: position row degree-elevated from 2 → 4.
    let row0 = elevate_bi2_to_bi4(&strip.pos);

    // ∂_u b(u, 0): degree-1 BB (derivative of degree-2 in u).
    // For BB-form p(u) = Σ p_i B_i^2(u), p'(u) = 2·Σ (p_{i+1} − p_i) B_i^1(u).
    let dpos_du: [Point3; 2] = [
        scale_pt(2.0, sub_pt(strip.pos[1], strip.pos[0])),
        scale_pt(2.0, sub_pt(strip.pos[2], strip.pos[1])),
    ];

    // j = 1 row: b(u) · dpos_du + a(u) · deriv, both products are
    // degree-4 polynomials in u.
    //
    // b(u) is degree-3 with BB-coefs beta.b[0..=3].
    // a(u) is degree-2 with BB-coefs beta.a[0..=2].
    // dpos_du is degree-1 with BB-coefs (as Point3) above.
    // deriv is degree-2 with BB-coefs strip.deriv[0..=2].
    let prod_b_dpos = bb_mul_scalar_vec(&beta.b, &dpos_du); // deg 3 × deg 1 = deg 4 → 5 coefs
    let prod_a_deriv = bb_mul_scalar_vec(&beta.a, &strip.deriv); // deg 2 × deg 2 = deg 4 → 5 coefs

    let row1: [Point3; 5] = [
        add_pt(prod_b_dpos[0], prod_a_deriv[0]),
        add_pt(prod_b_dpos[1], prod_a_deriv[1]),
        add_pt(prod_b_dpos[2], prod_a_deriv[2]),
        add_pt(prod_b_dpos[3], prod_a_deriv[3]),
        add_pt(prod_b_dpos[4], prod_a_deriv[4]),
    ];

    (row0, row1)
}

/// Degree-raise a bi-2 BB polynomial (3 coefs) to bi-4 (5 coefs).
///
/// Standard Bézier degree-elevation from degree 2 to degree 4:
/// two successive degree-1 elevations, or the direct formula
/// `q_i = Σ_{j=max(0,i-2)}^{min(2,i)} (C(2,j) C(2,i-j) / C(4,i)) · p_j`.
fn elevate_bi2_to_bi4(p: &[Point3; 3]) -> [Point3; 5] {
    // Coefficients of `C(2,j) · C(2,i-j) / C(4,i)` in each output
    // position. Direct derivation:
    //   q_0 = p_0
    //   q_1 = (1/2) p_0 + (1/2) p_1
    //   q_2 = (1/6) p_0 + (4/6) p_1 + (1/6) p_2
    //   q_3 = (1/2) p_1 + (1/2) p_2
    //   q_4 = p_2
    [
        p[0],
        Point3::new(
            0.5 * p[0].x + 0.5 * p[1].x,
            0.5 * p[0].y + 0.5 * p[1].y,
            0.5 * p[0].z + 0.5 * p[1].z,
        ),
        Point3::new(
            (1.0 / 6.0) * p[0].x + (4.0 / 6.0) * p[1].x + (1.0 / 6.0) * p[2].x,
            (1.0 / 6.0) * p[0].y + (4.0 / 6.0) * p[1].y + (1.0 / 6.0) * p[2].y,
            (1.0 / 6.0) * p[0].z + (4.0 / 6.0) * p[1].z + (1.0 / 6.0) * p[2].z,
        ),
        Point3::new(
            0.5 * p[1].x + 0.5 * p[2].x,
            0.5 * p[1].y + 0.5 * p[2].y,
            0.5 * p[1].z + 0.5 * p[2].z,
        ),
        p[2],
    ]
}

/// BB product of a scalar polynomial (`a`, degree `M = a.len() − 1`)
/// and a Point3 polynomial (`p`, degree `N = p.len() − 1`) at total
/// degree `M + N`.
///
/// Formula: `(a · p)_k = Σ_{i+j=k} w_{i,j,k} · a_i · p_j` where
/// `w_{i,j,k} = C(M, i) · C(N, j) / C(M + N, k)`.
///
/// Supports `M + N ∈ {2, 3, 4}` (enough for bi-2/bi-3/bi-4 outputs).
fn bb_mul_scalar_vec(a: &[f64], p: &[Point3]) -> Vec<Point3> {
    let m = a.len() - 1;
    let n = p.len() - 1;
    let mn = m + n;
    let mut out = vec![Point3::new(0.0, 0.0, 0.0); mn + 1];
    for (k, out_k) in out.iter_mut().enumerate().take(mn + 1) {
        for (i, &a_i) in a.iter().enumerate().take(m.min(k) + 1) {
            let j = k - i;
            if j > n {
                continue;
            }
            let w = (binom(m, i) * binom(n, j)) / binom(mn, k);
            let coef = w * a_i;
            out_k.x += coef * p[j].x;
            out_k.y += coef * p[j].y;
            out_k.z += coef * p[j].z;
        }
    }
    out
}

fn binom(n: usize, k: usize) -> f64 {
    // Small binomials only (up to C(4, 2) = 6 in our use case).
    if k > n {
        return 0.0;
    }
    let mut num = 1.0_f64;
    let mut den = 1.0_f64;
    for i in 0..k {
        num *= (n - i) as f64;
        den *= (i + 1) as f64;
    }
    num / den
}

#[inline]
fn add_pt(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

#[inline]
fn sub_pt(a: Point3, b: Point3) -> Point3 {
    Point3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

#[inline]
fn scale_pt(s: f64, p: Point3) -> Point3 {
    Point3::new(s * p.x, s * p.y, s * p.z)
}

#[cfg(test)]
#[path = "pns_boundary_tests.rs"]
mod tests;
