//! KP20 cubic τ reparameterisation and third-order guide composition.
//!
//! KP20 Appendix A publishes τ as a total-degree-3 triangular
//! Bernstein map into the sector wedge. This module builds that map
//! from equations (22) and (23), evaluates its derivatives through
//! order three, and composes a KP20 degree-6 guide jet with τ.

#![allow(clippy::needless_range_loop)]

use super::kp20_bi6_guide_sample::{
    Kp20GuideJet3, Kp20Tri6GuideSector, kp20_corner_block_from_jet,
};
use super::kp20_bi6_math::Kp20TauParams;
use crate::geometry::Point3;

/// Number of total-degree-3 triangular BB coefficients in a τ map.
pub const KP20_TAU_CP_COUNT: usize = 10;

/// Total-degree-3 KP20 τ map in triangular BB form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20TauMap {
    /// Two-dimensional BB coefficients in `flat_index(3, i, j)` order.
    pub cps: [[f64; 2]; KP20_TAU_CP_COUNT],
}

impl Kp20TauMap {
    /// Construct from explicit degree-3 triangular BB coefficients.
    pub const fn from_cps(cps: [[f64; 2]; KP20_TAU_CP_COUNT]) -> Self {
        Self { cps }
    }

    /// Build KP20 equation (22)'s sector map for valence `n` and sector `k`.
    pub fn for_sector(n: usize, k: usize) -> Result<Self, String> {
        let params = Kp20TauParams::for_valence(n)?;
        let alpha = core::f64::consts::TAU / n as f64;
        let a0 = k as f64 * alpha;
        let a1 = (k + 1) as f64 * alpha;
        let v_dot = [a0.cos(), a0.sin()];
        let v_ddot = [a1.cos(), a1.sin()];
        let mut cps = [[0.0_f64; 2]; KP20_TAU_CP_COUNT];

        cps[flat_index(3, 0, 0)] = [0.0, 0.0];
        cps[flat_index(3, 0, 1)] = scale2(v_dot, 1.0 / 3.0);
        cps[flat_index(3, 0, 2)] = scale2(v_dot, 2.0 / 3.0);
        cps[flat_index(3, 1, 0)] = scale2(v_ddot, 1.0 / 3.0);
        cps[flat_index(3, 2, 0)] = scale2(v_ddot, 2.0 / 3.0);
        cps[flat_index(3, 0, 3)] = scale2(v_dot, params.mu0);
        cps[flat_index(3, 3, 0)] = scale2(v_ddot, params.mu0);
        cps[flat_index(3, 1, 1)] = add2(scale2(v_dot, params.mu1), scale2(v_ddot, params.mu1));
        cps[flat_index(3, 1, 2)] = add2(scale2(v_dot, params.mu2), scale2(v_ddot, params.mu3));
        cps[flat_index(3, 2, 1)] = add2(scale2(v_dot, params.mu3), scale2(v_ddot, params.mu2));

        Ok(Self { cps })
    }

    /// Build a degree-3 triangular representation of an affine map
    /// `origin + u * du + v * dv`. Used as a compatibility bridge
    /// and as a test oracle for the cubic chain rule.
    pub fn affine(origin: [f64; 2], du: [f64; 2], dv: [f64; 2]) -> Self {
        let mut cps = [[0.0_f64; 2]; KP20_TAU_CP_COUNT];
        for i in 0..=3 {
            for j in 0..=(3 - i) {
                cps[flat_index(3, i, j)] = add2(
                    add2(origin, scale2(du, i as f64 / 3.0)),
                    scale2(dv, j as f64 / 3.0),
                );
            }
        }
        Self { cps }
    }

    /// Evaluate value and partial derivatives through order three at `(u, v)`.
    pub fn jet3(&self, u: f64, v: f64) -> Kp20TauJet3 {
        let mut out = [[[0.0_f64; 2]; 4]; 4];
        for a in 0..4 {
            for b in 0..(4 - a) {
                let weights = derivative_basis(3, a, b, u, v);
                for idx in 0..KP20_TAU_CP_COUNT {
                    out[a][b][0] += weights[idx] * self.cps[idx][0];
                    out[a][b][1] += weights[idx] * self.cps[idx][1];
                }
            }
        }
        Kp20TauJet3 { d: out }
    }
}

/// τ derivatives `d[p][q] = ∂u^p ∂v^q τ` for `p + q <= 3`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20TauJet3 {
    /// Two-dimensional derivative values.
    pub d: [[[f64; 2]; 4]; 4],
}

/// Compose `guide ∘ tau` and return its cubic Hermite jet.
pub fn kp20_compose_guide_tau_jet(
    guide: &Kp20Tri6GuideSector,
    tau: &Kp20TauMap,
    u: f64,
    v: f64,
) -> Kp20GuideJet3 {
    let tau_jet = tau.jet3(u, v);
    let guide_xy = tau_jet.d[0][0];
    let guide_jet = guide.jet3(guide_xy[0], guide_xy[1]);
    compose_guide_with_tau_jets(guide_jet, tau_jet)
}

/// Sample one bi-6 corner block from `guide ∘ tau` at patch-local `(u, v)`.
pub fn kp20_sample_corner_block_with_tau(
    guide: &Kp20Tri6GuideSector,
    tau: &Kp20TauMap,
    u: f64,
    v: f64,
) -> [[Point3; 4]; 4] {
    kp20_corner_block_from_jet(kp20_compose_guide_tau_jet(guide, tau, u, v))
}

fn compose_guide_with_tau_jets(g: Kp20GuideJet3, tau: Kp20TauJet3) -> Kp20GuideJet3 {
    let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for a in 0..4 {
        for b in 0..(4 - a) {
            out[a][b] = chain_rule_entry(g, tau, a, b);
        }
    }
    Kp20GuideJet3 { h: out }
}

fn chain_rule_entry(g: Kp20GuideJet3, tau: Kp20TauJet3, a: usize, b: usize) -> Point3 {
    match (a, b) {
        (0, 0) => g.h[0][0],
        (1, 0) => linear(g, tau.d[1][0]),
        (0, 1) => linear(g, tau.d[0][1]),
        (2, 0) => second(g, tau.d[1][0], tau.d[1][0], tau.d[2][0]),
        (1, 1) => second(g, tau.d[1][0], tau.d[0][1], tau.d[1][1]),
        (0, 2) => second(g, tau.d[0][1], tau.d[0][1], tau.d[0][2]),
        (3, 0) => third(
            g,
            tau.d[1][0],
            tau.d[1][0],
            tau.d[1][0],
            tau.d[2][0],
            tau.d[2][0],
            tau.d[2][0],
            tau.d[3][0],
        ),
        (2, 1) => third(
            g,
            tau.d[1][0],
            tau.d[1][0],
            tau.d[0][1],
            tau.d[2][0],
            tau.d[1][1],
            tau.d[1][1],
            tau.d[2][1],
        ),
        (1, 2) => third(
            g,
            tau.d[1][0],
            tau.d[0][1],
            tau.d[0][1],
            tau.d[1][1],
            tau.d[1][1],
            tau.d[0][2],
            tau.d[1][2],
        ),
        (0, 3) => third(
            g,
            tau.d[0][1],
            tau.d[0][1],
            tau.d[0][1],
            tau.d[0][2],
            tau.d[0][2],
            tau.d[0][2],
            tau.d[0][3],
        ),
        _ => Point3::new(0.0, 0.0, 0.0),
    }
}

fn linear(g: Kp20GuideJet3, a: [f64; 2]) -> Point3 {
    point_sum(&[(a[0], g.h[1][0]), (a[1], g.h[0][1])])
}

fn second(g: Kp20GuideJet3, a: [f64; 2], b: [f64; 2], ab: [f64; 2]) -> Point3 {
    point_sum(&[
        (a[0] * b[0], g.h[2][0]),
        (a[0] * b[1] + a[1] * b[0], g.h[1][1]),
        (a[1] * b[1], g.h[0][2]),
        (ab[0], g.h[1][0]),
        (ab[1], g.h[0][1]),
    ])
}

#[allow(clippy::too_many_arguments)]
fn third(
    g: Kp20GuideJet3,
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    ab: [f64; 2],
    ac: [f64; 2],
    bc: [f64; 2],
    abc: [f64; 2],
) -> Point3 {
    let third_order = [
        (a[0] * b[0] * c[0], g.h[3][0]),
        (
            a[0] * b[0] * c[1] + a[0] * b[1] * c[0] + a[1] * b[0] * c[0],
            g.h[2][1],
        ),
        (
            a[0] * b[1] * c[1] + a[1] * b[0] * c[1] + a[1] * b[1] * c[0],
            g.h[1][2],
        ),
        (a[1] * b[1] * c[1], g.h[0][3]),
    ];
    let second_order = [
        (ab[0] * c[0] + ac[0] * b[0] + bc[0] * a[0], g.h[2][0]),
        (
            ab[0] * c[1] + ac[0] * b[1] + bc[0] * a[1] + ab[1] * c[0] + ac[1] * b[0] + bc[1] * a[0],
            g.h[1][1],
        ),
        (ab[1] * c[1] + ac[1] * b[1] + bc[1] * a[1], g.h[0][2]),
    ];
    let first_order = [(abc[0], g.h[1][0]), (abc[1], g.h[0][1])];
    point_sum(&[
        third_order[0],
        third_order[1],
        third_order[2],
        third_order[3],
        second_order[0],
        second_order[1],
        second_order[2],
        first_order[0],
        first_order[1],
    ])
}

fn point_sum(terms: &[(f64, Point3)]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for &(w, p) in terms {
        out.x += w * p.x;
        out.y += w * p.y;
        out.z += w * p.z;
    }
    out
}

fn derivative_basis(n: usize, du: usize, dv: usize, u: f64, v: f64) -> Vec<f64> {
    if du == 0 && dv == 0 {
        return eval_basis(n, u, v);
    }
    debug_assert!(du + dv <= n);
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
const fn flat_index(n: usize, i: usize, j: usize) -> usize {
    i * (2 * n + 3 - i) / 2 + j
}

#[inline]
const fn tri_cp_count(n: usize) -> usize {
    (n + 1) * (n + 2) / 2
}

#[inline]
fn multinomial(n: usize, i: usize, j: usize, k: usize) -> f64 {
    factorial(n) / (factorial(i) * factorial(j) * factorial(k))
}

#[inline]
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, v| acc * v as f64)
}

#[inline]
fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn scale2(a: [f64; 2], s: f64) -> [f64; 2] {
    [a[0] * s, a[1] * s]
}

#[cfg(test)]
#[path = "kp20_bi6_tau_tests.rs"]
mod tests;
