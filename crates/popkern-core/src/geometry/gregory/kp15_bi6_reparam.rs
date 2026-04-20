//! Phase I.5 — Karčiauskas-Peters 2015 bi-6 reparameterisations
//! `ρ` and `β`.
//!
//! Shape spec:
//! `docs/journal/polar-kp-math.md` §2.1 (ρ), §2.4 (β + γ table).
//!
//! # Between-sector reparameterisation `ρ`
//!
//! Couples adjacent sectors via a 2D map whose Taylor expansion in
//! `v` up to second order is pinned by G² requirements across the
//! shared edge. From [KP15 §3.1 eq (4)]:
//!
//! ```text
//!   ρ(u, v) = ( u + b(u)·v + ½·e(u)·v²,   a(u)·v + ½·d(u)·v² )
//!   a(u) := −1
//!   b(u) := 2c · (1 − u)²
//!   d(u) := δ · (1 − u)² / (1 − u + τ·u),     δ := 24c(τ − 1)
//!   e(u) := b(u) · (b'(u) − d(u) / 2)
//! ```
//!
//! where `c := cos(2π/n)` and `τ` is the per-valence parameter
//! from [KP15 §3.6 step 2]:
//!
//! | n | τ_n |
//! |---|-----|
//! | 3 | 4/3 |
//! | 5 | 0.87 |
//! | 6 | 0.85 |
//! | 7 | 0.84 |
//! | 8 | 0.83 |
//! | 9 | 0.82 |
//! | 10 | 0.81 |
//! | >10 | 0.80 |
//!
//! # Tensor-border reparameterisation `β`
//!
//! Preserves position and first cross-derivative (`a_β(u) := 1,
//! b_β(u) := 0`) and reshapes only the second cross-derivative.
//! The shape scalars `d(u)` (cubic in `u`) and `e(u)` (quartic
//! in `u`) are specified as BB-coefficients [KP15 eq (12)]:
//!
//! ```text
//!   d₀ := 0,  d₁ := 0,  d₃ := 2c/3 + d₂
//!   e₀ := 0,  e₁ := 0,  e₂ := d₂/2,
//!   e₃ := (4c(1 − τ) + τ·d₂) / (4τ),   e₄ := 0
//! ```
//!
//! `d₂` is pinned by [KP15 eq (13)] for `n ≠ 3`:
//!
//! ```text
//!   d₂ := (1 − (3 + 2c)τ + (3 + c)τ²) / (2τ²)
//! ```
//!
//! For `n = 3`, `d₂ := 0.1` (paper convention, §3.5).
//!
//! # γ weights for tensor-border reparameterisation
//!
//! [KP15 eq (17)] provides the linear map `q_{k,2} = Σ γ^k_{ij} ·
//! b_{ij}` for `k ∈ 2..=6`, combining the depth-2 bi-3 tensor-
//! border `b` into the reparameterised row-2 CP at position `k`.
//! The full `γ^k_{ij}` table is tabulated inside [`gamma`].

use std::f64::consts::PI;

/// Per-valence and per-input parameters used throughout the K-P
/// 2015 bi-6 construction.
///
/// Construct via [`Bi6Reparams::for_valence`]; individual fields
/// can be overridden for experiments (e.g., forcing a non-default
/// `d_2` at `n = 3`).
#[derive(Debug, Clone, Copy)]
pub struct Bi6Reparams {
    /// Valence `n`.
    pub n: usize,
    /// `c := cos(2π / n)`.
    pub c: f64,
    /// `τ` per [KP15 §3.6 step 2]'s table.
    pub tau: f64,
    /// `d₂` per [KP15 eq (13)] (n ≠ 3) or the `0.1` default
    /// (n = 3).
    pub d2: f64,
}

impl Bi6Reparams {
    /// Default parameters for valence `n`. `τ` from the published
    /// table, `d₂` from (13) for `n ≠ 3` or the `0.1` paper
    /// default for `n = 3`.
    ///
    /// Panics if `n < 3`.
    pub fn for_valence(n: usize) -> Self {
        assert!(n >= 3, "K-P bi-6 needs valence ≥ 3; got {n}");
        let c = (2.0 * PI / n as f64).cos();
        let tau = tau_for_valence(n);
        let d2 = if n == 3 {
            0.1
        } else {
            // d₂ = (1 − (3 + 2c)τ + (3 + c)τ²) / (2τ²)
            (1.0 - (3.0 + 2.0 * c) * tau + (3.0 + c) * tau * tau) / (2.0 * tau * tau)
        };
        Self { n, c, tau, d2 }
    }

    /// `a_ρ(u) = −1` (constant).
    #[inline]
    pub fn rho_a(&self, _u: f64) -> f64 {
        -1.0
    }

    /// `b_ρ(u) = 2c · (1 − u)²`.
    #[inline]
    pub fn rho_b(&self, u: f64) -> f64 {
        let omu = 1.0 - u;
        2.0 * self.c * omu * omu
    }

    /// `b_ρ'(u) = −4c · (1 − u)`.
    #[inline]
    pub fn rho_b_prime(&self, u: f64) -> f64 {
        -4.0 * self.c * (1.0 - u)
    }

    /// `d_ρ(u) = δ · (1 − u)² / (1 − u + τ·u)` with
    /// `δ := 24c(τ − 1)`. Denominator is `1 − u + τu` — equals
    /// `1` at `u = 0`, equals `τ` at `u = 1`, never vanishes on
    /// `[0, 1]` for `τ > 0`.
    #[inline]
    pub fn rho_d(&self, u: f64) -> f64 {
        let delta = 24.0 * self.c * (self.tau - 1.0);
        let omu = 1.0 - u;
        let denom = omu + self.tau * u;
        delta * omu * omu / denom
    }

    /// `e_ρ(u) = b(u) · (b'(u) − d(u) / 2)`.
    #[inline]
    pub fn rho_e(&self, u: f64) -> f64 {
        self.rho_b(u) * (self.rho_b_prime(u) - self.rho_d(u) / 2.0)
    }

    /// `β`'s `d(u)` BB-CPs `[d₀, d₁, d₂, d₃]` (cubic in `u`).
    /// `d₀ = d₁ = 0`; `d₃ = 2c/3 + d₂`.
    pub fn beta_d_bb_coeffs(&self) -> [f64; 4] {
        [0.0, 0.0, self.d2, 2.0 * self.c / 3.0 + self.d2]
    }

    /// `β`'s `e(u)` BB-CPs `[e₀, e₁, e₂, e₃, e₄]` (quartic in `u`).
    /// `e₀ = e₁ = e₄ = 0`; `e₂ = d₂/2`;
    /// `e₃ = (4c(1 − τ) + τ·d₂) / (4τ)`.
    pub fn beta_e_bb_coeffs(&self) -> [f64; 5] {
        let e3 = (4.0 * self.c * (1.0 - self.tau) + self.tau * self.d2) / (4.0 * self.tau);
        [0.0, 0.0, self.d2 / 2.0, e3, 0.0]
    }
}

/// `τ` value for valence `n` per [KP15 §3.6 step 2]'s table.
///
/// Panics if `n < 3` (K-P bi-6 is not defined).
pub fn tau_for_valence(n: usize) -> f64 {
    assert!(n >= 3, "tau_for_valence needs n ≥ 3; got {n}");
    match n {
        3 => 4.0 / 3.0,
        4 => 4.0 / 3.0, // Not in the published table; Phase H
        // handles N = 4 rectangularly so this is
        // only exercised in testing / niche dispatches.
        5 => 0.87,
        6 => 0.85,
        7 => 0.84,
        8 => 0.83,
        9 => 0.82,
        10 => 0.81,
        _ => 0.80,
    }
}

/// `γ^k_{ij}` weight for `k ∈ 2..=6`, `(i, j)` within the valid
/// index set for that `k`. Zero for any `(k, i, j)` outside the
/// sparsity pattern documented in the module source.
///
/// Paper reference: [KP15 Appendix "γ table"] (math summary §2.4).
///
/// Indexing convention: `i` ∈ 0..=3 for a bi-3 tensor-border's
/// position index along the sector edge; `j` ∈ 0..=2 for the
/// depth (`j=0` is the position row, `j=1` the first cross-
/// derivative, `j=2` the second cross-derivative).
pub fn gamma(reparams: &Bi6Reparams, k: usize, i: usize, j: usize) -> f64 {
    let c = reparams.c;
    let tau = reparams.tau;
    let d2 = reparams.d2;
    match (k, i, j) {
        // k = 2
        (2, 0, 0) => 1.0 - d2 / 25.0,
        (2, 1, 0) => (6.0 + d2) / 50.0,
        (2, 2, 0) => 1.0 / 25.0,
        (2, 1, 1) => 9.0 / 25.0,
        (2, 2, 1) => 3.0 / 25.0,
        (2, 2, 2) => 1.0 / 25.0,

        // k = 3
        (3, 0, 0) => ((3.0 + 5.0 * c - 3.0 * d2) * tau - 6.0 * c) / (300.0 * tau),
        (3, 1, 0) => ((9.0 - 2.0 * c - 7.0 * d2) * tau + 2.0 * c) / (100.0 * tau),
        (3, 2, 0) => (9.0 + 3.0 * d2) / 100.0,
        (3, 3, 0) => 1.0 / 100.0,
        (3, 0, 1) => (18.0 + 2.0 * c + 3.0 * d2) / 600.0,
        (3, 1, 1) => (54.0 + 9.0 * d2) / 200.0,
        (3, 2, 1) => 27.0 / 100.0,
        (3, 3, 1) => 3.0 / 100.0,
        (3, 0, 2) => 1.0 / 100.0,
        (3, 1, 2) => 9.0 / 100.0,
        (3, 2, 2) => 9.0 / 100.0,
        (3, 3, 2) => 1.0 / 100.0,

        // k = 4
        (4, 1, 0) => ((6.0 + 6.0 * c - 5.0 * d2) * tau - 8.0 * c) / (150.0 * tau),
        (4, 2, 0) => ((9.0 - 4.0 * c - 5.0 * d2) * tau + 4.0 * c) / (75.0 * tau),
        (4, 3, 0) => (2.0 + d2) / 50.0,
        (4, 1, 1) => (18.0 + 2.0 * c + 3.0 * d2) / 150.0,
        (4, 2, 1) => (18.0 + 3.0 * d2) / 50.0,
        (4, 3, 1) => 3.0 / 25.0,
        (4, 1, 2) => 1.0 / 25.0,
        (4, 2, 2) => 3.0 / 25.0,
        (4, 3, 2) => 1.0 / 25.0,

        // k = 5
        (5, 2, 0) => ((3.0 + c - 2.0 * d2) * tau - 2.0 * c) / (30.0 * tau),
        (5, 3, 0) => ((3.0 - 2.0 * c - d2) * tau + 2.0 * c) / (30.0 * tau),
        (5, 2, 1) => (18.0 + 2.0 * c + 3.0 * d2) / 60.0,
        (5, 3, 1) => (6.0 + d2) / 20.0,
        (5, 2, 2) => 1.0 / 10.0,
        (5, 3, 2) => 1.0 / 10.0,

        // k = 6
        (6, 3, 0) => (6.0 - 2.0 * c - 3.0 * d2) / 30.0,
        (6, 3, 1) => (18.0 + 2.0 * c + 3.0 * d2) / 30.0,
        (6, 3, 2) => 1.0 / 5.0,

        _ => 0.0,
    }
}

#[cfg(test)]
#[path = "kp15_bi6_reparam_tests.rs"]
mod tests;
