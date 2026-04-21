//! Phase I.13 — Karčiauskas-Peters 2013 bi-5 reparameterisations
//! `ρ` and `β`.
//!
//! Paper: [KP13] Karčiauskas & Peters (2013), "Biquintic G² surfaces".
//! See `docs/papers/kk12bi5.pdf` §4.
//!
//! # Between-sector reparameterisation `ρ`
//!
//! For valence `n > 3` (n = 3 handled separately in §6), symmetry
//! pins the generic Taylor expansion [KP13 eq (9)] to the form
//! [KP13 eq (10)]:
//!
//! ```text
//!   ρ = ( u + b(u)·[v + ½(b'(u) - d(u)/2)·v²],   -v + ½·d(u)·v² )
//! ```
//!
//! so `a := -1` and `e(u) := b(u)·(b'(u) - d(u)/2)`.
//!
//! Two pieces of `ρ` cover one sector:
//!
//! * **`ρ_o`** (top of sector — near the star point `o`):
//!   `b(u) := 2c(1-u) + (2/3)c·u`, `d(u) := 0`.
//!   Paper: see paragraph following eq (10).
//! * **`ρ_ℓr`** (bottom of sector — near the mid-edge `m_k`):
//!   `b(u) := (2/3)c(1-u)²`, `d(u) := 0`.
//!
//! where `c := cos(2π/n)`.
//!
//! ρ_o and ρ_ℓr are not C²-connected, but the corresponding patches
//! `p^o` and `p^ℓ` / `p^r` are made C²-connected via the cap-assembly
//! equations eq (11)-(13).
//!
//! # Tensor-border reparameterisation `β`
//!
//! `β` is a two-piece bi-3 map that rereshapes the bi-3 input
//! tensor-border. Identity up to and including the first cross-
//! derivative (`a := 1`, `b(u) := 0` in eq (14)), so only the
//! second cross-derivative is reshaped. The map has two pieces
//! `β^m` (middle, star-near) and `β^r` (right), joined internally
//! with C² continuity; `β^m` is symmetric across its diagonal.
//!
//! [KP13 eq (15)] fixes the free BB-coefficients:
//!
//! ```text
//!   β^m_{22} := (1/3 + c/216,    1/3 + c/216)
//!   β^m_{32} := (1/2 + 7c/1296,  1/3 + c/72)
//!   β^r_{32} := (1,              1/3 + c/18)
//! ```
//!
//! These are data-independent and can be tabulated per-n.
//!
//! Several accessors below are not yet consumed — they land for
//! I.14 (cap assembly) which will call the eq (11)-(13) inter-
//! sector join machinery.

#![allow(dead_code)]

use std::f64::consts::PI;

/// Per-valence parameters for the KP13 bi-5 construction.
///
/// Unlike the bi-6 analogue, bi-5 has no `τ` or `d₂` free tuning
/// parameter — the paper pins all coefficients in closed form via
/// `c := cos(2π/n)`.
#[derive(Debug, Clone, Copy)]
pub struct Bi5Reparams {
    /// Valence `n`.
    pub n: usize,
    /// `c := cos(2π/n)`.
    pub c: f64,
}

impl Bi5Reparams {
    /// Default parameters for valence `n`.
    ///
    /// Panics if `n < 3`. (The n = 3 case uses the special-case
    /// construction in `kp13_bi5_cap_n3` but still needs `c`.)
    pub fn for_valence(n: usize) -> Self {
        assert!(n >= 3, "KP13 bi-5 needs valence ≥ 3; got {n}");
        let c = (2.0 * PI / n as f64).cos();
        Self { n, c }
    }

    // ---- ρ (between-sector) ----

    /// `a_ρ(u) = −1` (constant), from symmetry (KP13 §4.1).
    #[inline]
    pub fn rho_a(&self, _u: f64) -> f64 {
        -1.0
    }

    // ρ_o: top of sector (near star point `o`).

    /// `ρ_o`'s `b(u) = 2c(1-u) + (2/3)c·u`.
    #[inline]
    pub fn rho_o_b(&self, u: f64) -> f64 {
        2.0 * self.c * (1.0 - u) + (2.0 / 3.0) * self.c * u
    }

    /// `ρ_o`'s `b'(u) = -4c/3` (constant).
    #[inline]
    pub fn rho_o_b_prime(&self, _u: f64) -> f64 {
        -4.0 * self.c / 3.0
    }

    /// `ρ_o`'s `d(u) = 0` (KP13 §4.1 choice).
    #[inline]
    pub fn rho_o_d(&self, _u: f64) -> f64 {
        0.0
    }

    /// `ρ_o`'s `e(u) = b(u)·b'(u)` (since d = 0).
    #[inline]
    pub fn rho_o_e(&self, u: f64) -> f64 {
        self.rho_o_b(u) * self.rho_o_b_prime(u)
    }

    // ρ_ℓr: bottom of sector (near mid-edge `m_k`).

    /// `ρ_ℓr`'s `b(u) = (2/3)c(1-u)²`.
    #[inline]
    pub fn rho_lr_b(&self, u: f64) -> f64 {
        let omu = 1.0 - u;
        (2.0 / 3.0) * self.c * omu * omu
    }

    /// `ρ_ℓr`'s `b'(u) = -(4c/3)(1-u)`.
    #[inline]
    pub fn rho_lr_b_prime(&self, u: f64) -> f64 {
        -(4.0 / 3.0) * self.c * (1.0 - u)
    }

    /// `ρ_ℓr`'s `d(u) = 0`.
    #[inline]
    pub fn rho_lr_d(&self, _u: f64) -> f64 {
        0.0
    }

    /// `ρ_ℓr`'s `e(u) = b(u)·b'(u)` (since d = 0).
    #[inline]
    pub fn rho_lr_e(&self, u: f64) -> f64 {
        self.rho_lr_b(u) * self.rho_lr_b_prime(u)
    }

    // ---- β (tensor-border reparameterisation) ----

    /// `β^m_{22} = (1/3 + c/216, 1/3 + c/216)` (diagonal by
    /// symmetry). Paper: [KP13 eq (15)].
    pub fn beta_m_22(&self) -> (f64, f64) {
        let x = 1.0 / 3.0 + self.c / 216.0;
        (x, x)
    }

    /// `β^m_{32} = (1/2 + 7c/1296, 1/3 + c/72)`. Paper: [KP13
    /// eq (15)].
    pub fn beta_m_32(&self) -> (f64, f64) {
        (0.5 + 7.0 * self.c / 1296.0, 1.0 / 3.0 + self.c / 72.0)
    }

    /// `β^r_{32} = (1, 1/3 + c/18)`. Paper: [KP13 eq (15)].
    pub fn beta_r_32(&self) -> (f64, f64) {
        (1.0, 1.0 / 3.0 + self.c / 18.0)
    }

    // ---- inter-sector join weights (eq 12) ----

    /// `w_1 = -c/96` (KP13 eq 12).
    #[inline]
    pub fn w1(&self) -> f64 {
        -self.c / 96.0
    }

    /// `w_2 = c/96` (KP13 eq 12).
    #[inline]
    pub fn w2(&self) -> f64 {
        self.c / 96.0
    }

    /// `w_3 = c/30` (KP13 eq 12).
    #[inline]
    pub fn w3(&self) -> f64 {
        self.c / 30.0
    }

    /// Structural sanity check — verify fields are finite and
    /// `n ≥ 3`. [`for_valence`](Self::for_valence) output always
    /// passes.
    pub fn validate(&self) -> Result<(), String> {
        if self.n < 3 {
            return Err(format!("KP13 bi-5 valence n = {} < 3", self.n));
        }
        if !self.c.is_finite() {
            return Err(format!("KP13 bi-5 c = {} is not finite", self.c));
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "kp13_bi5_reparam_tests.rs"]
mod tests;
