//! Phase I.11.d.2 — Characteristic map `χ` for the KP13 bi-5 LSQ.
//!
//! KP13 §3 eq (4) minimises `Σ_k Σ_{t ∈ {m_k, v_k}} ‖H_t(g△^k ∘ χ^k)
//! − H_t(b^k)‖²` over free BB-coefs of `g△`. Evaluating the
//! composite `g△^k ∘ χ^k` at the corners of the unit square requires
//! the characteristic map `χ^k : [0, 1]² → Ω`.
//!
//! # Paper specification (KP13 §3)
//!
//! > `χ : [0..1]^{2×n} → R², normalized to edge-length 1 [...] so
//! > that `g△^k ∘ χ^k(1, 1) = v_k` (`χ^k(0, 0) = (0, 0)`, Fig. 5a)
//! > and `g△^k_{050} = m_k` (Fig. 5b).
//!
//! Standard reading: `χ` is the planar Catmull-Clark characteristic
//! map restricted to the n extraordinary-vertex sectors, normalised
//! so the first ring has edge length 1. Unit-square corners map to:
//!
//! | corner `(u, v)` | `χ^k(u, v)` in Ω |
//! |-----------------|-------------------|
//! | `(0, 0)`        | `o = (0, 0)` (star) |
//! | `(1, 0)`        | `m_k` (= `(c_k, s_k)`) |
//! | `(0, 1)`        | `m_{k+1}` (= `(c_{k+1}, s_{k+1})`) |
//! | `(1, 1)`        | `v_k` (outside the natural n-gon Ω) |
//!
//! # Current status — bilinear approximation
//!
//! This module ships a **bilinear** `χ^k(u, v) := u·m_k + v·m_{k+1}`
//! as a working approximation. This satisfies the first three corner
//! conditions exactly, and interprets `v_k := m_k + m_{k+1}` for the
//! last corner. Partial derivatives are constant (bilinear has no
//! curvature).
//!
//! The **true** CC characteristic map differs at interior points
//! (via CC subdivision eigenanalysis), and sets `v_k` to a CC-
//! specific position outside Ω. For the LSQ, bilinear gives a
//! well-defined boundary-fit objective that exercises the full
//! machinery; swapping in the proper CC map is additive work that
//! can land in I.11e or I.12 once boundary-fit quality is measured.
//!
//! See `docs/journal/polar-kp-math.md` §1.5 for the CC limit-point
//! formula (already implemented), and KP13 Fig 6b / §4.3 σ for the
//! related "characteristic parameterization" that σ approximates.

#![allow(dead_code)]

use super::kp13_bi5_gtri5::sector_corner;

/// Value of `χ^k` at a corner of the unit square, plus its first
/// and second partial derivatives. All quantities are in the
/// 2-D domain Ω (R²).
#[derive(Debug, Clone, Copy)]
pub struct ChiJet {
    /// `χ^k(u, v)` evaluated at the corner.
    pub pos: (f64, f64),
    /// `∂χ^k / ∂u` at the corner.
    pub d_u: (f64, f64),
    /// `∂χ^k / ∂v` at the corner.
    pub d_v: (f64, f64),
    /// `∂²χ^k / ∂u²` at the corner.
    pub d_uu: (f64, f64),
    /// `∂²χ^k / ∂v²` at the corner.
    pub d_vv: (f64, f64),
    /// `∂²χ^k / ∂u∂v` at the corner.
    pub d_uv: (f64, f64),
}

/// Bilinear approximation `χ^k(u, v) := u · m_k + v · m_{k+1}` to
/// the Catmull-Clark characteristic map.
///
/// Per `m_k = (cos(k · 2π/n), sin(k · 2π/n))`.
#[derive(Debug, Clone, Copy)]
pub struct BilinearChi {
    pub n: usize,
}

/// Which corner of the unit square to query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChiCorner {
    /// `(0, 0)` — star point.
    Star,
    /// `(1, 0)` — `m_k`.
    Mk,
    /// `(0, 1)` — `m_{k+1}`.
    Mk1,
    /// `(1, 1)` — `v_k`.
    Vk,
}

impl ChiCorner {
    #[inline]
    pub fn uv(self) -> (f64, f64) {
        match self {
            ChiCorner::Star => (0.0, 0.0),
            ChiCorner::Mk => (1.0, 0.0),
            ChiCorner::Mk1 => (0.0, 1.0),
            ChiCorner::Vk => (1.0, 1.0),
        }
    }
}

impl BilinearChi {
    /// New bilinear χ for valence `n ≥ 3`.
    pub fn new(n: usize) -> Self {
        assert!(n >= 3, "bilinear χ needs n ≥ 3, got {n}");
        Self { n }
    }

    /// Evaluate `χ^k(u, v) = u · m_k + v · m_{k+1}`.
    #[inline]
    pub fn eval(&self, k: usize, u: f64, v: f64) -> (f64, f64) {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        (u * mkx + v * mk1x, u * mky + v * mk1y)
    }

    /// Full `ChiJet` at the given corner of sector `k`. For the
    /// bilinear map, second partials are zero and first partials
    /// are constant.
    pub fn jet(&self, k: usize, corner: ChiCorner) -> ChiJet {
        let (mkx, mky) = sector_corner(self.n, k);
        let (mk1x, mk1y) = sector_corner(self.n, k + 1);
        let (u, v) = corner.uv();
        ChiJet {
            pos: (u * mkx + v * mk1x, u * mky + v * mk1y),
            d_u: (mkx, mky),
            d_v: (mk1x, mk1y),
            d_uu: (0.0, 0.0),
            d_vv: (0.0, 0.0),
            d_uv: (0.0, 0.0),
        }
    }
}

#[cfg(test)]
#[path = "kp13_bi5_chi_tests.rs"]
mod tests;
